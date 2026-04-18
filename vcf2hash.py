#!/usr/bin/env python3
import sys, gzip, argparse, hashlib, json, re

# --- I/O helpers ---
def open_maybe_gz(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

# --- DNA utils ---
RC = str.maketrans('ACGTWSMKRYNacgtwsmkryn', 'TGCAWSKMYRNtgcawskmyrn')
def revcomp(seq: str) -> str:
    return seq.translate(RC)[::-1]

def sha256_hex(s: str) -> str:
    return hashlib.sha256(s.encode('utf-8')).hexdigest()

# --- Chromosome normalization (matches your original approach) ---
def norm_chrom(s: str) -> str:
    s = str(s)
    s = re.sub(r'^glyma\.Wm82\.gnm2\.', '', s, flags=re.I)  # drop prefix if present
    s = re.sub(r'^(chromosome|chr)', 'Chr', s, flags=re.I)  # unify 'chr' casing
    # Zero-folding: Chr01 -> Chr1
    m = re.match(r'^Chr0*([1-9]\d*)$', s)
    if m: return f'Chr{m.group(1)}'
    return s

# --- FASTA ---
def read_fasta(path):
    seqs = {}
    name, buf = None, []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = ''.join(buf)
    return seqs

# --- VCF parsing: build per-chrom SNP edits (0-based positions, allele codes incl. IUPAC for het) ---
IUPAC = {
    'AA':'A','TT':'T','CC':'C','GG':'G',
    'AC':'M','CA':'M','GT':'K','TG':'K',
    'AG':'R','GA':'R','CT':'Y','TC':'Y',
    'AT':'W','TA':'W','CG':'S','GC':'S'
}
def encode_genotype(ref, alts, gt):
    if not gt or gt in {'.', './.', '.|.'} or gt[0] == '.':
        return 'N'
    parts = gt.replace('|','/').split('/')
    try:
        alleles = [ ([ref] + alts)[int(i)] for i in parts ]
    except Exception:
        return 'N'
    if any(len(a) != 1 for a in alleles):
        return 'N'
    key = ''.join(sorted(alleles))
    return IUPAC.get(key, 'N')

def vcf_snps(vcf_path):
    edits = {}  # chrom_norm -> list of (pos0, base)
    with open_maybe_gz(vcf_path) as f:
        sample_col = None
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                cols = line.strip().split('\t')
                sample_col = 9 if len(cols) > 9 else None
                continue
            if not line or line[0] == '#':
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8: continue
            chrom, pos, _id, ref, alt, qual, filt, info = cols[:8]
            # SNP only, skip INDELs
            if re.search(r'(^|;)INDEL(=|;|$)', info):
                continue
            alts = alt.split(',')
            if len(ref) != 1 or any(len(a) != 1 for a in alts):
                continue
            gt = None
            if sample_col is not None and len(cols) > sample_col:
                gt = cols[sample_col].split(':', 1)[0]
            allele = encode_genotype(ref, alts, gt)
            pos0 = int(pos) - 1
            k = norm_chrom(chrom)
            edits.setdefault(k, []).append((pos0, allele))
    # keep last edit if duplicates per position
    for k, lst in edits.items():
        d = {}
        for p, a in lst:
            d[p] = a
        edits[k] = sorted(d.items())  # list of (pos0, allele) sorted by pos
    return edits

# --- GFF3 parsing: collect genes (chrom, start, end, strand, gene_id) ---
def parse_attrs(attr_str: str):
    out = {}
    for item in attr_str.strip().strip(';').split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            out[k] = v
    return out

def gff_genes(gff_path):
    genes = []  # dicts with chrom_norm, start, end, strand, id
    with open_maybe_gz(gff_path) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9: continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            if ftype != 'gene':
                continue
            a = parse_attrs(attrs)
            gid = a.get('ID', a.get('Name', f"{seqid}:{start}-{end}"))
            genes.append({
                'chrom_norm': norm_chrom(seqid),
                'start': int(start),
                'end': int(end),
                'strand': strand,
                'gene_id': gid
            })
    return genes

# --- Apply SNP edits to FASTA sequences (in-place copy) ---
def apply_edits_to_genome(fasta_dict):
    # Build norm -> fasta key map
    norm2fa = {}
    for k in fasta_dict.keys():
        norm2fa[norm_chrom(k)] = k
    return norm2fa

def modify_chrom(seq: str, edits):
    if not edits: return seq
    arr = list(seq)
    n = len(arr)
    for pos0, base in edits:
        if 0 <= pos0 < n:
            arr[pos0] = base
    return ''.join(arr)

# --- Extract gene sequences from modified genome, strand-aware, and hash them ---
def gene_hashes(fasta_dict, edits_by_chr, genes):
    norm2fa = apply_edits_to_genome(fasta_dict)
    # precompute modified sequences lazily per chrom
    modified_cache = {}
    out = []
    for g in genes:
        cn = g['chrom_norm']
        if cn not in norm2fa:
            # try a fallback: maybe FASTA already uses the norm as key
            if cn in fasta_dict:
                fa_key = cn
            else:
                continue
        else:
            fa_key = norm2fa[cn]
        # build modified chrom once
        if fa_key not in modified_cache:
            edits = edits_by_chr.get(cn, [])
            modified_cache[fa_key] = modify_chrom(fasta_dict[fa_key], edits)
        chrom_seq = modified_cache[fa_key]
        s0 = g['start'] - 1
        e0 = g['end']   # slice end is exclusive
        if s0 < 0 or e0 > len(chrom_seq) or s0 >= e0:
            continue
        seq = chrom_seq[s0:e0]
        if g['strand'] == '-':
            seq = revcomp(seq)
        h = sha256_hex(seq)
        out.append({'gene_id': g['gene_id'], 'hash': h})
    return out

# --- CLI ---
def main():
    ap = argparse.ArgumentParser(description="Hash SNP-reflected gene sequences using FASTA+GFF3+VCF.")
    ap.add_argument('--fasta', required=True, help='Reference FASTA')
    ap.add_argument('--gff3', required=True, help='GFF3 annotations (genes)')
    ap.add_argument('--vcf',  required=True, help='Sample VCF/VCF.gz (single-sample preferred)')
    ap.add_argument('--out', default='-', help='Output JSONL (gene_id\\thash if --tsv)')
    ap.add_argument('--tsv', action='store_true', help='Write TSV (gene_id\\thash) instead of JSONL')
    args = ap.parse_args()

    fasta = read_fasta(args.fasta)
    edits = vcf_snps(args.vcf)
    genes = gff_genes(args.gff3)
    gh = gene_hashes(fasta, edits, genes)

    if args.tsv:
        out = sys.stdout if args.out == '-' else open(args.out, 'w')
        with out:
            for r in gh:
                out.write(f"{r['gene_id']}\t{r['hash']}\n")
    else:
        out = sys.stdout if args.out == '-' else open(args.out, 'w')
        with out:
            for r in gh:
                out.write(json.dumps(r, separators=(',', ':')) + "\n")

if __name__ == '__main__':
    main()
