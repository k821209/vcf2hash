#!/usr/bin/env python3
"""
vcf2hash: Gene-region allele profiling via SHA-256 hashing.

Converts a multi-sample VCF file into gene-region haplotypes by encoding
each gene's SNP constellation as an IUPAC string and assigning a
deterministic haplotype identifier via SHA-256 hashing.

Usage:
    python vcf2hash.py --fasta ref.fa --gff3 genes.gff3 --vcf sample.vcf.gz --out sample.tsv --tsv
"""
import sys, gzip, argparse, hashlib, json, re
from collections import defaultdict

# --- I/O helpers ---
def open_maybe_gz(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

# --- DNA utils ---
RC = str.maketrans('ACGTWSMKRYNacgtwsmkryn', 'TGCAWSKMYRNtgcawskmyrn')
def revcomp(seq: str) -> str:
    return seq.translate(RC)[::-1]

def sha256_hex(s: str) -> str:
    return hashlib.sha256(s.encode('utf-8')).hexdigest()

# --- Chromosome name matching ---
def _canon(s: str) -> str:
    """Reduce a chromosome name to a canonical form for matching.
    Strips common prefixes (chr, chromosome, scaffold, species.assembly.version.*)
    and normalizes case + leading zeros.
    """
    s = str(s).strip()
    # Strip dotted prefixes (e.g. glyma.Wm82.gnm2.Chr01 -> Chr01)
    # Pattern: word.word.word.REMAINDER — keep REMAINDER
    s = re.sub(r'^(\w+\.){2,}', '', s)
    # Strip chr/chromosome/scaffold prefix
    s = re.sub(r'^(chromosome|scaffold|chr)', '', s, flags=re.I)
    # Remove leading zeros from numeric part (01 -> 1)
    s = re.sub(r'^0+(\d)', r'\1', s)
    return s.lower()

def build_chrom_map(fasta_keys, gff_chroms, vcf_chroms):
    """Build a unified chromosome name mapping.
    Returns a dict: any raw name -> canonical key, and
    a dict: canonical key -> fasta key.
    """
    # canonical -> fasta key
    canon2fasta = {}
    for k in fasta_keys:
        c = _canon(k)
        canon2fasta[c] = k
        # Also store exact match
        canon2fasta[k] = k

    # Build raw -> canonical mapping for all sources
    raw2canon = {}
    for k in list(fasta_keys) + list(gff_chroms) + list(vcf_chroms):
        raw2canon[k] = _canon(k)

    return raw2canon, canon2fasta

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

# --- VCF parsing ---
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
        alleles = [([ref] + alts)[int(i)] for i in parts]
    except Exception:
        return 'N'
    if any(len(a) != 1 for a in alleles):
        return 'N'
    key = ''.join(sorted(alleles))
    return IUPAC.get(key, 'N')

def vcf_snps(vcf_path):
    """Parse VCF, return dict: raw_chrom -> list of (pos0, allele_base)."""
    edits = defaultdict(list)
    chroms_seen = set()
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
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info = cols[:8]
            chroms_seen.add(chrom)
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
            edits[chrom].append((pos0, allele))
    # Deduplicate: keep last edit per position
    deduped = {}
    for k, lst in edits.items():
        d = {}
        for p, a in lst:
            d[p] = a
        deduped[k] = sorted(d.items())
    return deduped, chroms_seen

# --- GFF3 parsing ---
def parse_attrs(attr_str: str):
    out = {}
    for item in attr_str.strip().strip(';').split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            out[k] = v
    return out

def gff_genes(gff_path):
    """Parse GFF3, return list of gene dicts and set of chroms seen."""
    genes = []
    chroms_seen = set()
    with open_maybe_gz(gff_path) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            if ftype != 'gene':
                continue
            chroms_seen.add(seqid)
            a = parse_attrs(attrs)
            gid = a.get('ID', a.get('Name', f"{seqid}:{start}-{end}"))
            genes.append({
                'chrom_raw': seqid,
                'start': int(start),
                'end': int(end),
                'strand': strand,
                'gene_id': gid
            })
    return genes, chroms_seen

# --- Apply SNP edits and hash ---
def modify_chrom(seq: str, edits):
    if not edits:
        return seq
    arr = list(seq)
    n = len(arr)
    for pos0, base in edits:
        if 0 <= pos0 < n:
            arr[pos0] = base
    return ''.join(arr)

def gene_hashes(fasta_dict, vcf_edits, genes, raw2canon, canon2fasta):
    """Extract gene sequences from modified genome and hash them."""
    # Remap VCF edits by canonical key
    edits_by_canon = {}
    for raw_chrom, edit_list in vcf_edits.items():
        c = raw2canon.get(raw_chrom, _canon(raw_chrom))
        edits_by_canon.setdefault(c, []).extend(edit_list)
    # Sort and deduplicate merged edits
    for c in edits_by_canon:
        d = {}
        for p, a in edits_by_canon[c]:
            d[p] = a
        edits_by_canon[c] = sorted(d.items())

    modified_cache = {}
    out = []
    for g in genes:
        c = raw2canon.get(g['chrom_raw'], _canon(g['chrom_raw']))
        # Find FASTA key
        fa_key = canon2fasta.get(c)
        if fa_key is None:
            # Try exact match
            if g['chrom_raw'] in fasta_dict:
                fa_key = g['chrom_raw']
            else:
                continue
        # Build modified chrom once
        if fa_key not in modified_cache:
            edits = edits_by_canon.get(c, [])
            modified_cache[fa_key] = modify_chrom(fasta_dict[fa_key], edits)
        chrom_seq = modified_cache[fa_key]
        s0 = g['start'] - 1
        e0 = g['end']
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
    ap = argparse.ArgumentParser(
        description="vcf2hash: Gene-region allele profiling via SHA-256 hashing.",
        epilog="Example: python vcf2hash.py --fasta ref.fa --gff3 genes.gff3 --vcf sample.vcf.gz --out sample.tsv --tsv"
    )
    ap.add_argument('--fasta', required=True, help='Reference genome FASTA')
    ap.add_argument('--gff3', required=True, help='Gene annotations (GFF3)')
    ap.add_argument('--vcf', required=True, help='Sample VCF (single-sample preferred)')
    ap.add_argument('--out', default='-', help='Output file (default: stdout)')
    ap.add_argument('--tsv', action='store_true', help='Write TSV (gene_id\\thash) instead of JSONL')
    args = ap.parse_args()

    # Load inputs
    sys.stderr.write("Reading FASTA...\n")
    fasta = read_fasta(args.fasta)

    sys.stderr.write("Reading GFF3...\n")
    genes, gff_chroms = gff_genes(args.gff3)

    sys.stderr.write("Reading VCF...\n")
    edits, vcf_chroms = vcf_snps(args.vcf)

    # Build chromosome name mapping
    raw2canon, canon2fasta = build_chrom_map(fasta.keys(), gff_chroms, vcf_chroms)

    # Report matching stats
    fasta_canons = {_canon(k) for k in fasta.keys()}
    vcf_canons = {_canon(k) for k in vcf_chroms}
    gff_canons = {_canon(k) for k in gff_chroms}
    matched = fasta_canons & vcf_canons & gff_canons
    sys.stderr.write(f"Chromosomes: {len(fasta_canons)} FASTA, {len(gff_canons)} GFF3, {len(vcf_canons)} VCF, {len(matched)} matched\n")

    # Hash
    sys.stderr.write(f"Hashing {len(genes)} genes...\n")
    gh = gene_hashes(fasta, edits, genes, raw2canon, canon2fasta)
    sys.stderr.write(f"Output: {len(gh)} gene hashes\n")

    # Write output
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
