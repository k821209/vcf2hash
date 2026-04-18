# vcf2hash

Gene-region allele profiling via SHA-256 hashing for scalable digital breeding.

## Overview

vcf2hash converts multi-sample VCF files into gene-region haplotypes by:

1. **Extracting** SNPs within each annotated gene region
2. **Encoding** diploid genotypes as IUPAC characters (heterozygous: W, S, M, K, R, Y)
3. **Hashing** the resulting gene sequence with SHA-256

Identical gene-region sequences always produce the same hash — no phasing, no LD estimation, no centralized database required.

## Installation

```bash
git clone https://github.com/k821209/vcf2hash.git
cd vcf2hash
# No dependencies beyond Python 3.6+ standard library
```

## Usage

```bash
# Single sample
python3 vcf2hash.py \
    --fasta reference.fa \
    --gff3 genes.gff3 \
    --vcf sample.vcf.gz \
    --tsv --out sample_hashes.tsv

# Multiple samples with GNU parallel
parallel -j 10 "python3 vcf2hash.py \
    --fasta ref.fa --gff3 genes.gff3 \
    --vcf {} --tsv --out {}.tsv" ::: vcf_dir/*.vcf.gz
```

### Output

**TSV mode** (`--tsv`):
```
GeneA    a1b2c3d4e5f6...
GeneB    f6e5d4c3b2a1...
```

**JSONL mode** (default):
```json
{"gene_id":"GeneA","hash":"a1b2c3d4e5f6..."}
```

### Arguments

| Flag | Description |
|------|-------------|
| `--fasta` | Reference genome FASTA |
| `--gff3` | Gene annotations (GFF3) |
| `--vcf` | Sample VCF (single-sample preferred, supports .gz) |
| `--out` | Output file (default: stdout) |
| `--tsv` | Write TSV instead of JSONL |

## How it works

```
VCF (SNPs)  ──┐
               ├──▶  IUPAC-encoded gene sequence  ──▶  SHA-256 hash
FASTA + GFF3 ─┘
```

- **SNPs only** — INDELs excluded to maintain consistent sequence length
- **Heterozygous sites** — encoded as IUPAC ambiguity codes (e.g., A/T → W)
- **Missing data** — positions without variant calls retain the reference allele
- **Strand-aware** — reverse complement applied for genes on the minus strand
- **Species-agnostic** — chromosome names are automatically matched across FASTA, GFF3, and VCF

## Test

```bash
bash test/run_test.sh
```

Runs vcf2hash on two synthetic samples and verifies that identical gene sequences produce the same hash while different sequences produce different hashes.

## Citation

Park JS, Kim KD, Lee Y, Kang YJ. vcf2hash: Gene-Region Allele Profiling via SHA-256 Hashing for Scalable Digital Breeding. (in preparation)

## License

MIT
