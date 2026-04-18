#!/bin/bash
# Test vcf2hash with synthetic data
# Expected: same SNPs = same hash, different SNPs = different hash

set -e
DIR="$(cd "$(dirname "$0")" && pwd)"
VCF2HASH="$DIR/../vcf2hash.py"

echo "=== vcf2hash test ==="
echo ""

echo "--- Sample 1 ---"
python3 "$VCF2HASH" \
    --fasta "$DIR/test_ref.fa" \
    --gff3 "$DIR/test_genes.gff3" \
    --vcf "$DIR/test_sample1.vcf" \
    --tsv --out "$DIR/out_sample1.tsv"
cat "$DIR/out_sample1.tsv"

echo ""
echo "--- Sample 2 ---"
python3 "$VCF2HASH" \
    --fasta "$DIR/test_ref.fa" \
    --gff3 "$DIR/test_genes.gff3" \
    --vcf "$DIR/test_sample2.vcf" \
    --tsv --out "$DIR/out_sample2.tsv"
cat "$DIR/out_sample2.tsv"

echo ""
echo "--- Comparison ---"
echo "Genes with SAME hash (identical allele profile):"
comm -12 <(sort "$DIR/out_sample1.tsv") <(sort "$DIR/out_sample2.tsv") || echo "(none)"
echo ""
echo "Genes with DIFFERENT hash:"
comm -3 <(sort "$DIR/out_sample1.tsv") <(sort "$DIR/out_sample2.tsv") | sed 's/^\t/  S2: /' | sed '/^[^\t ]/s/^/  S1: /'

# Cleanup
rm -f "$DIR/out_sample1.tsv" "$DIR/out_sample2.tsv"

echo ""
echo "=== Test complete ==="
