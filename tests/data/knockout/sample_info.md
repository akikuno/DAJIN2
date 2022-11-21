# sample infomation

- Stx2 2-cut deletion

## test_barcode25.fq.gz
- Three alleles of large deletions
- 500 reads/allele at 1:1:1 ratio

## test_barcode30.fq.gz
- Wild type control
- 2000 reads


# Design

- deletion: 2012-2739 (727 bases)

- deletion size mapped to deletion alleles
  - allele1: ~260 bp deletion
  - allele2: ~620 bp deletion
  - allele3: ~1800 bp deletion

To examine deletion size

```bash
cat tests/data/knockout/design_stx2.fa |
    grep -A 1 deletion > tmp_del.fa

minimap2 -ax map-ont tmp_del.fa tests/data/knockout/test_barcode25.fq.gz |
samtools sort > tmp_bc25.bam
samtools index tmp_bc25.bam
```
