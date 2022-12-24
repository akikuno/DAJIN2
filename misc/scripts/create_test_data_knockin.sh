#!/bin/bash

# - barcode31,flox
# - barcode32,1nt-left-substitution
# - barcode33,1nt-deletion
# - barcode34,1nt-insertion

# samtools view DAJINResults/Ayabe-Task1/BAM/barcode31/barcode31_allele1_flox_intact_31.906%.bam | wc -l
# samtools view DAJINResults/Ayabe-Task1/BAM/barcode32/barcode32_allele1_flox_intact_31.906%.bam | wc -l
# samtools view DAJINResults/Ayabe-Task1/BAM/barcode33/barcode33_allele1_flox_intact_36.029%.bam | wc -l
# samtools view DAJINResults/Ayabe-Task1/BAM/barcode34/barcode34_allele1_flox_intact_35.361%.bam | wc -l

# samtools view DAJINResults/Ayabe-Task1/BAM/barcode31/barcode31_allele1_flox_intact_31.906%.bam | grep -c GCTATACGAAGTTAT
# samtools view DAJINResults/Ayabe-Task1/BAM/barcode32/barcode32_allele1_flox_intact_31.906%.bam | grep -c GCTCTACGAAGTTAT
# samtools view DAJINResults/Ayabe-Task1/BAM/barcode33/barcode33_allele1_flox_intact_36.029%.bam | grep -c GCTTACGAAGTTAT
# samtools view DAJINResults/Ayabe-Task1/BAM/barcode34/barcode34_allele1_flox_intact_35.361%.bam | grep -c GCTGATACGAAGTTAT

# echo ATAACTTCGTATCAGC | tr ACGT TGCA | rev
# echo ATAACTTCGTAGAGC | tr ACGT TGCA | rev

# samtools view DAJINResults/Ayabe-Task1/BAM/barcode31/barcode31_allele1_flox_intact_31.906%.bam | wc -l
# samtools view DAJINResults/Ayabe-Task1/BAM/barcode32/barcode32_allele1_flox_intact_31.906%.bam | wc -l
# samtools view DAJINResults/Ayabe-Task1/BAM/barcode33/barcode33_allele1_flox_intact_36.029%.bam | wc -l
# samtools view DAJINResults/Ayabe-Task1/BAM/barcode34/barcode34_allele1_flox_intact_35.361%.bam | wc -l

cat <<EOF |
barcode31 GCTATACGAAGTTAT
barcode32 GCTCTACGAAGTTAT
barcode33 GCTTACGAAGTTAT
barcode34 GCTGATACGAAGTTAT
EOF
    while read -r barcode sequence; do
        samtools view DAJINResults/Ayabe-Task1/BAM/"$barcode"/"$barcode"_allele1_*.bam |
            grep "$sequence" |
            cut -f 1 |
            sort -u |
            head -n 100 >tmp_"$barcode"
        wc -l tmp_"$barcode"
        gzip -dc examples/flox-cables2/AyabeTask1/"$barcode".fq.gz |
            paste - - - - |
            grep -f tmp_"$barcode" |
            tr "\t" "\n" >tmp_"$barcode".fq
    done

# - barcode31,flox
# - barcode32,1nt-left-substitution
# - barcode33,1nt-deletion
# - barcode34,1nt-insertion

cat tmp_barcode31.fq tmp_barcode32.fq | gzip -c >test_flox_leftsub.fq.gz
cat tmp_barcode31.fq tmp_barcode33.fq | gzip -c >test_flox_del.fq.gz
cat tmp_barcode31.fq tmp_barcode34.fq | gzip -c >test_flox_ins.fq.gz
