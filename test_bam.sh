. library/maskFastq.sh
gzip -dc $control >tmp_control.fq
gzip -dc $control | maskFastq >tmp_control_mask.fq

gzip -dc $sample >tmp_sample.fq
gzip -dc $sample | maskFastq >tmp_sample_mask.fq

for fq in *fq; do
  minimap2 -ax map-ont .DAJIN_temp/fasta/control.fa $fq >${fq%.fq}.sam
done

for sam in *sam; do
  samtools sort $sam >${sam%.sam}.bam
  samtools index ${sam%.sam}.bam
done

minimap2 -ax map-ont .DAJIN_temp/fasta/control.fa tmp.fq >tmp.sam
minimap2 -ax map-ont .DAJIN_temp/fasta/control.fa tmp2.fq >tmp2.sam

head tmp.sam | cut -f 1-5
head tmp2.sam | cut -f 1-5

cat tmp.sam | awk '$2==0 || $2==16' | wc -l
cat tmp2.sam | awk '$2==0 || $2==16' | wc -l

cat tmp.sam | grep 01664beb-54e6-46e1-ad56-7bbec663abc2
cat tmp2.sam | grep 01664beb-54e6-46e1-ad56-7bbec663abc2
