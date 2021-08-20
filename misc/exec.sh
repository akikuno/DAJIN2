#!/bin/sh

[ "${ZSH_VERSION:-}" ] && setopt shwordsplit interactivecomments

if ! conda info -e | grep -q DAJIN2; then
  conda config --add channels defaults
  conda config --add channels bioconda
  conda config --add channels conda-forge
  conda update -y -n base conda
  conda create -y -n DAJIN2
  conda install -y -n DAJIN2 minimap2 samtools numpy pandas scikit-learn hdbscan joblib plotnine
  conda install -y -n DAJIN2 -c conda-forge r-base r-essentials r-languageserver
fi
#------------------------------------------------------------------------------
# Execute test
#------------------------------------------------------------------------------

sh misc/exec_test.sh

#------------------------------------------------------------------------------
# pm-tyr
# chr7:87490801-87494587
#------------------------------------------------------------------------------
conda activate DAJIN2
sample=./example/pm-tyr/barcode31.fq.gz
control=./example/pm-tyr/barcode32.fq.gz
alleles=./example/pm-tyr/design_tyr.fa
threads=$(getconf _NPROCESSORS_ONLN | awk '{print $0-2}')
genome=mm10
control_name="$(echo ${control##*/} | sed "s/\..*$//" | tr " " "_")"
sample_name="$(echo ${sample##*/} | sed "s/\..*$//" | tr " " "_")"
sh misc/build.sh
###
# rm -rf .DAJIN_temp DAJIN_results /tmp/.DAJIN_temp
# sh misc/build.sh && time bash ./DAJIN -a "$alleles" -c "$control" -s "$sample" -g "$genome" -t "$threads"

#------------------------------------------------------------------------------
# Stx2
## chr5:128,991,721-128,996,151
#------------------------------------------------------------------------------
conda activate DAJIN2
sample=./example/del-stx2/barcode25.fq.gz
control=./example/del-stx2/barcode30.fq.gz
alleles=./example/del-stx2/design_stx2.fa
threads=$(getconf _NPROCESSORS_ONLN | awk '{num=$0-2; if(num>0) print num; else print 1}')
genome=mm10
control_name="$(echo ${control##*/} | sed "s/\..*$//" | tr " " "_")"
sample_name="$(echo ${sample##*/} | sed "s/\..*$//" | tr " " "_")"
sh misc/build.sh
bash -n ./DAJIN -a "$alleles" -c "$control" -s "$sample" -g "$genome" -t "$threads"
###
rm -rf .DAJIN_temp DAJIN_results /tmp/.DAJIN_temp
sh misc/build.sh && time bash ./DAJIN -a "$alleles" -c "$control" -s "$sample" -g "$genome" -t "$threads"

#------------------------------------------------------------------------------
# Ayabe Task1: Cables2
## chr5:128,991,721-128,996,151
#------------------------------------------------------------------------------

conda activate DAJIN2
sample=./example/flox-cables2/barcode31.fq.gz
control=./example/flox-cables2/barcode42.fq.gz
alleles=./example/flox-cables2/alleles.fa
threads=$(getconf _NPROCESSORS_ONLN | awk '{num=$0-2; if(num>0) print num; else print 1}')
genome=mm10
control_name="$(echo ${control##*/} | sed "s/\..*$//" | tr " " "_")"
sample_name="$(echo ${sample##*/} | sed "s/\..*$//" | tr " " "_")"
sh misc/build.sh
###
# rm -rf .DAJIN_temp DAJIN_results /tmp/.DAJIN_temp
# sh misc/build.sh && time bash ./DAJIN -a "$alleles" -c "$control" -s "$sample" -g "$genome" -t "$threads"

# Cables2
sample=./example/cables2-simulated/target_aligned_reads.fasta.gz
control=./example/cables2-simulated/wt_aligned_reads.fasta.gz
alleles=./example/cables2-simulated/design_cables2.fa
# target=CCT:GTCCAGAGTGGGAGATAGCC,CCA:CTGCTAGCTGTGGGTAACCC
threads=4
genome=mm10
control_name="$(echo ${control##*/} | sed "s/\..*$//" | tr " " "_")"
sample_name="$(echo ${sample##*/} | sed "s/\..*$//" | tr " " "_")"

# time ./DAJIN batch -a "$alleles" -c "$control" -s "$sample" -g "$genome" -t "$threads"

# Tyr-simulation
sample=./example/tyr-simulated/target_aligned_reads.fasta.gz
control=./example/tyr-simulated/wt_aligned_reads.fasta.gz
alleles=./example/tyr-simulated/design_tyr140.fa
# target=CCC:TGCGGCCAGCTTTCAGGCAG
threads=12
genome=mm10
control_name="$(echo ${control##*/} | sed "s/\..*$//" | tr " " "_")"
sample_name="$(echo ${sample##*/} | sed "s/\..*$//" | tr " " "_")"

time ./DAJIN batch \
  -a "$alleles" \
  -c "$control" \
  -s "$sample" \
  -t "$target" \
  -g "$genome" \
  -@ "$threads"

###############################################################################
# batch
###############################################################################

time find ../nanosim_fasta/ -type f |
  sort |
  while read -r sample; do

    echo "========================"
    echo "$sample"
    echo "========================"

    alleles=../design_cables2.fa
    control=../control/barcode42.fastq
    target=CCT:GTCCAGAGTGGGAGATAGCC,CCA:CTGCTAGCTGTGGGTAACCC
    threads=12
    genome=mm10

    time ./DAJIN batch \
      -a "$alleles" \
      -c "$control" \
      -s "$sample" \
      -t "$target" \
      -g "$genome" \
      -@ "$threads"
  done
