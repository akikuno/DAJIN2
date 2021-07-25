#!/bin/sh

#----------------------------------------------------------
echo "$(date +'%Y-%m-%d %H:%M:%S') Generate SAM files" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/sam
control_name="$(basename "$control" | sed "s/\..*$//" | tr " " "_")"
sample_name="$(basename "$sample" | sed "s/\..*$//" | tr " " "_")"

find .DAJIN_temp/fasta/*.fa |
  while read -r allele; do
    allele_name="$(basename ${allele%.fa})"

    if [ -s /tmp/"$control_name"_"$allele_name".sam ]; then
      cp /tmp/"$control_name"_"$allele_name".sam \
        .DAJIN_temp/sam/"$control_name"_"$allele_name".sam
    else
      minimap2 -t "$threads" -ax map-ont "$allele" "$control" --cs=long 2>/dev/null |
        tee /tmp/"$control_name"_"$allele_name".sam |
        cat >.DAJIN_temp/sam/"$control_name"_"$allele_name".sam
    fi

    minimap2 -t "$threads" -ax map-ont "$allele" "$sample" --cs=long 2>/dev/null |
      cat >.DAJIN_temp/sam/"$sample_name"_"$allele_name".sam
  done
