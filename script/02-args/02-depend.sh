#!/bin/sh

################################################################################
# Chech dependencies
################################################################################

# curl or wget ------------------------------------------------

if type wget >/dev/null 2>&1; then
  CMD_CHECK='wget -q -O - --spider --tries=2 --wait=1 --timeout=5'
  CMD_GET='wget -q -O -'
elif type curl >/dev/null 2>&1; then
  CMD_CHECK='curl --retry 2 --retry-delay 1 -s -o /dev/null -w "%{http_code}"'
  CMD_GET='curl -s'
else
  error_exit 'No HTTP-GET/POST command found.'
fi

# Commands ------------------------------------------------

cat <<EOF |
  minimap2
  samtools
  gzip
  python
EOF
  while read -r CMD; do
    "$CMD" --version >/dev/null 2>&1 || error_exit "$CMD: command not found"
  done

# R packages ----------------------------------------------
# Rscript --slave --vanilla .DAJIN_temp/library/install_pkgs.R 1>/dev/null 2>&1

# Rscript --slave --vanilla -e "installed.packages()" >.DAJIN_temp/rpackages

# cat <<EOF |
#   dbscan
#   ggplot2
#   RColorBrewer
# EOF
#   while read -r RPKG; do
#     grep -q "$RPKG" .DAJIN_temp/rpackages || error_exit "$RPKG: package not found in R"
#   done

# rm .DAJIN_temp/rpackages

# Python packages ----------------------------------------

pip freeze >.DAJIN_temp/pypackages

cat <<EOF |
  numpy
  pandas
  scikit-learn
  hdbscan
  joblib
  plotnine
EOF
  while read -r PYPKG; do
    grep -q "$PYPKG" .DAJIN_temp/pypackages || error_exit "$PYPKG: package not found in Python"
  done

rm .DAJIN_temp/pypackages
