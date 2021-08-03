#!/bin/sh

[ "$#" -eq 0 ] && grep -v '```' .DAJIN_temp/document/usage.md && exit 0

while [ "$#" -gt 0 ]; do
  case "$1" in
  -h | --help)
    grep -v '```' .DAJIN_temp/document/usage.md
    exit 0
    ;;
  -v | --version)
    grep -v '```' .DAJIN_temp/document/version.md
    exit 0
    ;;

  #--- parse arguments
  -a | --alleles)
    alleles="$2"
    shift
    ;;
  -c | --control)
    control="$2"
    shift
    ;;
  -s | --sample)
    sample="$2"
    shift
    ;;
  -g | --genome)
    genome="$2"
    shift
    ;;
  -o | --output)
    output_dir="$2"
    shift
    ;;
  -t | --threads)
    threads="$2"
    shift
    ;;

  #--- error handling
  -*) echo "Unrecognized option: $1" && exit 1 ;;
  *) echo "Unrecognized argument: $1" && exit 1 ;;
  esac
  shift
done

################################################################################
# Check arguments
################################################################################

# Check mandatory arguments ------------

echo "$ARGS" |
  grep -q -e "\-a " -e "\--alleles " ||
  error_exit "-a|--alleles argument is required"

[ -r "$alleles" ] || error_exit "$alleles is not found"

grep -q -e '>wt' -e ">control" "$alleles" ||
  error_exit "$alleles must include '>control' or '>wt'"

echo "$ARGS" |
  grep -q -e "\-c " -e "\--control " ||
  error_exit "-c|--control argument is required"

[ -r "$control" ] || error_exit "$control is not found"

echo "$ARGS" |
  grep -q -e "\-s " -e "\--sample " ||
  error_exit "-s|--sample argument is required"

[ -r "$sample" ] || error_exit "$sample is not found"

control_name="$(basename "$control" | sed "s/\..*$//" | tr " " "_")"
sample_name="$(basename "$sample" | sed "s/\..*$//" | tr " " "_")"

[ "${control_name:-}" ] || error_exit "$control is an invalid file name."
[ "${sample_name:-}" ] || error_exit "$sample is an invalid file name"

# Check optional arguments ------------
