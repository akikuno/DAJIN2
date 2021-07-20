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

# Check whether mandatory options are inputted ------------

echo "$ARGS" |
  grep -q -e "\-a " -e "\--alleles " ||
  error_exit "-a|--alleles argument is required"

echo "$ARGS" |
  grep -q -e "\-c " -e "\--control " ||
  error_exit "-c|--control argument is required"

echo "$ARGS" |
  grep -q -e "\-s " -e "\--sample " ||
  error_exit "-s|--sample argument is required"
