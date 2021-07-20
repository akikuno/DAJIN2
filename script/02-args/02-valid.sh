#!/bin/sh

# Check whether mandatory options are inputted ------------

[ -z "$design" ] && error_exit "design argument is not specified"
[ -z "$input_dir" ] && error_exit "input_dir argument is not specified"
[ -z "$control" ] && error_exit "control argument is not specified"

echo "$ARGS" |
  grep -q -e "\-a " -e "\--alleles " ||
  error_exit "-a|--alleles argument is required"

echo "$ARGS" |
  grep -q -e "\-c " -e "\--control " ||
  error_exit "-c|--control argument is required"

echo "$ARGS" |
  grep -q -e "\-s " -e "\--sample " ||
  error_exit "-s|--sample argument is required"
