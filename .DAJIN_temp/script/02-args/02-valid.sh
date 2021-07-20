#!/bin/sh

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
