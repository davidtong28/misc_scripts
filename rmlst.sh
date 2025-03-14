#!/bin/bash



# Check if only one input file and one output is provided
if [ $# -lt 2 ]; then
  echo "takes in one fasta file run rMLST using API; then parses the json result file"
  echo "requires curl, R version 4.0+, R packages: jsonlite, dplyr, stringr, argparse"
  echo "Usage: $0 -i <input_fasta_file> -p <output_prefix> "
  exit 1
fi

# argument parsing
while getopts ":i:p:" opt; do
  case $opt in
    i)
      input="$OPTARG"
      ;;
    p)
      prefix="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done



if [ -z "$input" ]; then
  echo "Option -i is mandatory. Please provide it."
  exit 1
fi

if [ -z "$prefix" ]; then
  echo prefix is defaulted to $(realpath $input|sed -n 's/^\(.*\)\..*$/\1/p')
  prefix=$(realpath $input|sed -n 's/^\(.*\)\..*$/\1/p')
fi

# will not check for existance, defaults to overwrite

#perform rMLST on API

(echo -n '{"base64":true,"details":true,"sequence": "'; base64 $input ; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence" -d @- > $prefix.json





