#!/bin/bash

# Script to extract the LAST FINAL SINGLE POINT ENERGY from pull.out files,
# extract the arrangement name from the path, and print the energy and arrangement
# sorted by energy from lowest to highest.

if [ $# -ne 1 ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

dir="$1"

echo "Extracting last energies and sorting by arrangement..."

find "$dir" -path "*/pull/pull.out" -type f | while read -r file; do
  energy=$(grep "FINAL SINGLE POINT ENERGY" "$file" | tail -n 1 | awk '{print $NF}')
  if [ -n "$energy" ]; then
    # Extract the arrangement name (directory before /pull)
    arrangement=$(basename "$(dirname "$(dirname "$file")")")
    echo "$energy $arrangement"
  else
    echo "Warning: Energy not found in $file" >&2
  fi
done | sort -g