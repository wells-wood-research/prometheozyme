#!/bin/bash

# Check if exactly one argument (directory path) is provided
if [ $# -ne 1 ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

dir="$1"

# Function to process files with progress
process_files() {
  local category="$1"
  local find_path="$2"
  local exclude_paths="$3"

  echo "Processing $category files..."
  local total_files
  total_files=$(find "$dir" $find_path $exclude_paths -name "*.inp" | wc -l)
  local current_file=0

  while IFS= read -r -d '' file; do
    current_file=$((current_file + 1))
    echo "Processing file $current_file of $total_files: $file"
    /home/mchrnwsk/orca_6_0_1/orca "$file" > "${file%.inp}.out"
  done < <(find "$dir" $find_path $exclude_paths -name "*.inp" -print0)
}

# Process opt/*.inp files
process_files "opt/*.inp" "-path '*/opt/*.inp'"

# Process pull/*.inp files
process_files "pull/*.inp" "-path '*/pull/*.inp'"

# Process other *.inp files (excluding opt/ and pull/ subdirectories)
process_files "other *.inp" "" "! -path '*/opt/*' ! -path '*/pull/*'"

echo "All files processed."