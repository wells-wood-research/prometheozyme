#!/bin/bash

# Check if exactly one argument (directory path) is provided
if [ $# -ne 1 ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

dir="$1"

# Calculate total number of *.inp files across all categories
total_all_files=$(( \
  $(eval find "$dir" -path '*/opt/*.inp' -name "*.inp" | wc -l) + \
  $(eval find "$dir" -path '*/pull/*.inp' -name "*.inp" | wc -l) + \
  $(eval find "$dir" ! -path '*/opt/*' ! -path '*/pull/*' -name "*.inp" | wc -l) \
))
echo "Total files to process: $total_all_files"

# Function to process files with progress
process_files() {
  local category="$1"
  local find_path="$2"
  local exclude_paths="$3"
  local total_files="$4"

  echo "Processing $category files..."
  local current_file=0

  while IFS= read -r -d '' file; do
    current_file=$((current_file + 1))
    echo "Processing file $current_file of $total_files: $file"
    DISPLAY= /home/mchrnwsk/orca_6_0_1/orca "$file" > "${file%.inp}.out" 2>&1 < /dev/null
  done < <(eval find "$dir" $find_path $exclude_paths -name "*.inp" -print0)
}

# Process opt/*.inp files
process_files "opt/*.inp" "-path '*/opt/*.inp'" "" "$total_all_files"

# Process pull/*.inp files
process_files "pull/*.inp" "-path '*/pull/*.inp'" "" "$total_all_files"

# Process other *.inp files (excluding opt/ and pull/ subdirectories)
process_files "other *.inp" "" "! -path '*/opt/*' ! -path '*/pull/*'" "$total_all_files"

echo "All files processed."