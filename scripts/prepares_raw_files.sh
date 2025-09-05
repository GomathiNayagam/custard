#!/bin/bash

# Usage: ./rename_fastq.sh /path/to/your/fastq/files

if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/fastq/files"
    exit 1
fi

dir="$1"

if [ ! -d "$dir" ]; then
    echo "Error: Directory '$dir' does not exist."
    exit 1
fi

echo "Renaming FASTQ files in: $dir"

# Rename R1 files
for file in "$dir"/*_R1_001.fastq.gz; do
    [ -e "$file" ] || continue
    newname="${file/_R1_001.fastq.gz/_R1.fastq.gz}"
    mv "$file" "$newname"
done

# Rename R2 files
for file in "$dir"/*_R2_001.fastq.gz; do
    [ -e "$file" ] || continue
    newname="${file/_R2_001.fastq.gz/_R2.fastq.gz}"
    mv "$file" "$newname"
done

echo "Renaming completed."
