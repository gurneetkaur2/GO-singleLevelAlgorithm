#!/bin/sh

# Configuration stuff

fspec=links-anon.txt
num_files=2

# Work out lines per file.

total_lines=$(wc -l <${fspec})
( (lines_per_file = (total_lines + num_files - 1) / num_files) )

echo "Total lines     = ${total_lines}"
echo "Lines  per file = ${lines_per_file}"    
# Split the actual file, maintaining lines.

split --lines=${lines_per_file} ${fspec} links.

# Debug information

echo "Total lines     = ${total_lines}"
echo "Lines  per file = ${lines_per_file}"    
wc -l links.*
