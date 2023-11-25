#! /bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <executable> <number_of_runs> <output_file>"
    exit 1
fi

executable="$1"
num_runs="$2"
output_file="$3"

echo "" > "$output_file"

# r5301c7 is the hardware code for counting Floating-point operations on Intel Skylake architecture
for((i=1; i<=$num_runs; i++)); do
    echo "Run $i:" >> "$output_file"
    perf stat -e r5301c7 ./"$executable" >> "$output_file" 2>&1
    echo "" >> "$output_file"
done

echo "Performance data collection halting. Data stored in $output_file"
