#!/bin/bash

# Check if the first Python script exists, else exit immediately.
if [ ! -f "TPcuatrimestral.py" ]; then
    echo "Python script TPcuatrimestral.py not found!"
    exit 1
fi

# Run the first Python script to generate the FASTA file
python3 TPcuatrimestral.py

# Check if the first Python script executed successfully
if [ $? -eq 0 ]; then
    echo "First Python script executed successfully."
else
    echo "An error occurred while running the first Python script."
    exit 1
fi

# Define the generated FASTA file (assuming it's called orfs.fas)
input_fasta="orfs.fas"

# Verify if the input file exists
if [[ ! -f "$input_fasta" ]]; then
    echo "Error: The input file '$input_fasta' does not exist."
    exit 1
fi

############## End of the first part

# Output file for BLAST results
output_file="blast_results.xml"

# Verify if the output file already exists and ask whether to overwrite
if [[ -f "$output_file" ]]; then
    read -p "The output file '$output_file' already exists. Do you want to overwrite it? (y/n): " overwrite
    if [[ ! "$overwrite" =~ ^[yY]$ ]]; then
        echo "Exiting without overwriting the file."
        exit 1
    fi
fi

# Run the Python script for BLAST analysis
echo "Running BLAST, will stop at 10 hits or it will take a very long time."
echo "Running BLAST for the file '$input_fasta'..."

# Call the Python script to run BLAST
python3 TPCuatrimestralEj2.py "$input_fasta"

# Check if the BLAST Python script executed successfully
if [[ $? -ne 0 ]]; then
    echo "Error: BLAST analysis failed for '$input_fasta'."
    exit 1
fi

echo "BLAST analysis completed successfully for '$input_fasta'. Results saved to '$output_file'."

# Now, call the second Python script to align the sequences from BLAST results
echo "Running alignment for the BLAST results..."

# Call the Python script TPCuatrimestralEj3.py to align the sequences
python3 TPCuatrimestralEj3.py "$output_file"

# Check if the alignment Python script executed successfully
if [[ $? -ne 0 ]]; then
    echo "Error: Alignment failed using '$output_file'."
    exit 1
fi

echo "Alignment completed successfully. The aligned sequences have been saved."

############## New section to call the last Python script (TPCuatrimestralEj4.py)

# Finally, call the TPCuatrimestralEj4.py script after everything else is done
echo "Running final Python script TPCuatrimestralEj4.py..."

python3 TPCuatrimestralEj4.py

# Check if the final Python script executed successfully
if [[ $? -ne 0 ]]; then
    echo "Error: Final Python script TPCuatrimestralEj4.py failed."
    exit 1
fi

echo "Final analysis completed successfully by TPCuatrimestralEj4.py."

############## New section to call the Python_ej5.py script

# Finally, call the Python_ej5.py script
echo "Running Python script Python_ej5.py..."

python3 Python_ej5.py > primers_results.txt

# Check if the Python_ej5.py script executed successfully
if [[ $? -ne 0 ]]; then
    echo "Error: Python script Python_ej5.py failed."
    exit 1
fi

echo "Python script Python_ej5.py executed successfully."

