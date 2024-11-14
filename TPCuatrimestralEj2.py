from Bio.Blast import NCBIWWW
from Bio import SeqIO
import sys

def blast_protein_sequences(input_fasta, max_hits=10):
    """
    Runs BLAST for each protein sequence in the input FASTA file and saves the results in a fixed output XML file.
    Stops as soon as we collect a total of `max_hits` hits.

    Parameters:
    - input_fasta (str): Path to the FASTA file containing protein sequences (ORFs).
    - max_hits (int): Maximum number of hits to collect (default is 10).
    """
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    output_file = "blast_results.xml"  # Fixed output file name
    
    # Open the output file for writing
    with open(output_file, "w") as out_handle:
        all_hits = 0  # Track the number of hits
        for i, record in enumerate(sequences, start=1):
            print(f"Running BLAST for sequence {record.id} ({i}/{len(sequences)})...")
            
            # Perform the BLAST search using NCBI's `blastp` against the `swissprot` database
            result_handle = NCBIWWW.qblast("blastp", "swissprot", record.seq)
            
            # Write the BLAST result directly to the output file
            out_handle.write(result_handle.read())
            
            # Count the number of hits in this BLAST result
            result_handle.seek(0)  # Reset the result handle to start reading from the beginning
            from Bio.Blast import NCBIXML
            blast_records = NCBIXML.parse(result_handle)
            hits_in_record = 0
            for blast_record in blast_records:
                hits_in_record = len(blast_record.alignments)
            
            all_hits += hits_in_record
            print(f"Found {hits_in_record} hits for sequence {record.id}. Total hits so far: {all_hits}")
            
            # Stop if we have collected enough hits
            if all_hits >= max_hits:
                print(f"Collected {all_hits} hits. Stopping search.")
                break

    print(f"BLAST analysis completed. Results saved to {output_file}")

# Command-line execution
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_fasta>", file=sys.stderr)
        sys.exit(1)

    input_fasta = sys.argv[1]  # Input FASTA file (protein sequences)
    
    blast_protein_sequences(input_fasta)

