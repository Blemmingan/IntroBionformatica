from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

def extract_top_hits(blast_xml_file, num_hits=10):
    """
    Parses a BLAST XML file, extracts the top sequences (by score), 
    and returns them as a list of sequences.
    
    Parameters:
    - blast_xml_file (str): Path to the BLAST XML output file.
    - num_hits (int): Number of top sequences to retrieve (default is 10).
    
    Returns:
    - top_sequences (list of SeqRecord): List of top sequences extracted from BLAST results.
    """
    top_sequences = []

    with open(blast_xml_file, "r") as blast_handle:
        blast_records = NCBIXML.parse(blast_handle)
        
        for record in blast_records:
            sorted_hits = sorted(record.alignments, key=lambda x: x.hsps[0].score, reverse=True)
            
            for hit in sorted_hits[:num_hits]:  # Get top N hits
                seq = hit.hsps[0].sbjct
                top_sequences.append(seq)

    return top_sequences

def save_sequences_to_fasta(sequences, output_fasta):
    """
    Saves the list of sequences to a FASTA file.
    
    Parameters:
    - sequences (list): List of sequences to save.
    - output_fasta (str): Path to the output FASTA file.
    """
    with open(output_fasta, "w") as output_handle:
        for i, seq in enumerate(sequences):
            record = SeqRecord(Seq(seq), id=f"Hit_{i+1}", description="")
            SeqIO.write(record, output_handle, "fasta")
    
    print(f"Sequences have been saved to {output_fasta}.")

def align_sequences(input_fasta, output_fasta):
    """
    Aligns sequences using Clustal Omega via the command line.
    
    Parameters:
    - input_fasta (str): Input FASTA file with sequences to align.
    - output_fasta (str): Output file to save the aligned sequences.
    """
    clustalw_cline = f"clustalo -i {input_fasta} -o {output_fasta} --force"
    subprocess.run(clustalw_cline, shell=True)
    print(f"Alignment complete. Results saved to {output_fasta}.")

if __name__ == "__main__":
    # Define file paths
    blast_xml_file = "blast_results.xml"
    top_fasta_file = "top_10_sequences.fasta"
    aligned_fasta_file = "aligned_sequences.fasta"
    
    # Extract top hits from BLAST XML
    top_sequences = extract_top_hits(blast_xml_file)
    
    # Save the top sequences to a FASTA file
    save_sequences_to_fasta(top_sequences, top_fasta_file)
    
    # Perform multiple sequence alignment using Clustal Omega
    align_sequences(top_fasta_file, aligned_fasta_file)
