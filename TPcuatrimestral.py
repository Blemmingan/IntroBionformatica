from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

def find_orfs(sequence, min_length=100):
    """
    Finds open reading frames (ORFs) in a nucleotide sequence.
    The ORFs are translated into proteins and returned as a list.
    
    Parameters:
    - sequence (Seq): A nucleotide sequence object.
    - min_length (int): Minimum length of the ORF to be considered.
    
    Returns:
    - List of SeqRecords containing the ORF protein sequences.
    """
    orfs = []
    # Search for start and stop codons in the sequence
    start_codons = re.finditer(r'ATG', str(sequence))  # Standard start codon
    stop_codons = ['TAA', 'TAG', 'TGA']  # Standard stop codons

    # Find all ORFs in both strands
    for start in start_codons:
        start_pos = start.start()
        for stop_codon in stop_codons:
            stop_pos = sequence.find(stop_codon, start_pos + 3)
            if stop_pos != -1:
                # Extract the ORF between the start and stop codons
                orf = sequence[start_pos:stop_pos + 3]
                if len(orf) >= min_length:
                    # Translate the ORF into a protein sequence
                    protein = orf.translate()
                    orfs.append(SeqRecord(protein, id=f"orf_{start_pos}-{stop_pos}", description=""))
                break

    return orfs


def extract_orfs_from_genbank(genbank_file, output_fasta, output_whole_seq_fasta):
    """
    Extracts ORFs from a GenBank file and writes them in FASTA format.
    It also writes the whole sequence to a separate FASTA file.
    
    Parameters:
    - genbank_file (str): Path to the GenBank file.
    - output_fasta (str): Output FASTA file to store the protein sequences.
    - output_whole_seq_fasta (str): Output FASTA file to store the whole sequence.
    """
    with open(output_fasta, "w") as fasta_handle, open(output_whole_seq_fasta, "w") as whole_seq_handle:
        # Parse the GenBank file
        for record in SeqIO.parse(genbank_file, "genbank"):
            print(f"Processing GenBank record {record.id}...")
            
            # Extract and translate the ORFs from the sequence
            orfs = find_orfs(record.seq)
            
            # Write each ORF to the FASTA file
            SeqIO.write(orfs, fasta_handle, "fasta")
            print(f"{len(orfs)} ORFs found and written to {output_fasta}")
            
            # Write the whole sequence to a separate FASTA file
            whole_seq_record = SeqRecord(record.seq, id=record.id, description="Whole sequence")
            SeqIO.write(whole_seq_record, whole_seq_handle, "fasta")
            print(f"Whole sequence written to {output_whole_seq_fasta}")


# Example usage
if __name__ == "__main__":
    genbank_file = "sequence.gbk"
    output_fasta = "orfs.fas"  # Output FASTA file for ORFs
    output_whole_seq_fasta = "whole_sequence.fas"  # Output FASTA file for the whole sequence
    extract_orfs_from_genbank(genbank_file, output_fasta, output_whole_seq_fasta)

