import random
import json
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt


# Function to calculate GC content
def gc_content(seq):
    g = seq.count('G')
    c = seq.count('C')
    return (g + c) / len(seq) * 100


# Function to calculate melting temperature (Tm)
def calculate_tm(primer):
    return mt.Tm_NN(primer)


# Function to check if primer is valid
def is_valid_primer(primer, min_gc, max_gc, max_tm):
    # Check if GC content is within the required range
    if not (min_gc <= gc_content(primer) <= max_gc):
        return False

    # Check if the terminal bases are not GC
    if primer[0] in "GC" or primer[-1] in "GC":
        return False

    # Check if the melting temperature is within the required range
    tm = calculate_tm(primer)
    if tm > max_tm:
        return False

    return True


# Function to design primers
def design_primers(sequence, num_primers, min_len, max_len, min_gc, max_gc, max_tm):
    primers = []
    for _ in range(num_primers):
        while True:
            # Select a random region of the sequence for primer design
            start = random.randint(0, len(sequence) - max_len)
            end = start + random.randint(min_len, max_len)
            primer = sequence[start:end]

            # Check if it's a valid primer
            if is_valid_primer(primer, min_gc, max_gc, max_tm):
                primers.append(str(primer))
                break

    return primers


# Function to read configuration file (JSON format)
def read_config(config_file):
    with open(config_file) as f:
        config = json.load(f)
    return config


# Main function
def main():
    # Read config file
    config = read_config('primer_config.json')

    # Input sequence (for testing purposes)
    sequence = Seq(config["sequence"])

    # Design primers
    primers = design_primers(sequence, num_primers=config["num_primers"], min_len=config["min_len"], max_len=config["max_len"], min_gc=config["min_gc"], max_gc=config["max_gc"], max_tm=config["max_tm"])

    # Print the designed primers
    for primer in primers:
        print(f"Primer: {primer}")


if __name__ == "__main__":
    main()
