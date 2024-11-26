import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def run_getorf(input_fasta, output_fasta):
    """
    Llama al programa getorf de EMBOSS para encontrar ORFs en una secuencia de nucleótidos
    y los escribe en un archivo de salida en formato FASTA.

    Parameters:
    - input_fasta (str): archivo FASTA de entrada con la secuencia de nucleótidos.
    - output_fasta (str): archivo FASTA de salida para las secuencias de proteínas (ORFs).
    """
    # Ejecutar EMBOSS getorf
    
    subprocess.run(['getorf', '-sequence', input_fasta, '-outseq', output_fasta, '-minsize', '100'], check=True)

def run_patmatmotifs(protein_fasta, output_file, prosite_db_path):
    """
    Realiza el análisis de dominios en las secuencias de proteínas utilizando patmatmotifs de EMBOSS,
    usando la base de datos PROSITE.
	

    Parameters:
    - protein_fasta (str): archivo FASTA con las secuencias de proteínas.
    - output_file (str): archivo de salida para los resultados del análisis de dominios.
    - prosite_db_path (str): ruta al archivo prosite.dat con la base de datos PROSITE.
    """
    with open(output_orfs, 'r') as f:
    	sequences = f.read().split('>')[1:]
    	for seq in sequences:
       	 header, sequence = seq.split('\n', 1)
       	 sequence = ''.join(sequence.split())
         print(f"Sequence length: {len(sequence)}")

    # Ejecutar EMBOSS patmatmotifs para buscar dominios en las proteínas
    subprocess.run(['patmatmotifs', '-sequence', protein_fasta, '-outfile', output_file, '-prune=false'], check=True)


def main(input_fasta, output_orfs, output_patmatmotifs_results, prosite_db_path):
    """
    Función principal para coordinar el análisis de ORFs y la búsqueda de dominios
    usando EMBOSS.

    Parameters:
    - input_fasta (str): archivo FASTA de entrada con la secuencia de nucleótidos.
    - output_orfs (str): archivo FASTA de salida con las secuencias de proteínas (ORFs).
    - output_patmatmotifs_results (str): archivo de salida para los resultados del análisis de dominios.
    - prosite_db_path (str): ruta al archivo prosite.dat.
    """
    # Paso 1: Ejecutar getorf para encontrar los ORFs y obtener las proteínas
    print("Encontrando ORFs y traduciendo las secuencias...")
    run_getorf(input_fasta, output_orfs)
    
    # Paso 2: Ejecutar patmatmotifs para buscar los dominios en las secuencias de proteínas
    print("Buscando dominios en las proteínas usando PROSITE...")
    run_patmatmotifs(output_orfs, output_patmatmotifs_results, prosite_db_path)
    
    print(f"Análisis completado. Los resultados de dominios se guardaron en: {output_patmatmotifs_results}")

if __name__ == "__main__":
    # Ruta al archivo FASTA de entrada con la secuencia de nucleótidos
    input_fasta = "whole_sequence.fas"
    
    # Archivo de salida para las secuencias de proteínas (ORFs)
    output_orfs = "orfs_proteinas.fasta"
    
    # Archivo de salida para los resultados del análisis de dominios (PROSITE)
    output_patmatmotifs_results = "dominios_resultados_patmatmotifs.txt"
    
    # Ruta al archivo prosite.dat (base de datos PROSITE)
    prosite_db_path = "PROSITE/prosite.dat" 
    
    # Ejecutar el análisis
    main(input_fasta, output_orfs, output_patmatmotifs_results, prosite_db_path)

