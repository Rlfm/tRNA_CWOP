import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import pandas as pd
import numpy as np

aa_dict = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z',
    'XAA': 'X', 'UNK': 'X'
}

def read_tRNA_Scan(file):
    
    with open(file,"r") as f:
        lines = f.readlines()
        output = []
        output.append(f"Sequence Name,tRNA #,tRNA Begin,Bounds End,tRNA Type,AntiCodon,Intron Begin,Bounds End,Inf Score,Isotype CM,Isotype Score\n")
        for line in lines[3:]:
            line = [x.strip() for x in line.split("\t") if x != " "]
            to_output = []
            for x in line:
                to_output += x.split("\t")
            #print(to_output)
            output.append(f'{",".join(to_output[:-1])}\n')

    with open("tRNA.csv","w") as f:
        f.writelines(output)

    tRNA_df = pd.read_csv("tRNA.csv")
    tRNA_df['tRNA Type'] = tRNA_df['tRNA Type'].str.upper().map(aa_dict)
    tRNA_df.to_csv("tRNA.csv")
    return tRNA_df

def aa_codon_seq(aa,CDS):
    codon_df = pd.read_csv("codons.csv")
    aa_codons = codon_df[codon_df["aminoacid"] == aa]["codon"].to_list()

    codon_seq = list()

    for CDS in SeqIO.parse(CDS,"fasta"):
        codon_seq += [str(CDS.seq[i*3:(i+1)*3]) for i in range(round(len(CDS.seq)/3)) if CDS.seq[i*3:(i+1)*3] in aa_codons]

    return codon_seq

def build_transition_matrix(tRNA_abundance):
    keys = list(tRNA_abundance.keys())
    abundances = np.array([tRNA_abundance[k] for k in keys])
    
    # Create matrix where each row is the abundance vector
    T = np.tile(abundances, (len(keys), 1))
    
    # Normalize each row (so rows sum to 1)
    T = T / T.sum(axis=1, keepdims=True)
    
    # Return as pandas DataFrame for readability (optional)
    return pd.DataFrame(T, index=keys, columns=keys)
