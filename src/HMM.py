import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from hmmlearn import hmm
import pandas as pd
from numpy import linalg as lg
from utils import read_tRNA_Scan,base_paring
from Bio import SeqIO
from configs import aa_dict
import dataframe_image as dfi

def aa_codon_seq(CDS,aa=False):
    codon_df = pd.read_csv("data/codons.csv")
    
    if aa:
        aa_codons = codon_df[codon_df["aminoacid"] == aa]["codon"].to_list()
    else:
        aa_codons = codon_df["codon"].to_list()
    
    codon_seq = list()
    for CDS in SeqIO.parse(CDS,"fasta"):
        codon_seq += [str(CDS.seq[i*3:(i+1)*3]) for i in range(round(len(CDS.seq)/3)) if CDS.seq[i*3:(i+1)*3] in aa_codons]

    return codon_seq

def build_transition_matrix(tRNA_abundance):
    keys = list(tRNA_abundance.keys())
    
    # Create matrix where each row is the abundance vector
    T = np.tile(list(tRNA_abundance.values()), (len(keys), 1))
    
    # Normalize each row (so rows sum to 1)
    T = T / T.sum(axis=1, keepdims=True)
    
    # Return as pandas DataFrame for readability (optional)
    transition_matrix = pd.DataFrame(T, index=keys, columns=keys)
    return transition_matrix

def build_emission_matrix(states,observations):
    
    """ shape (n_states,n_observations)
        codons| ACC | ACG | ACT
    tRNA       |     |     |    
    -----------------------------
    UGG        |     |  0.5| 0.5   
    -----------------------------
    UGC        |  1  |     |    

    """
    bp = base_paring(states)

    n_states = len(states)
    n_observations = len(observations)

    emission_matrix = np.zeros(shape=(n_states,n_observations))
    for i in range(n_states):
        anticodon = states[i]
        for j in range(n_observations):
            codon = observations[j]
            
            if codon in bp[anticodon]:
                emission_matrix[i][j] = 1
    
    emission_matrix += 1e-5
    emission_matrix = emission_matrix / np.sum(emission_matrix,axis=1,keepdims=True)
    #print(f"{self.emission_matrix=}")
    return emission_matrix


def HMM(CDS,aa=False,bias=False):

    # --- PARAMETERS ---

    try:
        tRNA_df = pd.read_csv("data/tRNA.csv")
    except FileNotFoundError:
        tRNA_df = read_tRNA_Scan("data/tRNA.txt")

    if aa:
        tRNA_df = tRNA_df[tRNA_df["tRNA Type"] == aa]
        states = tRNA_df["AntiCodon"].unique() # tRNAs
        n_states = len(states)
        #print('Number of hidden states :',n_states)

        codon_df = pd.read_csv("data/codons.csv")
        observations = codon_df[codon_df["aminoacid"] == aa]["codon"].unique() #Syn codons (regardless interval)
        n_observations = len(observations)
        #print('Number of observations  :',n_observations)
    else:
        states = tRNA_df["AntiCodon"].unique() # tRNAs
        n_states = len(states)
        #print('Number of hidden states :',n_states)

        codon_df = pd.read_csv("data/codons.csv")
        observations = codon_df["codon"].unique() #Syn codons (regardless interval)
        n_observations = len(observations)
        #print('Number of observations  :',n_observations)

    tRNA_abundance = {k:tRNA_df["AntiCodon"].to_list().count(k) for k in tRNA_df["AntiCodon"].to_list()}
        
    codon_seq = aa_codon_seq(CDS,aa)
    observations_sequence = np.array([np.where(observations == x) for x in codon_seq]).reshape(-1, 1)

    if bias:
        model = hmm.CategoricalHMM(n_components=n_states,init_params="")
        model.emissionprob_ = build_emission_matrix(states,observations)
        print(f"Before: {model.emissionprob_=}")

    else:
        model = hmm.CategoricalHMM(n_components=n_states,init_params="e")
        print("Before: em not set")

    model.transmat_ = build_transition_matrix(tRNA_abundance)
    model.startprob_ = model.get_stationary_distribution()

    print(f"Before: {model.transmat_=}")
    print(f"Before: {model.startprob_=}\n")

    model.fit(observations_sequence)
    
    print(f"After: {model.emissionprob_=}")
    print(f"After: {model.startprob_=}")
    print(f"After: {model.transmat_=}\n\n")

    return pd.DataFrame(data = model.emissionprob_,index=states,columns=observations)

def highligh_bp(df):

    assert isinstance(df, pd.DataFrame), f"Expected DataFrame, got {type(df)}"

    tRNAs = df.index.to_list()
    CODONS = df.columns.to_list()
    styles = pd.DataFrame('', index=df.index, columns=df.columns)

    std_bp = base_paring(tRNAs)
    ncb_bp = base_paring(tRNAs,ncb=True)

    for anticodon in ncb_bp:
        for codon in ncb_bp[anticodon]:
            if codon in CODONS:
                styles.loc[anticodon,codon] = 'background-color: red'

    for anticodon in std_bp:
        for codon in std_bp[anticodon]:
            if codon in CODONS:
                styles.loc[anticodon,codon] = 'background-color: cyan'

    return styles

def table_output_html(title,df):

    np.nan_to_num(df,copy=False)

    print(f"{df=}")

    # ---- Vizualisation ----    
    try:
        tRNA_df = pd.read_csv("data/tRNA.csv")
    except FileNotFoundError:
        tRNA_df = read_tRNA_Scan("data/tRNA.txt")

    """
    tRNAs = tRNA_df["AntiCodon"].unique()
    tRNAs = tRNAs[tRNAs != "NNN"]
    """    

    df = df.mask(df.abs() <= 1e-5, "")
    df = df.style.apply(highligh_bp,axis=None)

    html = df.to_html()
    #dfi.export(df, f"{file_name}", dpi=500)
    return f"""
    <div style="margin-bottom:30px;">
        <h3 style="font-family:sans-serif;">--- {title} ----</h3>
        {html}
    </div> 
    """



if __name__ == "__main__":

    std_emissions = list()
    biased_emissions = list()

    std_html = ""
    bias_html = ""

    amino_acids = list(aa_dict.values())[:3]
    inv_aa_dict = {v:k for k,v in aa_dict.items()}

    for aa in amino_acids:
        print(f"------------------- {aa} -------------------")
        std_em  = HMM("data/cds_from_genomic.fna",aa)
        std_em.to_csv(f"output/aa/std/{aa}.csv")
        
        std_html += table_output_html(inv_aa_dict[aa],std_em)
        std_emissions.append(std_em)

        bias_em = HMM("data/cds_from_genomic.fna",aa,bias=True)
        bias_em.to_csv(f"output/aa/bias/{aa}.csv")
        
        bias_html += table_output_html(inv_aa_dict[aa],bias_em)
        biased_emissions.append(bias_em)
    
    std_concat = pd.concat(std_emissions)
    std_concat.to_csv("output/std_concat.csv")
    std_concat_html = table_output_html("Concat all std init (Em not defined)",std_concat)
    std_html += std_concat_html

    with open(f"output/std.html","w") as f:
        f.write(std_html)

    bias_concat = pd.concat(biased_emissions)
    bias_concat.to_csv("output/bias_concat.csv")
    bias_concat_html = table_output_html("Concat all bias init (Em init based on wooble)",bias_concat)
    bias_html += bias_concat_html

    with open(f"output/bias.html","w") as f:
        f.write(bias_html)

    delta_std_bias = round((std_concat - bias_concat) / std_concat * 100)
    delta_std_bias.to_csv("output/delta_aa.csv")
    delta_std_bias_html = std_concat_html + bias_concat_html + table_output_html("AA model: delta between Em default or biased init (in %)",delta_std_bias)

    with open(f"output/delta_aa.html","w") as f:
        f.write(delta_std_bias_html)

    raise Exception("Not going full model")
    # ---- FULL MODEL ----

    all_std = HMM("data/cds_from_genomic.fna",aa=False,bias=False)
    all_std.to_csv("output/full_model_std.csv")
    all_std_html = table_output_html("Full model std, all tRNA + codons (Em not defined)",all_std)

    with open(f"output/full_model_std.html","w") as f:
        f.write(all_std_html)

    all_bias = HMM("data/cds_from_genomic.fna",aa=False,bias=True)
    all_bias.to_csv("output/full_model_bias.csv")
    all_bias_html = table_output_html("Full model bias, all tRNA + codons (Em not defined)",all_bias)

    with open(f"full_model_bias.html","w") as f:
        f.write(all_bias_html)

    full_delta_std_bias = round((all_std - all_bias) / all_std * 100)
    full_delta_std_bias.to_csv("output/delta_full.csv")
    full_delta_std_bias_html = all_std_html + all_bias_html + table_output_html("FULL model: delta between Em default or biased init (in %)",full_delta_std_bias)

    with open(f"output/delta_full.html","w") as f:
        f.write(full_delta_std_bias_html)
