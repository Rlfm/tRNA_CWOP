import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from hmmlearn import hmm
import pandas as pd
from numpy import linalg as lg
from utils import aa_codon_seq,read_tRNA_Scan,build_transition_matrix



def HMM(states,observations,Tr_matrix,Em_matrix,codon_seq):
    
    n_states = len(states)

    # Transpose P because we're looking for left eigenvector: πᵗ P = πᵗ → (Pᵗ πᵗ = πᵗ)
    eigvals, eigvecs = np.linalg.eig(Tr_matrix.T)
    # Find the eigenvector corresponding to eigenvalue 1
    stationary_index = np.argmin(np.abs(eigvals - 1))
    pi = np.real(eigvecs[:, stationary_index])  # get the real part in case of numerical noise
    # Normalize to sum to 1
    pi /= np.sum(pi)

    print("Stationary distribution (π):", pi)

    state_probability = pi
    print("State probability: ", state_probability)
    
    print("\nEmission probability:\n", Em_matrix)

    # --- MODEL ---
    model = hmm.CategoricalHMM(n_components=n_states)
    model.startprob_ = state_probability
    model.transmat_ = Tr_matrix
    model.emissionprob_ = Em_matrix
    
    observations_sequence = np.array([np.where(observations == x) for x in codon_seq]).reshape(-1, 1)
    observations_sequence

    # Predict the most likely sequence of hidden states
    
    log_probability, hidden_states = model.decode(observations_sequence,
                                                lengths = len(observations_sequence),
                                                algorithm ='viterbi' )

    print('Log Probability :',log_probability)
    print("Most likely hidden states:", hidden_states)

    return log_probability,hidden_states

def main(aa,CDS):

    # --- PARAMETERS ---

    try:
        tRNA_df = pd.read_csv("tRNA.csv")
    except FileNotFoundError:
        tRNA_df = read_tRNA_Scan("tRNA.txt")

    tRNA_df = tRNA_df[tRNA_df["tRNA Type"] == aa]

    states = tRNA_df["AntiCodon"].unique() # tRNAs
    n_states = len(states)
    print('Number of hidden states :',n_states)

    codon_df = pd.read_csv("codons.csv")

    observations = codon_df[codon_df["aminoacid"] == aa]["codon"].unique() #Syn codons (regardless interval)
    n_observations = len(observations)
    print('Number of observations  :',n_observations)

    # ------------- I - Converge transtion matrix

    tRNA_abundance = dict.fromkeys(tRNA_df["AntiCodon"].to_list())
    tRNA_abundance = {k:tRNA_df["AntiCodon"].to_list().count(k) for k in tRNA_abundance.keys()}

    transition_probability = build_transition_matrix(tRNA_abundance)
    print("Initial Transition probability:\n", transition_probability)
    
    """ shape (n_states,n_observations)
        codons| ACC | ACG | ACT
    tRNA       |     |     |    
    -----------------------------
    UGG        |     |  0.5| 0.5   
    -----------------------------
    UGC        |  1  |     |    

    """
    emission_probability = np.ones(shape=(n_states,n_observations)) / n_observations
    codon_seq = aa_codon_seq(aa,CDS)

    delta = 10**99
    i = 0
    while delta > 0.001:
        print(f" ------------------------ iteration {i} ------------------------")
        log_probability,hidden_states = HMM(states,observations,transition_probability,emission_probability,codon_seq)
        print(f"{hidden_states=}")

        states_transitions = [(hidden_states[i],hidden_states[i+1]) for i in range(len(hidden_states) - 2)]
        new_T = np.zeros(shape=(n_states,n_states))
        for state0,state1 in states_transitions:
            new_T[state0][state1] += 1
        
        new_T = new_T / new_T.sum(axis=1,keepdims=True)
        np.nan_to_num(new_T,copy=False)
        
        assert new_T.shape == transition_probability.shape, f"{new_T.shape=} {transition_probability.shape=}"

        for i in range(new_T.shape[0]):
            for j in range(new_T.shape[1]):
                delta += abs(new_T[i][j] - transition_probability.values[i][j]) / transition_probability.values[i][j]

        delta = delta / (new_T.shape[0] * new_T.shape[1])
        print(f"{new_T=}")
        print(f"{transition_probability=}")
        print(f"{delta=}")

        transition_probability = pd.DataFrame(new_T, index=states, columns=states)

if __name__ == "__main__":
    main("A","cds_from_genomic.fna")