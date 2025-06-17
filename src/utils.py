import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from itertools import combinations as find_combi
from configs import aa_dict,Wobble_pair,Mod_bases

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

def base_paring(ANTICODONS,ncb=False):

    codons_reco = dict.fromkeys(ANTICODONS)

    for anticodon in ANTICODONS:
        try: 
            codon_base_pair = "".join([Wobble_pair[x][0] for x in anticodon][::-1])
            """
            print(f"A 5' - {anticodon} - 3'")
            print(f"C 3' - {codon_base_pair[::-1]} - 5'")
            print(f"C 5' - {codon_base_pair} - 3'\n")
            """
            # Normal wooble 
            WT_codons_reco  = [f"{codon_base_pair[:-1]}{x}" for x in Wobble_pair[anticodon[0]]]
            
            if ncb:
                #Non-canonanical bases wooble 
                BM_codons_reco = []

                if Mod_bases["TadA"] and anticodon[0] == "A":
                    BM_codons_reco  += [f"{codon_base_pair[:-1]}{x}" for x in Wobble_pair["I"]]
                elif Mod_bases["TilS"] and anticodon[0] == "C":
                    BM_codons_reco  += [f"{codon_base_pair[:-1]}{x}" for x in Wobble_pair["k2C"]]

                elif anticodon[0] == "T":
                    if Mod_bases["mnmGorE"]:
                        BM_codons_reco  += [f"{codon_base_pair[:-1]}{x}" for x in Wobble_pair["xm5s2U"]]
                    if Mod_bases["TrmL"]:
                        BM_codons_reco  += [f"{codon_base_pair[:-1]}{x}" for x in Wobble_pair["xm5Um"]]
                    if Mod_bases["trmJ"]:
                        BM_codons_reco  += [f"{codon_base_pair[:-1]}{x}" for x in Wobble_pair["Um"]]
                    if Mod_bases["trmA"]:
                        BM_codons_reco  += [f"{codon_base_pair[:-1]}{x}" for x in Wobble_pair["xm5U"]]
                    if Mod_bases["trhO"]:
                        BM_codons_reco  += [f"{codon_base_pair[:-1]}{x}" for x in Wobble_pair["xo5U"]]
                
                codons_reco[anticodon] = list(dict.fromkeys(WT_codons_reco + BM_codons_reco))

            else:
                codons_reco[anticodon] = list(dict.fromkeys(WT_codons_reco))

        except KeyError:
            print(f"[ERROR] with anticodon {anticodon}")

    return codons_reco


def delta(a,b):

    assert a.shape == b.shape, f"{a.shape=} {b.shape=}"
    
    _delta = 0
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            delta_ij = abs(a[i][j] - b[i][j]) / b[i][j]
            _delta += delta_ij
    
    return _delta / (a.shape[0] * a.shape[1])