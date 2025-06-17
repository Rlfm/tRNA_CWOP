aa_dict = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    #'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z',
    #'XAA': 'X', 'UNK': 'X'
}

Wobble_pair = {
    "A":["T","C","G"],
    "T":["A","G","T"], # G not efficient + depended on base modif of U34
    "C":["G"], # confirmed
    "G":["C","T"], # confirmed 
    
    "I":["C","T"], # exceptionnaly A
    "k2C":["A"], # C34 of CAU tRNA Ile always modif
    "xm5s2U":["A"],
    "xm5Um":["A"],
    "Um":["A"],
    "xm5U": ["A"],
    "xo5U": ["T","A","G"]
    }

Mod_bases = {
    # I - Inosine (modif A)
    "TadA":True,
    # k2C - Lysidine
    "TilS":True,
    # xm5s2U: 5-(Carboxymethylaminomethyl)-2-thiouridine
    "mnmGorE": True,
    # xm5Um: 5-(Carboxymethylaminomethyl)-2′-O-methyluridine
    "TrmL" : True,
    # Um: 2′-O-methyluridine
    "trmJ" : True,
    # xm5U: 5-(Carboxymethylaminomethyl)uridine
    "trmA": True,
    # xo5U: 5-Hydroxyuridine
    "trhO": True
}

