import sys
import pandas as pd
from typing import Dict, List, Sequence
from collections import Counter
from Bio import SeqIO
from pathlib import Path

def read_edgelist(fName):
    """Reads a network in edgelist format"""
    Network = set()
    with open(fName, 'r') as f1:
        for line in f1:
            Ary = line.strip().split()
            if f"{Ary[0]}\t{Ary[1]}" in Network or f"{Ary[1]}\t{Ary[0]}" in Network:
                continue
            
            Network.add(f"{Ary[0]}\t{Ary[1]}")
    return(Network)

def write_edgelist(Network, fName, OutDir):
    """Takes a newtork and write it in edgelist format"""
    fOut = open(OutDir + "/" + fName + '_ex.str', 'w')
    for i in Network:
        print(i, file = fOut)
    fOut.close()

def read_list(fName):
    """Read a file as a list where each line represents one of its elements"""
    try:
        with open(fName, 'rt') as f1:
            return f1.read().strip().split('\n')
    except OSError as e:
        print("I/O error({0}): {1} {2}".format(e.errno, e.strerror, fName))
        sys.exit(1)

def write_dict(dict_, fOut = 'Output.txt'):
    with open(fOut, 'wt') as f1:
        for k, v in dict_.items():
            print(k, v, sep = '\t', file = f1)

def progress_bar(i, total):
    # write progress----
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s]" % ('='*int(50* i / total)))
    sys.stdout.flush()
    
def orca_input(Network: Dict[str, List[str]], fName: str) -> None:
    Ncount = 0
    NMapping: Dict[str, int] = dict()
    Net = dict()
    for i in Network:
        Ary = i.strip().split()
        
        if Ary[0] not in NMapping:
            NMapping[Ary[0]] = Ncount
            Ncount += 1

        if Ary[1] not in NMapping:
            NMapping[Ary[1]] = Ncount
            Ncount += 1
        
        Net[f'{NMapping[Ary[0]]} {NMapping[Ary[1]]}\n'] = [len(Network[i])]
    
    write_dict(NMapping, 'NMapping.txt')
    with open(fName, 'w') as fOut:
        fOut.write(f'{len(NMapping)} {len(Net)}\n')
        fOut.write("".join(Net.keys()))
    
    return NMapping

def read_annots(S2Q: Dict[str, str], FAnnots: str, NMapping: Dict[str, int]) -> Dict:
    Annots: Dict[str, List[str]] = dict()
    # Collect the annotaions for uids in S2Q
    with open(FAnnots, 'rt') as f1:
        for line in f1:
            Ary = line.strip().split('\t')
            if Ary[0] not in S2Q or S2Q[Ary[0]] not in NMapping: continue
            if S2Q[Ary[0]] not in Annots:
                Annots[S2Q[Ary[0]]] = []
                Annots[S2Q[Ary[0]]].append(Ary[1])
                Annots[S2Q[Ary[0]]].append(Ary[2].replace(' ', ';;'))
                Annots[S2Q[Ary[0]]].append(Ary[3])
                Annots[S2Q[Ary[0]]].append(Ary[4].replace(' ', ';;'))
            else:
                Annots[S2Q[Ary[0]]][0] += f';;{Ary[1]}' # pathways: delim ';;'
                Annots[S2Q[Ary[0]]][1] += f';;{Ary[2].replace(" ", ";;")}'  # GO terms: delim ';;'
                Annots[S2Q[Ary[0]]][2] += f';;{Ary[3]}' # SL: delim ';;'
                Annots[S2Q[Ary[0]]][3] += f';;{Ary[4].replace(" ", ";;")}'  # domain: delim ';;'
    return Annots

def ext_selected_fasta(uniids: Sequence, FName: str):
    seq_fasta = []
    with open(FName, 'rt') as f1:
        proteome = SeqIO.parse(f1, "fasta")
        for record in proteome:
            rid = record.id.split('|')[1] if '|' in record.id else record.id
            if rid in uniids:
                seq_fasta.append(record.format("fasta"))
    
    with open(Path(FName).stem + '_sel.fasta', 'wt') as fOut:
        fOut.write("".join(seq_fasta))

if __name__ == '__main__':
    #print(read_list('SpNames_ex.txt'))
    ext_selected_fasta(['A0A078FUP6', 'A0A078IFS9'], sys.argv[1])
