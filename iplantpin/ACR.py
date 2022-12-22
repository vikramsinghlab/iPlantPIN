#-*- coding:utf-8 -*-

from typing import Iterator, Iterable, ByteString
import re, sys, os
import numpy as np
import pandas as pd
import timeit

from multiprocessing import cpu_count
from functools import partial
from joblib import Parallel, delayed
import iplantpin as ipp

def resplit(chunks_of_a_file: Iterator, split_pattern: re.Pattern) -> Iterable[ByteString]:
    """Reads chunks of a file one chunk at a time, splits them into data rows by
    `split_pattern` and joins partial data rows across chunk boundaries."""
    
    partial_line = None
    for chunk in chunks_of_a_file:
        if partial_line:
            partial_line = "".join((partial_line, chunk))
        else:
            partial_line = chunk
        if not chunk:
            break
        lines = split_pattern.split(partial_line)
        partial_line = lines.pop()
        yield from lines
    if partial_line:
        yield partial_line

def clean_fasta(FSeqIn):
    """Load fasta file and clean unkown amino acids by replacing them with -. It
    also removes sequences smaller than threshold value (lag + 1)"""
    
    SeqDict = dict()
    with open(FSeqIn, mode="r") as file_descriptor:
        buffer_size = 1024 * 8
        sentinel = ""
        chunks = iter(partial(file_descriptor.read, buffer_size), sentinel)
        data_rows_delimiter = re.compile(">")   
        lines = resplit(chunks, data_rows_delimiter)
        next(lines)
        
        for line in lines:
            Ary = line.strip().split("\n")
            ind = 0
            if '|' in Ary[0]: ind = 1
            name, sequence = Ary[0].split('|')[ind], re.sub('[^ARNDCQEGHILKMFPSTWYV]', '', ''.join(Ary[1:]).upper())
            if len(sequence) < 31: continue 
            SeqDict[name] = sequence
    
    return SeqDict

def moran_acr(sequence, AAidx, lag = 30):
    """Computes moran autocorrelation for a protein and return lag * features
    variables"""
        
    Sprop = AAidx.loc[list(sequence),:]
    xmean = Sprop.mean()
    N = len(sequence)
    
    acr_m = []
    for d in range(1, lag + 1):
        t1 = Sprop.iloc[:N - d].reset_index(drop = True) - xmean
        t2 = Sprop.iloc[d:].reset_index(drop = True) - xmean
        s1 = (t1 * t2).mean()
        s2 = (Sprop - xmean).pow(2).mean()
        s3 = s1 / s2
        acr_m.append(s3)
    
    return(pd.DataFrame(acr_m).to_numpy().flatten(order = 'F').tolist())
    
def read_ppi(file, sep = '\n'):
    """Returns a list of protein interactions"""
    
    ppi = [];
    with open(file) as f1:
        for line in f1:
            Ary = line.strip().split()
            ppi.append(Ary[0:2])
    return ppi

def acr(FSeq, FInt):
    """Computes moran autocorrelation for each protein involved in interaction
    and returns a concatinated vector of length 2 * Nprops * 30"""
    
    # Read sequences
    SeqDict = clean_fasta(FSeq)
            
    # read amino acid indices 
    AAidx = pd.read_csv('Data/AAindex.txt', header = None, index_col = 0, sep = "\t")
    AAidx.columns = ['Pb', 'Pl', 'Hb', 'SCV', 'Hl', 'SASA', 'NCI']
    
    # standardization
    AAidx = (AAidx - AAidx.mean()) / AAidx.std()
    
    sys.stdout.write('# Extracting Features\n')
    acr_out = dict(); Total = len(SeqDict)
    for index, p in enumerate(SeqDict):
        acr_out[p] = moran_acr(SeqDict[p], AAidx)
        
        ipp.progress_bar(index + 1, Total)
    
    sys.stdout.write('\n# Done Extracting Features\n')
    # read PPI pair
    plantppi = read_ppi(FInt)

    fOut = open(ipp.Path(FInt).stem + '_ACR.txt', 'w')
    for p1, p2 in plantppi:
        if len(acr_out[p1]) < 2 or len(acr_out[p2]) < 2: continue
        print(f"{p1} {p2}", * acr_out[p1] + acr_out[p2], file = fOut)
    fOut.close()

if __name__ == '__main__':
    #'Data/Ensemble_PPIs_0.99.txt' 
    acr('Data/Ensemble_Expt_0.99.fasta', 'EnsemblePlant_NPPIs_0.99.txt')
