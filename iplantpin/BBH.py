from __future__ import annotations
from typing import List, Union, Tuple, Set, Dict
import sys
import subprocess
from Bio import SeqIO
import pandas as pd
from io import StringIO
from pathlib import Path

class BBH:
    def __init__(self, # Required directories
                   fbhitdir: str = 'bHits/fbhitdir',          # forward blast hits output directory
                   fbhitseqdir: str = 'bHits/fbhitseqdir',    # forward blast hits extracted subject sequences output directory
                   bbhitdir: str = 'bHits/bbhitdir',          # backward blast hits output directory
                   
                   # Diamond blastp parameters 
                   evalue: float = '1e-5',
                   pident: int = '40',
                   qcovs: int = '80',
                   num_threads: int = '28',
                   sensitivity: str = '--very-sensitive') -> None:
        self.fbhitdir = fbhitdir
        self.fbhitseqdir = fbhitseqdir
        self.bbhitdir = bbhitdir
        self.evalue = evalue
        self.pident = pident
        self.qcovs = qcovs
        self.num_threads = num_threads
        self.sensitivity = sensitivity
        
    def create_directories(self) -> None:
        """Make the required directory structure to store the intermediate results"""
        Path(self.fbhitdir).mkdir(exist_ok = True, parents = True)
        Path(self.fbhitseqdir).mkdir(exist_ok = True, parents = True)
        Path(self.bbhitdir).mkdir(exist_ok = True, parents = True)

    def exec_diamond(self, query_seq: str, subject_bDB: str,
      bOutDir: str) -> List[Union[str, float]]:
        r"""query query_sequences file against subject_sequences database using
        the diamond blastp program to identify homologous sequences in query
        """
        # the names of output files for database and blastp output
        subject_name: str = Path(subject_bDB).stem
        query_name: str = Path(query_seq).stem
        
        Bout = subprocess.check_output(['bin/diamond', 'blastp', '--query', query_seq, '--db', subject_bDB, '--evalue', self.evalue, '--threads', self.num_threads, '--max-target-seqs', '1', '--query-cover', self.qcovs, '--id', self.pident, self.sensitivity, "--outfmt", '6', 'qseqid', 'sseqid', 'pident', 'bitscore'])

        Bout = pd.read_csv(StringIO(str(Bout, 'utf-8')),    # Byte to string conversion
                             header = None, sep = "\t")
        # Column names for diamond blastp output
        Bout.columns = ['qseqid', 'sseqid', 'pident', 'bitscore']
        
        # Remove duplicates
        Bout = Bout.drop_duplicates(subset = 'sseqid')
        Bout.to_csv(bOutDir + '/' + query_name + '_' + subject_name + '.out',
                     sep = "\t", index = False)
        
        return Bout
        
    def extract_fasta(self, blastout: List[Union[str, float]], subject_fa: str) -> Dict[str, List[float]]:
        r"""Extract fasta sequences corresponding to the best hits in the 
        forward blast that will be used as query in the reverse blast
        """
        
        bHits: Set[Tuple(str, str)] = set()    
        bHits_fasta: List[str] = []

        # Extract fasta sequences for best forward blast hits
        try:
            handle = open(subject_fa, "r")
            subject_proteome = SeqIO.parse(handle, "fasta")
            for record in subject_proteome:
                if record.id in blastout.sseqid.values:
                    bHits_fasta.append(record.format("fasta"))
            handle.close()
        except IOError:
            print("Failed to open " + subject_fa)
            sys.exit(1)

        # Write a FASTA formated output file
        try:
            writeFile = open(self.fbhitseqdir + '/' + Path(subject_fa).stem + '.fasta', "w")
            writeFile.write("".join(bHits_fasta))
            writeFile.close()
        except IOError:
            print("Failed to create query.fasta")
            sys.exit(1)

    def reciprocal_blast(self, query_seq: str, subject_seq: str, query_bDB: str, 
      subject_bDB: str) -> Tuple[Tuple(str, str), ...]:
        """Perfomrs blastp on the best hits identified using forward blast and find
        best hits by querying subject sequences against query sequence database. Top
        hits common to both forward and reverse blastp ouputs will considered as
        one-to-one orthologs"""
      
        fblast_out = self.exec_diamond(query_seq, subject_bDB, self.fbhitdir)
        bfbHits = self.extract_fasta(fblast_out, subject_seq)
        bblast_out = self.exec_diamond(self.fbhitseqdir + '/' + Path(subject_seq).stem + '.fasta', 
          query_bDB,self.bbhitdir)
        
        return(set(zip(fblast_out.qseqid, fblast_out.sseqid)) & set(zip(bblast_out.sseqid, bblast_out.qseqid)))
