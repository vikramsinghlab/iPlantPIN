from .BBH import BBH
import sys
from .utils import progress_bar, orca_input, write_dict, read_annots, ext_selected_fasta, read_list

def _str2uni(fName, S2Q):
    S2U = dict()
    with open(fName, 'rt') as f1:
        for line in f1:
            Ary = line.strip().split(',')
            if Ary[0] not in S2Q: continue
            S2U[Ary[1]] = Ary[0]
    
    diff = set(S2Q.keys()).difference(S2U.values())
    if diff:
        for i in diff: 
            print(i)
            del S2Q[i]
        
    return S2U, S2Q
    
def interolog_bbh(QuerySeq: str, qbDb: str):
    """Read the BBH output and construct an interologous network using STRING
    data"""
    
    # Create required directories and objects
    bbhobj = BBH()
    bbhobj.create_directories()
    
    SpNames = read_list('Data/SpNames_ex.txt')
    TEMPLATES = len(SpNames)
        
    sys.stdout.write('# Step 1: Identifying orthologues in all templates\n')
    # Step 1: Identify Subject to query mapping for all the targets
    Subj2Query = dict()
    for index, fName in enumerate(SpNames):
        orthologs = bbhobj.reciprocal_blast(QuerySeq, '../Examples_filter/' + fName + '_ex.fasta', qbDb, f"../sbDb/{fName}/{fName}.dmnd")
        
        for qid, sid in orthologs:
            if '|' in qid: qid = qid.split('|')[1]
            if '|' in sid: sid = sid.split('|')[1]
            Subj2Query[sid] = qid
        progress_bar(index + 1, TEMPLATES)

    sys.stdout.write('\nDone #### mapping orthologues\n')

#    write_dict(Subj2Query, 'Subj2Query.txt')

    # StepII: Read the String to Uniprot mapping
    Str2Uni, Subj2Query = _str2uni('../Examples_filter/Uni2Str_ex.txt', Subj2Query)

    sys.stdout.write('# Step    2: Identifying orthologus pairs in all templates\n')
    Network = dict(); Query2Subj = dict()
    # StepIII: Map the string ppi data of template to uniprot qury data and construct network
    for index, fName in enumerate(SpNames):
        with open('../ExamplesPIN_filter/' + fName + '_ex.txt', 'rt') as f1:
            for line in f1:
                Ary = line.strip("\n").split(',')

                if Ary[0] not in Str2Uni or Ary[1] not in Str2Uni: continue
                if Str2Uni[Ary[0]] not in Subj2Query or Str2Uni[Ary[1]] not in Subj2Query: continue

                src = Subj2Query[Str2Uni[Ary[0]]]
                trgt = Subj2Query[Str2Uni[Ary[1]]]                
                taxid = Ary[0].split(".")[0]
                
                if f"{src} {trgt}" in Network:
                    Network[f"{src} {trgt}"].append(taxid)
                elif f"{trgt} {src}" in Network:
                    Network[f"{trgt} {src}"].append(taxid)
                else:
                    Network[f"{src} {trgt}"] = []
                    Network[f"{src} {trgt}"].append(taxid)
                
        progress_bar(index + 1, TEMPLATES)
    sys.stdout.write('\nDone #### mapping orthologus pairs\n')
    return Network, Subj2Query

if __name__ == '__main__':
    Network, Subj2Query = interolog_bbh(sys.argv[1], sys.argv[2])
