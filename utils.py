"""
Utility functions for analyses
"""




import pandas as pd




def read_fasta(fasta):
    """
    Read the protein sequences in a fasta file
    
    Parameters
    -----------
    fasta : str
    	Filename of fasta file containing protein sequences
    
    Returns
    ----------
    (list_of_headers, list_of_sequences) : tuple
    	A tuple of corresponding lists of  protein descriptions and protein sequences 
        in the fasta file
    	
    """
    with open(fasta, 'r') as fastaobject:
        headers, sequences = [], []
        for line in fastaobject:
            if line.startswith('>'):
                head = line.replace('>','').strip()
                headers.append(head)
                sequences.append('')
            else :
                seq = line.strip()
                if len(seq) > 0:
                    sequences[-1] += seq
    return (headers, sequences)




def read_fasta_as_df(fasta):
    '''Read a multiple sequence alignment file in fasta format and return a Pandas 
    dataframe with the sequences as indices and the sites as columns.'''
    
    (headers, sequences) = read_fasta(fasta)
    df = pd.DataFrame([list(seq) for seq in sequences], index=headers)
    return df




def write_fasta(headers, sequences, savefasta):
    '''Write sequences to fasta file
    
    Parameters
    -----------
    headers : list
        List of headers/descriptions for sequences
    sequences : list
        List of sequences corresponding to the descriptions in headers
    savefasta : str
        Path of fasta file to write sequences to
        
    '''
    
    with open(savefasta, 'w') as fastapath:
        for i in range(len(headers)):
            fastapath.write('>' + headers[i] + '\n' + sequences[i] + '\n')
