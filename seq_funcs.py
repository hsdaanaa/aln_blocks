#TODO: cut unnecessary codes
import sys, os
import pandas as pd
import numpy as np
from utils import map_dict_vals
#------------------------------------------------------------------
def has_N_or_gap(seq): 
    if 'N' in seq:
        return True
    if '-' in seq:
        return True
    
    return False
#------------------------------------------------------------------
def aln_to_df(path, site_type = 'codon'):
    """reads a fasta formatted alignment file into a dataframe. 

    flow
    ----
    - checks specification for seq column. if site_col_sep is 
     False, then the entire sequences (for each entry) sequence 
     appears in a single column, if True, then each site (nucleotide
     or codon; appears in separate columns. this is specified by site_type) 
    -reads FASTA-formatted cds/aln file into a dictionary 
    - checks input site_type and uses map_dict_vals to 
      convert sequences into the desired sequence format
    - converts cds/aln dictionary into dataframe.
    
    parameters
    ----------
    path: str
        path to FASTA-formatted alignment/cds file
    site_col_sep: bool
        specifies how to present seqinfo
        
    site_type: 'codon', 'nucleotide'
        specifies format of sequences in the alignment
        if 'codons' each column/row is a nucleotide 
        triplet. if 'nc', each column/row is a single
        nucleotide. 

    returns
    -------
    DataFrame"""
    
   # if site_col_sep == False:
   #    return fasta_to_df(path, columns = ['ID', 'seq'])

    # converts alignment/cds into a dictionary
    cds_dict = scaffold_data_to_dict(path) #fasta_to_dict(path)

    # checks value of site_type
    if site_type == 'nucleotide':
        # converts each dictionary value into a list of strings (single nucleotides)
        cds_dict = map_dict_vals(cds_dict, list) 

    elif site_type =='codon':
        # converts each dictionary value into a list of string triplets (codons)
        cds_dict = map_dict_vals(cds_dict, get_codons)

    # converts dictionary to a dataframe and outputs dataframe
    df = pd.DataFrame.from_dict(cds_dict, orient = 'index') 

    return df 
#------------------------------------------------------------------
def df_to_fasta(df, path_to_outputfile_dir, filename, cds_name = 'columns', seq_name_col = 0):
    """takes an input alignment dataframe and converts it 
    into a fasta file. path_to_outputfile_dir and filename specify 
    the dir path and name of the fasta file respectively. 
    cds_name specifies the position (row or column) of the
    sample_names in the of the input dataframe. 
    seq_name_col specifies the index position of the first sample name
    in the alignment
    
    flow
    ----
    - converts alignment into a dictionary 
        - the sequence columns (or rows) of each sample in 
        the alignment are joined to form a single continous 
        nucleotide sequence. the sample name and the sequence
        are stored in a dictionary (key = sample name, 
        value = sequence)
        
    - converts dictionary into FASTA-file using dict_to_fasta()
    
    parameters 
    ----------
    df: pandas.core.DataFrame
        alignment dataframe 
    path_to_outputfile_dir: str
        path to directory of output FASTA file
    filename: str
        name of output FASTA file
    cds_name: 'columns', 'rows'
        specifies position of sample names in the alignment dataframe
        'columns' if sample names are dataframe columns. 'rows' if
        sample names are dataframe rows. 
    seq_name_col: int
        index position of first sample name, (if samples are columns)
    
    returns
    -------
    None"""
    
    fasta_dict = {} # dictionary to store sample names and sequences
    
    # loops over columns/rows, concatenates single nucleotides
    # to form full sequence. assigns sample name as key and sequence
    # as value. 
    if cds_name == 'columns':
        for column in df.columns:
            seq = ''.join(list(df[column]))
                      
            fasta_dict[column] = seq # assign column name as dict key 
        
        dict_to_fasta(fasta_dict, path_to_outputfile_dir, filename)
        
    elif cds_name == 'row' or cds_name == 'index':
        df = df.reset_index(False, inplace = False)
        for row in df.values:
            cds_name = row[seq_name_col]
            seq = ''.join(row[1:])
            fasta_dict[cds_name] = seq # assign column name as dict key 
        
        dict_to_fasta(fasta_dict, path_to_outputfile_dir, filename) 
#------------------------------------------------------------------------
def dict_to_fasta(dict_var, path_to_output_file, output_file_name):
    """outputs a dictionary to a fasta formatted file. cds names should be
    keys and values should be sequences"""

    # make file name
    file_name = os.path.join(path_to_output_file, output_file_name)
    # open file 
    file_obj = open(file_name, 'w')
     # loop over dictionary keys and writes cds_name/sequence lines.
    for cds_name in dict_var.keys():
            seq = ''.join(dict_var[cds_name]) # gets cds_sequence from dictionary of cds
            # writes data into file
            cds_name = cds_name.rstrip('\n')
            file_obj.write('>{}\n{}\n'.format(cds_name, seq))
            
    file_obj.close() 
#------------------------------------------------------------------------   
def scaffold_data_to_dict(input_path, delim1 ='>', delim2 = 'lastpos', verbose = 0): # change function name
    """"""
    descs = []
    seqs = []
    temp_line = []
    with open(input_path, 'r') as file:
        
            for line in file:
                if line.startswith(delim1):
                    descs.append(line.lstrip('>').rstrip('\n'))
                    if ''.join(temp_line) == '':
                        continue
                    else:
                        seqs.append(''.join(temp_line))
                    
                    temp_line = []
                else:
                    temp_line.append(line.replace('\n', ''))
                

            seqs.append(''.join(temp_line))
            
            return dict(zip(descs, seqs))
#------------------------------------------------------------------------
def fasta_to_dict(input_path):
    """takes an input path to a fasta file and delimiters between 
    an identifier, creates a dictionary with identifier: sequence"""
    
    cds_seq_dict = {}
    cds_names    = []
    with open(input_path, 'r') as file:
        for line in file:
            if line.startswith(">"):
                cds_name = line[1:]
                cds_name = cds_name.rstrip('\n')
                key = line[1:].rstrip('\n')
                cds_seq_dict[cds_name] = ''
                cds_names.append(cds_name)
            else:
                cds_seq_dict[cds_name] += line.rstrip('\n')
            
    assert len(cds_names) == len(cds_seq_dict), 'some sequences in the alignment had the same name.sequence names in alignment {}'.format(cds_names)
    return cds_seq_dict
#------------------------------------------------------------------------
def get_codons(sequence):
    """Returns a codon sequence of a given DNA/RNA seq input.
     A codon is three nucleotides long; so the function 
     groups seq input into threes
    
    Parameters
    ----------
    sequence : str
    
                nucelotide sequence
    Returns
    -------
    list"""
    
    codons = []
    
    assert len(sequence) %3 == 0, 'input sequence length was not a multiple of 3. got: {}'.format(len(sequence))
    
    for index in range(0,len(sequence),(3)):
        codons.append(sequence[index:index+3])

    return codons    