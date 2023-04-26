# #!/usr/bin/python3
#--------------------------Description---------------------------------
# functions to extract and calculate statistics for blocks of aligned 
# sites in pairwise or multiple sequence alignments
#--------------------------Dependencies--------------------------------
import sys, os, numpy as np, pandas as pd
from seq_funcs import fasta_to_dict, get_codons, aln_to_df
#------------------------------Codes-----------------------------------
def get_aln_block_stats(pass): 
    pass

#----------------------------------------------------------------------
def get_aln_block_list(path_to_fasta_aln, gap_char = '-', verbose =  0): 
    """outputs blocks of aligned sites (i.e groups of aligned sites that 
    do not containin a gap in any of the sample sequences.

    for algorithm to get blocks, see get_blocks function.

    parameters
    ----------
    path_to_fasta_aln: str
     path to FASTA formatted alignment
    gap_char: str
        character that indicates gap presence
    verbose : int 0,1,2
        for debugging output

    returns
    -------
    list"""

    # read alignment as dataframe
    aln_df = aln_to_df(path_to_fasta_aln, site_type = 'nucleotide').T
    if verbose in [1,2]: 
        print('aln len: {}'.format(len(aln_df)))
        if verbose == 2: 
            print('fetching list of aligned sites')
    
    # enlist nucleotides at each aligned site
    aln_df['aln_sites'] = aln_df.apply(lambda x: tuple(x), axis = 1)
    if verbose in [1,2]:
        print('converting sites to list')

    aln_sites = aln_df['aln_sites'].to_list()
    if verbose in [1,2]:
        print('getting aligned blocks')  
    
    # get aligned blocks
    block_list = get_blocks(aln_sites, gap_char = gap_char)
    if verbose in [1,2]:
        print('block num: {}'.format(len(block_list)))

    return block_list
#------------------------------------------------------------------------
def get_blocks(aln_sites_list, gap_char): 
    """"outputs blocks of aligned sites (i.e groups of aligned sites 
    that do not containin a gap in any of the sample sequences. 

    algorithm 
    ---------
    creates list to store aligned sites for a block (block_var) and to store blocks (block_list_var)
    loops over each aligned site
        checks for gap char
            if gap char exists, 
                appends any previous info to block_var
                appends block_var to block_list_var
                restart loop starting at next column
            otherwise 
                appends aligned site to block_var
        
    parameters
    ----------
    aln_sites: list of tuples
        each tuple corresponds to sample sequences for one aligned site
        e.g. [(A, A, A)] would mean all sample sequences have 'A' at the 
        aligned site
    gap_char: str
        character that indicates a gap in the alignment
    
    returns 
    -------
    list of lists"""

    count  = 0  # to keep track of aln pos
    blocks = [] # to store all blocks
    block  = [] # to store one block

    # loop over aligned sites
    while count < len(aln_sites_list):
        # extact aln site
        aln_site = aln_sites_list[count]
        # check for gap
        if not gap_char in aln_site: 
            # if no gap, then append to block var
            block.append(aln_site)
        else:
            # otherwise, check if block is empty 
            if not block == []: 
                blocks.append(block) # append block
                block = []           # empty for next block
        count+=1                     # increase count to move to next aln site

    # after above process append remaining of alignment if (there's no gap)
    if not block == []:
        blocks.append(block)

    return blocks