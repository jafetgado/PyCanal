"""
pycanal: A Python package for conservation analysis of homologous proteins
"""




import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pycanal.utils as utils
import warnings

warnings.filterwarnings('ignore')
pd.set_option('use_inf_as_na', True)  #Treat inf as NaN 




class Canal():
    '''
    Class for conservation analysis of amino acid sites in aligned protein sequences.
    
    Parameters
    --------------
    fastafile : str
        Text file containing aligned protein sequences for conservation analysis in fasta
        format. Sequences must be aligned before analysis so that all sequences have the 
        same number of sites (including gaps). Unless otherwise specified, the first 
        sequence is used as the reference sequence.
    ref : int (default=0)
        Position of reference sequence in the mulitple sequence alignment file. The first
        sequence is used as reference if ref=0.
    startcount : int (default=1)
        The position numbering of the first residue in the reference sequence. Default is
        1, meaning that the first residue is labeled as position 1.
    verbose : bool (default=True)
        If True, print out analyses details.
        
    
    Example
    ----------
    Import Canal
    
    >>> from pycanal import Canal
    
    Create an instance of the Canal class
    
	>>> canal = Canal(fastafile='alignment.fasta',ref=0, startcount=-16, verbose=True)
    
    Compute conservation scores of each site in reference sequence with relative entropy 
    method
    
	>>> cons_scores = canal.analysis(include=None, method='relative')
    
    Plot the distribution of amino acids in at position 77 and save image as 
    position77.png
    
	>>> canal.plotSiteDistribution(site=77, saveplot='position77')
    
    Determine consensus sequence from the alignment
    
	>>> consensus_sequence = canal.getConsensusSequence(savefasta='consensus_sequence.fasta')
	
    '''
    
    
    
    
    def __init__(self, fastafile, ref=0, startcount=1, verbose=True):
        
        # Ensure that fastafile path is correct
        if not os.path.isfile(fastafile):
            raise OSError(f'{fastafile} not found. Specify the correct path to the ' \
                          'aligned sequences file.')
        
        # Read protein sequences and descriptions (heads)
        (headers, sequences) = utils.read_fasta(fastafile)
        
        # Ensure equal number of sequences and descriptions
        head_count, seq_count = len(headers), len(sequences)
        assert head_count == seq_count, f'There are {seq_count} sequences in {fastafile}'\
                                        f' but {head_count} descriptions. There must be '\
                                        f'an equal number of sequences and descriptions.'
        
        # Ensure equal lengths of sequences (i.e. they are aligned)
        seq_lengths = list(set([len(seq) for seq in sequences]))
        assert len(seq_lengths) == 1, f'Sequences are of varying length. Ensure that seq'\
                                      f'uences in fasta file have been aligned so that'\
                                      f' they are of the same length.'
        
        # Summary of sequence data
        reference_sequence = sequences[ref].replace('-', '')
        num_sites = len(reference_sequence)
        indices = np.arange(5) - min(startcount, 1) + 1
        startsites = np.array(list(reference_sequence))[indices]
        startnums = indices + startcount
        startlabel = [f'{startsites[i]}-{startnums[i]}' for i in range(len(startnums))]
        startlabel = ', '.join(startlabel)
        if verbose:
            print(f'\n{len(sequences)} sequences in fasta file')
            print(f'\nMultiple sequence alignment has {len(sequences[ref])} total '\
                   'positions')
            print(f'\nReference sequence is {headers[ref]}')
            print(f'\nReference sequence has {num_sites} residues (sites) and the first'\
                  f' residue is labeled as position {startcount}')
            print(f'\nReference sequence residues are labeled as {startlabel}')
        
        # Initialize
        self.fastafile = fastafile
        self.ref = ref
        self.startcount = startcount
        self.verbose = verbose
        self.reference_sequence = reference_sequence
        self.reference_header = headers[ref]
        self.site_freqs = self.msa_freqs = self.alignment_df = None
        
        
        
        
    def calcFrequencies(self, include=None):
        '''
        Calculate the frequencies of the amino acids in each site of the reference 
        sequence and in the entire alignment (i.e, all sites). By default, only the 
        frequencies of the 20 canonical amino acids are calculated, i.e A, C, D, E, F, G, 
        H, I, K, L, M, N, P, Q, R, S, T, V, W, Y. Use the include option to specify 
        non-canonical residues or gaps (i.e., B, J, O, X, Z, -) if you 
        want them to be included in the analysis.
                
        Parameters
        -------------
        include : list or None (default=None)
            List of characters to include in analysis. Ignored if `None`.
        
        Returns
        ---------
        (site_freqs, msa_freqs) : tuple
            A tuple of two Pandas dataframes containing the frequencies of amino acids in
            each site and in the entire alignment, respectively.
        site_freqs : Pandas dataframe
            Indices are amino acids, columns are sites in reference sequence.
        msa_freqs : Pandas dataframe
            Indices are amino acids, column is the frequency.
        '''
        
        # Print progress
        if self.verbose:
            print('\nCalculating amino acid frequencies...')
        
        # Amino acids to calculate frequencies
        self.aminoacid_letters = list('ACDEFGHIKLMNPQRSTVWY')
        if include is not None:
            self.aminoacid_letters.extend(include)
        
        # Read alignment as dataframe
        fasta_df = utils.read_fasta_as_df(self.fastafile)
        
        # Only use sites in reference sequence
        refseq_residue_positions = [i for i in range(fasta_df.shape[1]) \
                                    if fasta_df.iloc[self.ref,i].isalpha()]
        fasta_df = fasta_df.iloc[:, refseq_residue_positions]
        
        # Renumber sites using reference sequence numbering
        fasta_df.columns = np.arange(self.startcount, self.startcount + fasta_df.shape[1])
        
        # Calculate site frequencies
        site_counts = pd.DataFrame(index=self.aminoacid_letters, dtype=float)
        for col in fasta_df.columns:
            counts = np.unique(fasta_df[col], return_counts=True)
            counts = pd.Series(counts[1], index=counts[0])
            site_counts[col] = counts.reindex(np.array(self.aminoacid_letters))
        site_counts = site_counts.fillna(0)
        site_freqs = site_counts / site_counts.sum(axis=0)
        #self.site_frequencies = site_freqs
        
        # Calculate MSA frequencies
        all_sites = fasta_df.values.reshape(-1)
        msa_counts = np.unique(all_sites, return_counts=True)
        msa_counts = pd.Series(msa_counts[1], index=msa_counts[0])
        msa_freqs = pd.DataFrame(index=self.aminoacid_letters, dtype=float)
        msa_freqs['msa_freqs'] = msa_counts[np.array(self.aminoacid_letters)]
        msa_freqs['msa_freqs'] = msa_freqs['msa_freqs'] / msa_freqs['msa_freqs'].sum()
        
        if self.verbose:
            print('Done.')
        
        self.alignment_df = fasta_df
        self.site_freqs = site_freqs
        self.msa_freqs = msa_freqs
        
        return site_freqs, msa_freqs
    
    
    
        
    def calcConservationScores(self, site_freqs, msa_freqs, method='relative_entropy'):
        '''
        Calculate conservation scores from the amino acid distribution at each site in the
        reference sequence.
        
        Parameters
        ------------
        site_freqs : Pandas dataframe
            Pandas dataframe containing the frequencies of amino acids in each site of
            the reference sequence. The indices of the dataframe are the amino acids and
            the columns are the sites.
        msa_freqs : Pandas dataframe
            A single-column Pandas dataframe containing the frequencies of amino acids 
            in all sites (i.e, the entire alignment). The indices are the amino acids and
            the column is frequency.
        method : str {'shannon', 'relative', 'lockless'}
            Method for calculating conservation scores. 
            
            If 'shannon', calculate the Shannon entropy. Higher values indicate lower 
            conservation and greater variability at the site.
            
            If 'relative', calculate the relative entropy i.e., the Kullback-Leibler 
            divergence. Higher values indicate greater conservation and lower variability
            at the site.
            
            If 'lockless', calculate the evolutionary conservation parameter defined by 
            Lockless and Ranganathan (1999). Higher values indicate greater conservation
            and lower variability at the site.
        
        Returns
        --------
        cons_scores : Pandas dataframe
            A Pandas dataframe of conservation scores. Indices are the positions in the 
            reference sequence.
        
        References
        ------------
        .. [1] Capra, J.A. and Singh M. (2007). Predicting functionally important residues
           from sequence conservation.
           [2] Lockless, S.W. and Ranganathan R. (1999). Evolutionarily conserved pathways
           of energetic connectivity in protein families.
        
        '''
        
        # Print progress
        if self.verbose:
            print(f'\nCalculating conservation scores with {method} method...')
            
        # Ensure that site_freqs and msa_freqs have identical indices
        assert list(site_freqs.index) == list(msa_freqs.index), 'site_freqs and msa_freqs'\
                                                    ' must have identical indices'
        
        # Calculate conservation scores          
        cons_scores = []                                      
        prob_msa = msa_freqs.iloc[:,0]  # Probability of amino acid occurence in MSA
        for col in site_freqs.columns:
            prob_site = site_freqs[col] # Probability of amino acid occurence in site
            
            if method=='shannon':
                score = -np.sum((prob_site * np.log(prob_site)).fillna(0))
            elif method=='relative':
                score = np.sum((prob_site * np.log(prob_site / prob_msa)).fillna(0))
            elif method=='lockless':
                score = np.sqrt(np.sum(((np.log(prob_site / prob_msa)) ** 2).fillna(0)))
                
            cons_scores.append(score)
        cons_scores = pd.DataFrame(cons_scores, index=site_freqs.columns, 
                                   columns=[method])
        
        # Print progress
        if self.verbose:
            print('Done')
            
        return cons_scores
    
        
        
        
    def analysis(self, include=None, method='relative'):
        '''Carry out conservation analysis by calculating conservation scores from the
        alignment for each site in the reference sequence. By default, only the 
        frequencies of the 20 canonical amino acids are used in the analysis, i.e A, C, D,
        E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y. Use the include option to 
        specify non-canonical residues or gaps (i.e., B, J, O, X, Z, -) if you want them 
        to be included in the analysis.
        
        Parameters
        ------------
        include : list or None (default=None)
            List of characters to include in analysis. Ignored if `None`.
        method : str {'shannon', 'relative', 'lockless', 'all'}
            Method for calculating conservation scores. 
            
            If 'shannon', calculate the Shannon entropy. Higher values indicate lower 
            conservation and greater variability at the site.
            
            If 'relative', calculate the relative entropy i.e., the Kullback-Leibler 
            divergence. Higher values indicate greater conservation and lower variability
            at the site.
            
            If 'lockless', calculate the evolutionary conservation parameter defined by 
            Lockless and Ranganathan (1999). Higher values indicate greater conservation
            and lower variability at the site.
            
            If 'all', calculate all the above mentioned conservation scores.
            
            
        Returns
        ---------
        cons_scores : Pandas dataframe
            A Pandas dataframe of conservation scores. Indices are the positions in the 
            reference sequence. Columns are the conservation analysis methods.
        '''
        
        # Calculate amino acid frequencies
        site_freqs, msa_freqs = self.calcFrequencies(include=include)
        self.site_freqs = site_freqs
        self.msa_freqs = msa_freqs
        
        # Calculate conservation scores
        if method != 'all':
            cons_scores = self.calcConservationScores(site_freqs, msa_freqs, 
                                                      method=method)
        else:
            cons_scores = pd.DataFrame(index=site_freqs.columns, dtype=float)
            for each_method in ['shannon', 'relative', 'lockless']:
                each_cons_scores = self.calcConservationScores(site_freqs, msa_freqs, 
                                                               method=each_method)
                cons_scores[each_method] = each_cons_scores.iloc[:,0]
        
        return cons_scores
                
        
    
    
    def plotSiteDistribution(self, site, saveplot=None):
        '''Plot the amino acid distribution at a specific site in the alignment. The 
        amino acid at that site in the reference sequence is colored red in the plot.
        All other amino acids are colored gray.
        
        Parameters
        -------------
        site : int
            The position in the reference sequence to be plotted
        saveplot : str or None (default=None)
            If not None, save the plot with the name specified by saveplot.
        
        Examples
        -----------
        
        Plot the distribution of amino acids at site 100 and save the plot as
        'position100.png'
        
        >>> canal.plotSiteDistribution(site=100, saveplot='position100')
        
        '''
        
        # Ensure that site_freqs is not None
        if (self.site_freqs is None) or (self.msa_freqs is None):
            raise NotImplementedError("The amino acid frequencies have not been "\
                                      "calculated because 'calcFrequencies' or 'analysis'"\
                                      " methods have not been implemented.")
        
        # Get non-zero amino acid percentages in site
        this_site_dist = self.site_freqs.loc[:,site] * 100
        this_site_dist = this_site_dist.iloc[(this_site_dist > 0).values]
        this_site_dist = this_site_dist.sort_values(ascending=False)
        
        # Plotting parameters 
        plt.rcParams['figure.figsize'] = [6,4]
        axis_font = {'fontname':'Arial', 'size':'18'}
        ticks_font = {'fontname':'Arial', 'size':'9'}
        xvalues = np.arange(this_site_dist.shape[0])
        
        # Color the amino acid in reference sequence red, and all others grey
        colors = ['dimgrey'] * len(xvalues)
        refseq_residue  = self.reference_sequence[site - self.startcount]
        refseq_residue_pos = list(this_site_dist.index).index(refseq_residue)
        colors[refseq_residue_pos] = 'firebrick'
        
        # Plot distribution
        plt.bar(xvalues, this_site_dist.values, width=0.75, linewidth=0.75, color=colors, 
                edgecolor='black')
        _ = plt.xticks(xvalues, this_site_dist.index, **ticks_font)
        _ = plt.yticks(**ticks_font)
        plt.ylabel('Frequency (%)', **axis_font)
        plt.title(f'Position {site}', **axis_font)
        plt.tight_layout()
        
        # Save plot
        if saveplot is not None:
            plt.savefig(f'{saveplot}.jpg', format='jpg', dpi=600)
        plt.show(); plt.close()
        
      
        
      
    def getConsensusSequence(self, savefasta=None):
        '''Obtain the consensus protein sequence from the alignment. The amino acid at
        each site in the consensus sequence is the majority amino acid at that site in 
        the alignment. Only sites without gaps in the reference sequence are considered
        
        Parameters
        ------------
        savefasta : str or None
            If str, write the reference and consensus sequences to a savefasta. Ignored
            if None.
        
        Returns
        ---------
        consensus_sequence : str
            Consensus sequence from the alignment.
        
        '''
        
        # Ensure that site_freqs is not None
        if self.alignment_df is None:
            raise NotImplementedError("Cannot determine consensus sequence since neither"\
                                      " 'calcFrequencies' nor 'analysis' methods have"\
                                      " been implemented.")
        
        # Get consensus sequence
        consensus_sequence = ''
        for col in self.site_freqs:
            majority_amino_acid = self.site_freqs[col].sort_values(ascending=False).index[0]
            consensus_sequence += majority_amino_acid
        
        if savefasta is not None:
            headers = ['consensus_sequence', self.reference_header]
            sequences = [consensus_sequence, self.reference_sequence]
            utils.write_fasta(headers, sequences, savefasta)
        
        return consensus_sequence
        

    
    
  
