**PyCanal: A Python package for conservation analysis of homologous protein sequences**
===========================================================================================

PyCanal is a python tool for evaluating amino acid conservation in sites of homologous proteins. Homologous proteins may be obtained via tools such as BLAST, HMM, etc. PyCanal takes as input a multiple sequence alignment of homologous proteins in fasta format, and then calculates conservation scores for each site in a reference sequence chosen by the user from the alignment.

Installation
-------------

Download from GitHub and install:

.. code:: shell-session

    git clone https://github.com/jafetgado/pycanal.git
    cd pycanal
    python setup.py install


Examples
----------
An example showing how to compute conservation scores with PyCanal for 1,748 sequences of family 7 glycoside hydrolases.

.. code:: python

    from pycanal import Canal

    # Create an instance of the Canal class
    canal = Canal(fastafile='alignment.fasta', #Multiple sequence alignment (MSA) of homologous sequences
                  ref=0, #Position of reference sequence in MSA, use first sequence
                  startcount=-16, #Position label of first residue in reference sequence
                  verbose=True # Print out progress
                  )

    # Compute conservation scores of each site in reference sequence with relative entropy method
    cons_scores = canal.analysis(include=None, method='relative')

    # Plot the distribution of amino acids in at position 77 and save image as position77.png
    canal.plotSiteDistribution(site=77, saveplot='position77')

    # Determine consensus sequence from the alignment and save as fasta file
    consensus_sequence = canal.getConsensusSequence(savefasta='consensus_sequence.fasta')


Amino acid distribution at position 77

.. image:: https://github.com/jafetgado/pycanal/blob/master/example/position77.png

Relative entropy conservation scores

.. image:: https://github.com/jafetgado/pycanal/blob/master/example/position77.png







Citation
--------------
If you find PyCanal useful, please cite:

    Gado, J.E. (2020). Data-driven computational study of enzymes for a bio-based circular economy.

