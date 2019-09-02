From Command Line
==========================

First install the software and export SCHORDINGER env variable:

- conda install -n growai python=3.7 -c NostrumBioDiscovery growai

- export SCHRODINGER=/opt/schrodinger2018-1/

Growing methodology
------------------------

- git clone https://github.com/danielSoler93/growai.git

- cd growai/growai/examples/toluene

- python -m growai.grow --pdb toluene.pdb --resname TOL --only_grow --grow_iterations 8 
