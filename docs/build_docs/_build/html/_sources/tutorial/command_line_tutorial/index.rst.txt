From Command Line
==========================

First install the software, export SCHORDINGER env variable and download the example folder:

- conda install -c NostrumBioDiscovery -n growai python=3.7 growai

- export SCHRODINGER=/opt/schrodinger2018-1/

- git clone https://github.com/danielSoler93/growai.git

- cd growai/growai/examples/toluene


Growing methodology
------------------------

- python -m growai.grow --pdb toluene.pdb --resname TOL --only_grow --grow_iterations 8 


Growing and ranking methodology (commercial schrodinger is need it)
---------------------------------------------------------------------------

- python -m growai.grow --pdb toluene.pdb --resname TOL --grow_iterations 8
