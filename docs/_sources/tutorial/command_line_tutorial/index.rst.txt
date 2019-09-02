From Command Line
==========================

First install the software, export SCHORDINGER env variable and download the example folder:

- conda install -c NostrumBioDiscovery -n growai python=3.7 growai

- export SCHRODINGER=/opt/schrodinger2018-1/

- git clone https://github.com/danielSoler93/growai.git

- cd growai/growai/examples/toluene


Growing methodology
++++++++++++++++++++++++

::

    python -m growai.grow --pdb toluene.pdb --resname TOL --only_grow --grow_iterations 8 

Plot druglikeness
-------------------------

::

    python -m growai.analysis.analise  --sdf_files round*/round*.sdf --druglikeness --title Druglikeness --xlabel LSTM --ylabel QED --ylim 0 0.5 --output drug_like.png


Asses validity
-------------------------

::

    python -m growai.analysis.analise  --sdf_files round*/round*.sdf --validity --iterations 3 --growing_sites 4

Asses chemical space
-----------------------

::

    python -m growai.analysis.analise  --sdf_files round*/round*.sdf --pca --title PCA --xlabel PCA1 --ylabel PCA2 --ylim 0 0.5 --output pca.png


Growing and ranking methodology (commercial schrodinger is need it)
---------------------------------------------------------------------------

- python -m growai.grow --pdb toluene.pdb --resname TOL --grow_iterations 8 (Still on trial)
