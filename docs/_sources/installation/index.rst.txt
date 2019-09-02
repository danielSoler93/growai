Installation
====================

From Conda (recommended)
--------------------------

::

 conda create -c NostrumBioDiscovery -n grow_test python=3.7  growai
 source activate grow_test
 export SCHRODINGER=/opt/schrodinger2018-1/
 python -m growai.gorw -h


From PyPi
-----------

::

 pip install growai

 conda install rdkit (or build it from source as it has no pip support)

 
From Source Code
---------------------

::

 git clone https://github.com/danielSoler93/growai.git
 
 cd growai

 python setup.py install


