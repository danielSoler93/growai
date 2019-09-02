import os
import sys

try:
    SCHRODINGER = os.environ["SCHRODINGER"]
except KeyError:
    sys.exit("Please set SCHRODINGER env variable to /path/schrodinger2018-1/ by doing \
    export SCHRODINGER=/path/schrodinger2018-1/")
    

if not SCHRODINGER:
    sys.exit("Please set SCHRODINGER env variable to /path/schrodinger2018-1/ by doing \
    export SCHRODINGER=/path/schrodinger2018-1/")
else:
    pass
