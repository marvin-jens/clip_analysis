from glob import glob
import os
here = os.path.abspath(os.getcwd())
for fname in glob(os.path.join(here,"*.py")):
    name = os.path.basename(fname).split('.py')[0]
    __import__(name)