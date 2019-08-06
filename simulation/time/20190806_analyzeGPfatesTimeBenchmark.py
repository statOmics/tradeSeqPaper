import pandas as pd
import numpy as np
import sys
sys.path.append("/accounts/campus/hector.rouxdebezieux/R/x86_64-pc-linux-gnu-library/3.5/GPfates/venv/lib/python3.7/site-packages")
from GPfates import GPfates
import pickle as pickle
from GPfates.gp_utils import bifurcation_statistics
from memory_profiler import profile


with open('m.pkl', 'rb') as input:
    m = pickle.load(input)

p = m.identify_bifurcation_point()
bif_stats = bifurcation_statistics(m.fate_model, m.e)
