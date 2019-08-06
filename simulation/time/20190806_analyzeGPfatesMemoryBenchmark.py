import pandas as pd
import numpy as np
from GPfates import GPfates
import pickle as pickle
from GPfates.gp_utils import bifurcation_statistics
from memory_profiler import profile


with open('m.pkl', 'rb') as input:
    m = pickle.load(input)

## bifurcation analysis
@profile
def bifStats_kvdb(m, fate_model, e):
    p = m.identify_bifurcation_point()
    bif_stats = bifurcation_statistics(fate_model, e)
    return bif_stats

bifStats_kvdb(m, m.fate_model, m.e)


