import pandas as pd
import numpy as np
from GPfates import GPfates
import pickle as pickle

def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


logexp = pd.read_table('./timeBenchLogCpm.txt', index_col=0, sep=' ')
sInfo = pd.read_table('./timeBenchSampleInfo.txt', index_col=0, sep=' ')


# run GPfates by feeding it the true pseudotime.
m = GPfates.GPfates(sInfo, logexp, pseudotime_column=1)
# perform GPLVM
m.dimensionality_reduction()
# store reduced dimensions
m.store_dr()
m.model_fates(t='global_pseudotime')

save_object(m, "m.pkl")


## identify bifurcation
#p = m.identify_bifurcation_point()
##print(p)
#
## get output
#weights = m.fate_model.phi
#from GPfates.gp_utils import bifurcation_statistics
#bif_stats = bifurcation_statistics(m.fate_model, m.e)
#
## write output
#np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/GPfatesWeights.txt', weights)
#np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/GPfatesBifStats.txt', bif_stats)
