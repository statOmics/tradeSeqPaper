import pandas as pd
import numpy as np
from GPfates import GPfates

logexp = pd.read_table('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_multifurcating_4/simDyntoyLogCpm.txt', index_col=0, sep=' ')
sInfo = pd.read_table('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_multifurcating_4/sampleInfoSimDyntoy.txt', index_col=0, sep=' ')

# run GPfates by feeding it the true pseudotime.
m = GPfates.GPfates(sInfo, logexp, pseudotime_column=1)
# perform GPLVM
m.dimensionality_reduction()
# store reduced dimensions
m.store_dr()
m.model_fates(t='global_pseudotime')
# plot
m.make_fates_viz()
m.fates_viz.plot()
GPfates.plt.show()

# identify bifurcation
p = m.identify_bifurcation_point()
print(p)

# get output
weights = m.fate_model.phi
from GPfates.gp_utils import bifurcation_statistics
bif_stats = bifurcation_statistics(m.fate_model, m.e)

# write output
np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_multifurcating_4/GPfatesWeights.txt', weights)
np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_multifurcating_4/GPfatesBifStats.txt', bif_stats)

### GPfates without true pseudotime
#m = GPfates.GPfates(sInfo, logexp)
#m.dimensionality_reduction()
#m.store_dr()
#m.infer_pseudotime()
#m.model_fates()
#m.make_fates_viz()
#m.fates_viz.plot()
#GPfates.plt.show()
#p = m.identify_bifurcation_point()
#
## get output
#weights = m.fate_model.phi
#from GPfates.gp_utils import bifurcation_statistics
#bif_stats = bifurcation_statistics(m.fate_model, m.e)
#
## write output
#np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/simulation/sim2_dyntoy_multifurcating_4/GPfatesWeights_sim2Dataset2_noTruePseudoT.txt', weights)
#np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/simulation/sim2_dyntoy_multifurcating_4/GPfatesBifStats_sim2Dataset2_noTruePseudoT.txt', bif_stats)
