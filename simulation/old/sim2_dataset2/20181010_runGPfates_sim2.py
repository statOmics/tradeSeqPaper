import pandas as pd
import numpy as np
from GPfates import GPfates

logexp = pd.read_table('/Users/koenvandenberge/logCpmDataset2.txt', index_col=0, sep=' ')
sInfo = pd.read_table('/Users/koenvandenberge/sampleInfoDataset2.txt', index_col=0, sep=' ')

m = GPfates.GPfates(sInfo, logexp, pseudotime_column=4)
# perform GPLVM
m.dimensionality_reduction()
# store reduced dimensions
m.store_dr()
m.model_fates(t='pseudotime')
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
np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/simulation/sim2_dataset2/GPfatesWeights_sim2Dataset2.txt', weights)
np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/simulation/sim2_dataset2/GPfatesBifStats_sim2Dataset2.txt', bif_stats)

## GPfates without true pseudotime
m = GPfates.GPfates(sInfo, logexp)
m.dimensionality_reduction()
m.store_dr()
m.infer_pseudotime()
m.model_fates()
m.make_fates_viz()
m.fates_viz.plot()
GPfates.plt.show()
p = m.identify_bifurcation_point()

# get output
weights = m.fate_model.phi
from GPfates.gp_utils import bifurcation_statistics
bif_stats = bifurcation_statistics(m.fate_model, m.e)

# write output
np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/simulation/sim2_dataset2/GPfatesWeights_sim2Dataset2_noTruePseudoT.txt', weights)
np.savetxt('/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/simulation/sim2_dataset2/GPfatesBifStats_sim2Dataset2_noTruePseudoT.txt', bif_stats)