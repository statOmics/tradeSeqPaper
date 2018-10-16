import pandas as pd
import numpy as np
from GPfates import GPfates


### on simulated data
logexp = pd.read_table('/Users/koenvandenberge/logCpmDataset1.txt', index_col=0, sep=' ')
sInfo = pd.read_table('/Users/koenvandenberge/sampleInfoDataset1.txt', index_col=0, sep=' ')

m = GPfates.GPfates(sInfo, logexp)
# perform GPLVM
m.dimensionality_reduction()
# store reduced dimensions
m.store_dr()
# get pseudotime
m.infer_pseudotime(s_columns=['bgplvm_0', 'bgplvm_1'])
# model trajectory
m.model_fates(maxiter=3000)
# plot
m.make_fates_viz()
m.fates_viz.plot()
GPfates.plt.show()

# identify bifurcation
p = m.identify_bifurcation_point()
print(p)
m.s.pseudotime.min()


m.calculate_bifurcation_statistics()
