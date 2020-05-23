import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

dir = os.environ["RAREVARDIR"]  # upper-level directory

v6_outlier_file = os.path.join(dir,"archive/data/outliers_medz_picked_counts_per_ind.txt")
v8_outlier_file = os.path.join(dir,"archive/data/v8/outliers_medz_picked_counts_per_ind.txt")
fig_file = os.path.join(dir,"paper_figures/v6_v8_outlier_hist.pdf")

# read in outlier data
v6_outlier = pd.read_csv(v6_outlier_file, sep='\t', names=["indiv", "count"])
v8_outlier = pd.read_csv(v8_outlier_file, sep='\t', names=["indiv", "count"])

# make histogram
plt.figure(figsize=(10,5))
plt.hist([v6_outlier["count"], v8_outlier["count"]], bins=range(1, 60, 10), color=['r','b'], alpha=0.8, label=["v6p", "v8"])
plt.legend()
plt.title("Comparison of Outliers between GTEx v6p and GTEx v8")
plt.xlabel("Number of multi-tissue outlier genes")
plt.ylabel("Number of individuals")
plt.savefig(fig_file, bbox_inches='tight')
plt.close()
