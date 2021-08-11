import sys
import pandas as pd
from sklearn.neighbors import LocalOutlierFactor
import argparse

################################################################################
# Input
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--control", default=".DAJIN_temp/classif/tmp_control")
parser.add_argument("-s", "--sample", default=".DAJIN_temp/classif/tmp_control")
parser.add_argument("-t", "--threads", default="1")
args = parser.parse_args()

df_control = pd.read_csv(args.control, header=None)
df_sample = pd.read_csv(args.sample, header=None)
threads = int(args.threads)

################################################################################
# LOF
################################################################################

lof = LocalOutlierFactor(
    n_neighbors=20,
    algorithm="auto",
    leaf_size=30,
    metric="euclidean",
    contamination="auto",
    novelty=True,
    n_jobs=threads,
)

lof.fit(df_control)

output = pd.Series(lof.predict(df_sample))

################################################################################
# Output
################################################################################

output.to_csv(sys.stdout, index=False, header=False)
