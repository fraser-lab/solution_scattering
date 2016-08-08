import pandas
import numpy as np

TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
VARIANTS = ["WT", "S99T"]
# VARIANTS = ["WT"]
for variant in VARIANTS:
  protein_off = pandas.read_csv("datfiles/{}_protein_off.dat".format(variant), sep=" ", header=None, names=["q", "I", "sigI"])
  buffer_off = pandas.read_csv("datfiles/{}_buffer_off.dat".format(variant), sep=" ", header=None, names=["q", "I", "sigI"])
  subtracted = protein_off - buffer_off
  subtracted.q = protein_off.q
  subtracted.sigI = np.sqrt(protein_off.sigI**2 + buffer_off.sigI**2)
  subtracted.to_csv("datfiles/{}_subtracted_off.dat".format(variant), sep=" ", header=False, index=False)
  for time in TIMES:
    protein_on = pandas.read_csv("datfiles/{}_protein_{}_on.dat".format(variant, time), sep=" ", header=None, names=["q", "I", "sigI"])
    buffer_on = pandas.read_csv("datfiles/{}_buffer_{}_on.dat".format(variant, time), sep=" ", header=None, names=["q", "I", "sigI"])
    subtracted = protein_on - buffer_on
    subtracted.q = protein_on.q
    subtracted.sigI = np.sqrt(protein_on.sigI**2 + buffer_on.sigI**2)
    subtracted.to_csv("datfiles/{}_subtracted_{}_on.dat".format(variant, time), sep=" ", header=False, index=False)