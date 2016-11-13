import pandas
import numpy as np

# TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms"] # 
VARIANTS = ["WT", "S99T"]
CONCENTRATIONS = ["PC0", "PC1", "PC2"]
TEMPERATURES = [14, 21, 28]
# VARIANTS = ["WT"]
for temp in TEMPERATURES:
  for variant in VARIANTS:
    buffer = pandas.read_csv("static_dats/{}_B_{}.dat".format(variant, temp), sep=" ", header=None, names=["q", "I", "sigI"])
    for concentration in CONCENTRATIONS:
      protein = pandas.read_csv("static_dats/{}_{}_{}.dat".format(variant, concentration, temp), sep=" ", header=None, names=["q", "I", "sigI"])
      subtracted = protein - buffer
      subtracted.q = protein.q
      subtracted.sigI = np.sqrt(protein.sigI**2 + buffer.sigI**2)
      subtracted.to_csv("static_dats/{}_{}_{}_subtracted.dat".format(variant, concentration, temp), sep=" ", header=False, index=False)
