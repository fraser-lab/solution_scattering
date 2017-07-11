import plotly
from plotly.graph_objs import Scatter, Layout
import pandas as pd
import numpy as np
import pickle as pkl
from parse import parse_tpkl_2
import pathlib



working_dir = pathlib.Path().cwd()
data_directories = [item for item in working_dir.iterdir() if item.is_dir()]
d1 = data_directories[0]
d1.name
d1_files = list(d1.glob(pattern='**/*.tpkl'))
len(d1_files)
first_file = d1_files[0]
name = first_file.name
samp, rep, time = name.split('_')
time = time.replace('.tpkl','')

print("Working directory: {}\n".format(working_dir))
print("Data directories:")
for item in data_directories:
    print(item)
print("\n")
print("Directory\t\t# of Files\tSample\t\tRepetition\tTimepoint")
print("-"*80)
print("{}\t\t{}\t\t{}\t{}\t\t{}".format(d1.name,len(d1_files),samp,rep,time))