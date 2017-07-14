### use python 3 error & statement
import datetime
import math

now = datetime.datetime.now()
release = datetime.datetime(2008, 12, 3)
deltaT = now-release
years  = math.floor((deltaT.days) / 365)
days = deltaT.days - years*365

import sys
if sys.version_info[0] < 3:
    raise Exception("Please use Python 3, it has been out for {} years {} days, and 2.7 may soon be unsupported (http://www.python3statement.org/)".format(years,days))

### import needed modules & functions
from time import clock
from parse import parse
import pathlib


script, directory = sys.argv
reference = parse("/Volumes/beryllium/saxs_waxs_tjump/Trypsin/Trypsin-BA-Buffer-1/xray_images/Trypsin-BA-Buffer-1_26_-10us-10.tpkl")





def sample_map(direct):
    data_dir = pathlib.Path(direct)
    data_directories = [item for item in data_dir.iterdir() if item.is_dir()]
    buffer_datasets = []
    protein_datasets = []
    for item in data_directories:
        if 'static' in item.name:
            pass
        elif "Buffer" in item.name:
            buffer_datasets.append(item)
        else:
            protein_datasets.append(item)

    buffer_d = [item.name for item in buffer_datasets]
    buffer_dd = [item.replace("-Buffer-","-") for item in buffer_d]
    protein_d = [item.name for item in protein_datasets]

    ind_dict = dict((k,i) for i,k in enumerate(buffer_dd))
    matched_samples = set(protein_d).intersection(set(buffer_dd))
    indices = [ ind_dict[x] for x in matched_samples ]
    indices = sorted(indices)

    good_protein = [item for item in data_directories if item.name in matched_samples]
    good_buffers = [buffer_datasets[i] for i in indices]

    good_map = dict(zip(good_protein,good_buffers))

    return good_map


# print(sample_map(directory))
for protein, buff in sample_map(directory).items():
    # print(protein)
    # print(buff)

    protein_files = list(protein.glob(pattern='**/*.tpkl'))
    buffer_files = list(buff.glob(pattern='**/*.tpkl'))
    t0 = clock()
    for file in protein_files:
        on = parse(str(file))
        on.alg_scale(reference)
    t1 = clock()
    print("parsing & scaling for {} took {} seconds".format(protein.name, t1-t0))
    for file in buffer_files:
        on = parse(str(file))
        on.alg_scale(reference)
    t2 = clock()
    print("parsing & scaling for {} took {} seconds".format(buff.name, t2-t1))
    




