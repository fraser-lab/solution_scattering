import sys
if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")
import pathlib
from sys import argv
import parse

script, data_d = argv


working_dir = pathlib.Path().cwd()
data_dir = pathlib.Path(data_d)
data_directories = [item for item in data_dir.iterdir() if item.is_dir()]


print("\n\nWorking directory: {}\n".format(working_dir))
print("Data directories:")
for item in data_directories:
    print(item)
print("\n")
print("Directory\t\t\t\t# of Files")
print("*"*80)

for item in data_directories:
    # if "static" in item.name:
    #     pass
    # else:
    files = list(item.glob(pattern='**/*.tpkl'))
    # len(files)
    # first_file = files[0]
    # name = first_file.name
    # samp, rep, time = name.split('_')
    # time = time.replace('.tpkl','')
    print("{:40}{:>}".format(item.name,len(files)))


# print("Buffer Datasets:")
buffer_datasets = []
for item in data_directories:
    if "Buffer" in item.name:
        # print(item)
        buffer_datasets.append(item)

# print("\n")
# print("Protein Datasets:")
protein_datasets = []
for item in data_directories:
    if "Buffer" not in item.name:
        if "static" not in item.name:
            # print(item)
            protein_datasets.append(item)

buffer_d = [item.name for item in buffer_datasets]
buffer_dd = [item.replace("-Buffer-","-") for item in buffer_d]
protein_d = [item.name for item in protein_datasets]

matched_samples = set(protein_d).intersection(set(buffer_dd))

# print(matched_samples)
print("\n\nWriting the following parameters to 'params.py':\n")
print("*"*80)

matched_datasets = []
for item in data_directories:
    if item.name in matched_samples:
        matched_datasets.append(item)

REPS = []
TIMES = []
for item in matched_datasets:
    files = list(item.glob(pattern='**/*.tpkl'))
    for file in files:
        name = file.name
        samp, rep, time = name.split('_')
        time = time.replace('.tpkl','')
        REPS.append(rep)
        TIMES.append(time)



REPS = sorted(list(set(REPS)), key=float)
TIMES = list(set(TIMES))
OFFS = [item for item in TIMES if "-10us" in item]
ONS = [item for item in TIMES if "-10us" not in item]

tup =  [parse.unit_sort(item) for item in ONS]
tup_sort = sorted(tup, key=lambda item: (item[0],parse.natural_keys(item[1])))
clean_ons = [item[1] for item in tup_sort]

print("REPS = {}\n".format(REPS))
print("OFFS = {}\n".format(sorted(OFFS, key=parse.alphanum_key)))
print("ONS = {}\n".format(clean_ons))
print("*"*80)

with open('params.py','w') as f:
    f.write("REPS = {}\n\n".format(REPS))
    f.write("OFFS = {}\n\n".format(sorted(OFFS, key=parse.alphanum_key)))
    f.write("ONS = {}\n\n".format(clean_ons))


print("\n\nChecking internal consistency of data:")
print("\nONS & OFFS Match = {}".format(len(OFFS)==len(clean_ons)))
print(dict(zip(clean_ons,sorted(OFFS, key=parse.alphanum_key))))



times_scaled = [parse.times_numeric(time) for time in clean_ons]
print("\nmapping times")
print(dict(zip(clean_ons,times_scaled)))
# TIMES = sorted(TIMES, key=alphanum_key)


# print("TIMES = {}".format(TIMES))



