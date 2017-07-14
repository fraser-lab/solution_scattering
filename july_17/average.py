from sys import argv
from time import clock
from parse import parse


directories = argv[1:]
reference = parse("/Volumes/beryllium/saxs_waxs_tjump/Trypsin/Trypsin-BA-Buffer-1/xray_images/Trypsin-BA-Buffer-1_26_-10us-10.tpkl")

t0 = clock()
import pathlib
data_dir = pathlib.Path(directories[0])
files = list(data_dir.glob(pattern='**/*.tpkl'))
for file in files:
    on = parse(str(file))
    on.alg_scale(reference)

t1 = clock()
print("parsing & scaling took {} seconds".format(t1-t0))




