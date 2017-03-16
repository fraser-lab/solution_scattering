from imp import load_source
from numpy import load, column_stack, savetxt
from glob import glob
import sys

script, main_dir = sys.argv

### tell script where to find table.py
table_script = '/Users/student/Desktop/solution_scattering/table.py'

def tpkl2dat(oldfile, newfile):
    table = load_source('table', table_script)
    from table import table
    convertme = load(oldfile)
    ally = column_stack((convertme.q, convertme.SA, convertme.sigSA))
    savetxt(newfile, ally, fmt='%f', delimiter='    ', newline='\n', header='q    I    sigI')


list_of_dirs = glob(main_dir+"*/")
holo = dict.fromkeys(list_of_dirs)
for item in list_of_dirs:
    addme = glob(item+"xray_images/*.tpkl")
    holo[item] = addme


for file_dir in list_of_dirs:
    i = 1
    print("Handling data in {}".format(file_dir))
    for item in holo[file_dir]:
        print("Converting file #{} of {}.".format(i, len(holo[file_dir])))
        unfiled = item[:-5]
        filed = unfiled+".dat"
        tpkl2dat(item,filed)
        i+=1
