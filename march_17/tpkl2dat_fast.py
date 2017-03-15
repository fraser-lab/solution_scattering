from imp import load_source
from numpy import load
from glob import glob
import sys

script, main_dir = sys.argv

### tell script where to find table.py
table_script = '/Users/student/Desktop/solution_scattering/table.py'

def tpkl2dat(oldfile, newfile):
    table = imp.load_source('table', table_script)
    from table import table
    convertme = np.load(oldfile)
    ally = np.column_stack((convertme.q, convertme.SA, convertme.sigSA))
    np.savetxt(newfile, ally, fmt='%f', delimiter='    ', newline='\n', header='q    I    sigI')


list_of_dirs = glob(main_dir+"*/")
holo = dict.fromkeys(list_of_dirs)
for item in list_of_dirs:
    addme = glob(item+"xray_images/*.tpkl")
    holo[item] = addme

i = 1
for item in holo[list_of_dirs[8]]:
    print("Converting file #{} of {}.".format(i, len(holo[list_of_dirs[8]])))
    unfiled = item[:-5]
    filed = unfiled+".dat"
    tpkl2dat(item,filed)
    sys.stdout.flush()
    i+=1
