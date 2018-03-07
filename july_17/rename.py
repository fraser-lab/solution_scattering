import subprocess
import glob
from sys import argv
import os


script, dir = argv


files = glob.glob(dir+"/**.tpkl")



for item in files:
    # name, meg, n, time = item.split("_")
    filename = os.path.basename(item)
    dirname = os.path.dirname(item)
    # samp, meg, rep, time, onoff = filename.split("_")
    samp, rep, time, onoff = filename.split("_")

    # combo = int(meg)*int(rep)
    oldname = dirname+'/'+filename
    # newname = dirname+'/'+samp+'_'+str(combo)+'_'+time+'_'+onoff
    newname = dirname+'/'+samp+'_'+rep+'_'+time+'.tpkl'



    subprocess.call(['mv', oldname, newname])

    # print(name)