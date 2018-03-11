import glob
import matplotlib.pyplot as pyplot
import pandas as pd

diffs = glob.glob('./*diff*')

low = [item for item in diffs if "3C" in item]
med = [item for item in diffs if "11C" in item]
high = [item for item in diffs if "19C" in item]

subselection = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us"]


ii = -1
for item in low:
    for sub in subselection:
        if sub in item:
            if ii < 0:
                ii=0
                nii = ii
            else:
                N = plt.cm.inferno.N
                ii += int(N/len(subselection))
                nii = N-ii    
            data = pd.read_table(item,skiprows=1,names=['q','I','sigI'],delim_whitespace=True,engine='python')
            plt.plot(data.q,data.I,label=item.split('_')[-1].replace('.dat',''),color=plt.cm.inferno(nii))
plt.legend()
plt.xscale('log')
plt.show()