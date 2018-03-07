import pandas as pd
import matplotlib.pyplot as plt
import parse
import pathlib
from new_analysis import *


data_dir = "/Volumes/beryllium/solution_scattering/july_17"
samp_dir = pathlib.Path(data_dir)
samp_files = list(samp_dir.glob(pattern='./CypA-S99T*.dat'))
# samp_files = list(samp_dir.glob(pattern='./CypA-WT*.dat'))

def name_to_time(name):
    names = name.split('_')

    time = str(names[-1])
    time = time.replace(".dat","")
    # t = int(list(filter(str.isfloat, time))[0])

    if "us" in time:
        
        time=time.replace("us","")
        t = float(time)
        t=t*1000
    elif "ms" in time:
        
        time=time.replace("ms","")
        t = float(time)
        t=t*1000000

    else:
        time=time.replace("ns","")
        t = float(time)
    
    return t


in_files = [(str(item), name_to_time(str(item)),parse.parse(str(item)))for item in samp_files]

# in_files.sort(key=lambda x: x[1])

# on_data = pd.read_table("CypA-WT-1_100us.dat", delim_whitespace=True, engine='python', skiprows=1, names=['q','SA','sigSA'])
# off_data = pd.read_table("CypA-WT-1_-10.1us.dat", delim_whitespace=True, engine='python', skiprows=1, names=['q','SA','sigSA'])

# for file in samp_files:
# 	data = parse.parse(str(file))
# 	plt.plot(data.q,data.SA)
# plt.plot(on_data.q[2:],on_data.SA[2:], label="100us", lw=3)
# plt.plot(off_data.q[2:],off_data.SA[2:], label="off", lw=1.5)
# plt.legend()
# plt.xscale('log')
# # plt.savefig("cypa_wt1_corrected_curves.png")
# plt.show()

# real_space_plotter(in_files)


def rg_vs_t(samples):
	
    if isinstance(samples,list):
        pass
    else:
        samples = [samples]
    
    fig=plt.figure(figsize=(6,6),dpi=100)
    fig.suptitle("Real Space Analysis")
    gs=GridSpec(2,2) # 2 rows, 2 columns
    ax1=fig.add_subplot(gs[0,0]) # First row, first column
    ax2=fig.add_subplot(gs[0,1]) # First row, second column
    ax3=fig.add_subplot(gs[1,0]) # Second row, first column
    ax4=fig.add_subplot(gs[1,1]) # Second row, second column
    ii = -1

    samples.sort(key=lambda x: x[1])
    for sample in samples:

        if ii < 0:
            ii=0
            nii = ii
        else:
            N = plt.cm.inferno.N
            ii += int(N/len(samples))
            nii = N-ii

        # names = sample[0].split('_')

        # time = str(names[-1])
        # time = time.replace(".dat","")
        # # t = int(list(filter(str.isfloat, time))[0])

        # if "us" in time:
            
        #     time=time.replace("us","")
        #     t = float(time)
        #     t=t*1000
        # if "ms" in time:
            
        #     time=time.replace("ms","")
        #     t = float(time)
        #     t=t*1000000

        # t = int(filter(str.isdigit, time))
        t = sample[1]
        x2 = sample[2].q**2
        y2 = np.log(sample[2].SA)
        mask2 = np.array(x2, dtype=bool)
        mask2[x2>0.01]=False
        mask2[x2<0.0010]=False
        x2 = x2[mask2]
        y2 = y2[mask2]
        fit = np.polyfit(x2,y2,1)
        fit_fxn = np.poly1d(fit)
        m,b = np.poly1d(fit)
        rg = np.sqrt(-3*m)
        ax2.scatter(t,rg, color=plt.cm.inferno(nii))
        # ax2.plot(x2,fit_fxn(x2), color=plt.cm.inferno(nii))
        ax2.set_xlabel("t (nanosec)")
        ax2.set_ylabel("RG")
        ax2.set_xscale('log')
        # ax2.set_xlim(0.0,0.008)
        ax2.set_title("RG vs Time")

        ax1.scatter(t,b, color=plt.cm.inferno(nii))
        # ax2.plot(x2,fit_fxn(x2), color=plt.cm.inferno(nii))
        ax1.set_xlabel("t (nanosec)")
        ax1.set_ylabel("I_0")
        ax1.set_xscale('log')
        # ax2.set_xlim(0.0,0.008)
        ax1.set_title("I_0 vs Time")


        x3 = sample[2].q
        y3 = sample[2].SA
        ax3.plot(x3[x3>0.03],y3[x3>0.03], color=plt.cm.inferno(nii))
        # ax2.plot(x2,fit_fxn(x2), color=plt.cm.inferno(nii))
        ax3.set_xlabel("q")
        ax3.set_ylabel("I")
        ax3.set_xscale('log')
        # ax2.set_xlim(0.0,0.008)
        ax3.set_title("Raw Scattering")




    plt.legend()
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    # plt.savefig("S99T_TR_RS.png", dpi=300)
    plt.show()
    return

rg_vs_t(in_files)


# import parse

# test = parse.parse("./cypa_wt1_diff_march2017/CypA-WT-1_diff_10us.dat")
# plt.plot(test.q,test.SA)
# plt.xscale("log")
# plt.show()

# print(test)