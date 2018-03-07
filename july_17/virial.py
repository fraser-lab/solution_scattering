import parse
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
plt.style.use('ggplot')
from trace import Trace
from scipy.stats import linregress
import os
import pandas as pd
import numpy as np
from glob import glob


def second_virial_calc(samples, full_conc):
    plt.figure(num=1, figsize=(10,6), dpi=100)
    n=0
    sv_x = []
    sv_y = []
    I_0z = []
    temps = set([item[2] for item in samples])
    temps = sorted(list(temps))
    concs = [full_conc/1000, full_conc/3000, full_conc/9000]
    output = []
    for temp in temps:
        for item in samples:
            if item[2] == temp:
                x = item[3].q**2
                y = np.log(item[3].SA)
                data_mask = np.array(x, dtype=bool)
                data_mask[x>0.008]=False
                x = x[data_mask]
                y = y[data_mask]
                I_0 = np.exp(linregress(x,y)[1])
                I_0z.append(I_0)
                sv_y.append(1/I_0)
                sv_x.append(concs[int(item[1])])
                n+=1
            else:
                pass
        sv_xp = np.array(sv_x)
        sv_yp = np.array(sv_y)
        fit_I0 = np.polyfit(sv_xp,sv_yp,1)
        fit_fxnI0 = np.poly1d(fit_I0)
        plt.scatter(sv_xp,sv_yp, label=temp)
        plt.plot(sv_xp,fit_fxnI0(sv_xp), label=temp+"_fit")
        vir_stats = linregress(sv_xp,sv_yp)
        I_0_0 = 1/vir_stats[1]
        slope = vir_stats[0]
        MW = 18500
        A = slope*I_0_0/(2*MW)
        print("\nStats for virial fit:\n{}\n".format(vir_stats))
        print("I(0,0) = {}".format(I_0_0))
        print("A = {}".format(A))
        print("I(c,0) for pc0 = {}".format(I_0z[0]))
        print("I(c,0) for pc1 = {}".format(I_0z[1]))
        print("I(c,0) for pc2 = {}".format(I_0z[2]))
        output.append((temp,A))
    plt.xlabel("$concentration (g/mL)$")
    plt.ylabel("$1/I(c,0)$")
    plt.title("Second Virial Plot")
    plt.legend()
    # plt.savefig("second_vir.png")
    plt.show()

    return output

def iter_vir(samples, full_conc):
    n=0
    
    temps = set([item[2] for item in samples])
    temps = sorted(list(temps))
    concs = [full_conc/1, full_conc/3, full_conc/9]
    I0q_output = []
    a2_output = []
    spf_output = []
    for temp in temps:
        sv_x = []
        sv_y = []
        I_0z = []
        for item in samples:
            if item[2] == temp:
                # x = item[3].q
                y = item[3].SA
                # data_mask = np.array(x, dtype=bool)
                # data_mask[x>0.008]=False
                # x = x[data_mask]
                # y = y[data_mask]
                # I_0 = np.exp(linregress(x,y)[1])
                # I_0z.append(I_0)
                sv_y.append(1/y)
                sv_x.append(concs[int(item[1])])
                n+=1
            else:
                pass
        sv_xp = np.array(sv_x)
        # print(sv_xp.shape)
        sv_yp = np.stack(sv_y)
        # print(sv_yp.shape)
        m,b = np.polyfit(sv_xp,sv_yp,1)
        mw = 18500
        conc = 50
        # print(fit_I0)
            # fit_fxnI0 = np.poly1d(fit_I0)
            # plt.scatter(sv_xp,sv_yp, label=temp)
            # plt.plot(sv_xp,fit_fxnI0(sv_xp), label=temp+"_fit")
            # vir_stats = linregress(sv_xp,sv_yp)
            # I_0_0 = 1/vir_stats[1]
            # slope = vir_stats[0]
            # MW = 18500
            # A = slope*I_0_0/(2*MW)
            # print("\nStats for virial fit:\n{}\n".format(vir_stats))
            # print("I(0,0) = {}".format(I_0_0))
            # print("A = {}".format(A))
            # print("I(c,0) for pc0 = {}".format(I_0z[0]))
            # print("I(c,0) for pc1 = {}".format(I_0z[1]))
            # print("I(c,0) for pc2 = {}".format(I_0z[2]))
        I0q_output.append((temp,1/b))
        a2_output.append((temp,m*(1/b)/2*mw))
        spf_output.append((temp, 1/(1+m/b*conc)))

    return spf_output

# reference = parse.parse("/Volumes/beryllium/saxs_waxs_tjump/cypa/statics/CypA-S99T-static_0_14.dat")


s99 = glob("/Volumes/beryllium/saxs_waxs_tjump/cypa/statics/**S99T**.dat")
wt = glob("/Volumes/beryllium/saxs_waxs_tjump/cypa/statics/**WT**.dat")


s99_map = []
wt_map = []
for item in s99:
    file = os.path.basename(item)
    samp,conc,temp = file.split("_")
    temp = temp.replace(".dat","")
    on_data = pd.read_table(item, delim_whitespace=True, engine='python', skiprows=1, names=['q','SA','sigSA'])
    s99_map.append((samp,conc,temp,on_data))

for item in wt:
    file = os.path.basename(item)
    samp,conc,temp = file.split("_")
    temp = temp.replace(".dat","")
    on_data = pd.read_table(item, delim_whitespace=True, engine='python', skiprows=1, names=['q','SA','sigSA'])
    wt_map.append((samp,conc,temp,on_data))


# s99_svc = second_virial_calc(s99_map, 50)
# wt_svc = second_virial_calc(wt_map, 50)
print(wt_map)
wt_full = iter_vir(wt_map, 50)

# print(wt_full)

plt.plot(wt_map[0][3].q,wt_full[0][1], label=wt_full[0][0])
plt.plot(wt_map[0][3].q,wt_full[1][1], label=wt_full[1][0])
plt.plot(wt_map[0][3].q,wt_full[2][1], label=wt_full[2][0])
# plt.xticks(wt_map[0][3].q)
# plt.xscale('log')
plt.yscale('log')
plt.xlim(0,0.15)
plt.ylim(10,40)
plt.legend()
plt.show()


def vir_temp(svcs, names):
    n=0
    for svc in svcs:
        x_vir = np.array([float(item[0]) for item in svc])
        y_vir = np.array([float(item[1]) for item in svc])
        y_vir = y_vir/10e-4
        plt.figure(num=1, figsize=(10,6), dpi=100)
        plt.scatter(x_vir, y_vir, label=names[n])
        fit_virtemp = np.polyfit(x_vir,y_vir,1)
        fit_fxn_virtemp = np.poly1d(fit_virtemp)
        plt.plot(x_vir,fit_fxn_virtemp(x_vir), label=names[n]+"_fit")
        # from scipy.stats import linregress
        temp_stats = linregress(x_vir,y_vir)
        print("{}".format(temp_stats))
        n+=1
    plt.ylabel("$A (10^{-4}*mol*ml*g^2)$")
    plt.xlabel("$T (^o C)$")
    plt.title("Second Virial Coefficients")
    plt.legend()
    plt.show()

    return


# vir_temp([s99_svc,wt_svc], ["S99", "WT"])


