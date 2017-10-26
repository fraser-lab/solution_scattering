import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
plt.style.use("ggplot")

def kratky_plotter(sample):
    plt.figure(num=1, figsize=(4,6), dpi=300)
    x = sample.q
    y = sample.SA*sample.q**2
    data_mask = np.array(x, dtype=bool)
    data_mask[x>0.3]=False
    x = x[data_mask]
    y = y[data_mask]
    plt.plot(x,y)
    plt.xlabel("$q$")
    plt.ylabel("$I*q^2$")
    plt.title("Kratky Analysis")
    plt.show()
    return

def guinier_plotter(sample):
    plt.subplots(figsize=(4,6), dpi=300)
    plt.subplot(211)
    x = sample.q**2
    y = np.log(sample.SA)
    data_mask = np.array(x, dtype=bool)
    data_mask[x>0.008]=False
    x = x[data_mask]
    y = y[data_mask]
    fit = np.polyfit(x,y,1)
    fit_fxn = np.poly1d(fit)
    plt.scatter(x,y)
    plt.plot(x,fit_fxn(x))
    plt.xlabel("$q^2$")
    plt.ylabel("$\ln(I)$")
    plt.xlim(0.0,0.008)
    plt.title("Guinier Analysis")
    plt.subplot(212)
    plt.scatter(x,y-fit_fxn(x))
    plt.xlabel("$q^2$")
    plt.ylabel("$\ln(I)$")
    plt.xlim(0.0,0.008)
    plt.title("Guinier Residuals")
    plt.tight_layout()
    plt.show()
    return

def real_space_plotter(samples):

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

    for sample in samples:

        if ii < 0:
            ii=0
            nii = ii
        else:
            N = plt.cm.inferno.N
            ii += int(N/len(samples))
            nii = N-ii

        x1 = sample.q
        y1 = sample.SA
        mask1 = np.array(x1, dtype=bool)
        mask1[x1>0.5]=False
        mask1[x1<0.03]=False
        x1 = x1[mask1]
        y1 = y1[mask1]
        ax1.plot(x1,y1, color=plt.cm.inferno(nii))
        ax1.set_yscale('log')
        ax1.set_xlabel("$q$")
        ax1.set_ylabel("$\ln(I)$")
        ax1.set_title("Raw Scattering")

        x2 = sample.q**2
        y2 = np.log(sample.SA)
        mask2 = np.array(x2, dtype=bool)
        mask2[x2>0.008]=False
        mask2[x2<0.00125]=False
        x2 = x2[mask2]
        y2 = y2[mask2]
        fit = np.polyfit(x2,y2,1)
        fit_fxn = np.poly1d(fit)
        ax2.scatter(x2,y2, color=plt.cm.inferno(nii))
        ax2.plot(x2,fit_fxn(x2), color=plt.cm.inferno(nii))
        ax2.set_xlabel("$q^2$")
        ax2.set_ylabel("$\ln(I)$")
        ax2.set_xlim(0.0,0.008)
        ax2.set_title("Guinier Analysis")

        ax4.scatter(x2,y2-fit_fxn(x2), color=plt.cm.inferno(nii))
        ax4.set_xlabel("$q^2$")
        ax4.set_ylabel("$\ln(I)$")
        ax4.set_xlim(0.0,0.008)
        ax4.set_title("Guinier Residuals")

        x3 = sample.q
        y3 = sample.SA*sample.q**2
        mask3 = np.array(x3, dtype=bool)
        mask3[x3>0.3]=False
        x3 = x3[mask3]
        y3 = y3[mask3]
        ax3.plot(x3,y3, color=plt.cm.inferno(nii))
        ax3.set_xlabel("$q$")
        ax3.set_ylabel("$I*q^2$")
        ax3.set_title("Kratky Analysis")
    
    plt.legend()
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()
    return

from scipy.stats import linregress

def second_virial_calc(sample, sample_conc):
    plt.figure(num=1, figsize=(4,6), dpi=300)
    n=0
    sv_x = []
    sv_y = []
    I_0z = []
    for item in sample:
        x = item.q**2
        y = np.log(item.SA)
        data_mask = np.array(x, dtype=bool)
        data_mask[x>0.008]=False
        x = x[data_mask]
        y = y[data_mask]
        I_0 = np.exp(linregress(x,y)[1])
        I_0z.append(I_0)
        sv_y.append(1/I_0)
        sv_x.append(sample_conc[n])
        n+=1
    sv_xp = np.array(sv_x)
    sv_yp = np.array(sv_y)
    fit_I0 = np.polyfit(sv_xp,sv_yp,1)
    fit_fxnI0 = np.poly1d(fit_I0)
    plt.scatter(sv_xp,sv_yp)
    plt.plot(sv_xp,fit_fxnI0(sv_xp))
    plt.xlabel("$concentration (g/mL)$")
    plt.ylabel("$1/I(c,0)$")
    plt.title("Second Virial Plot")
    plt.legend()
    plt.show()
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
    return


