

### for use with pkl or dat files
### input can be single profile or list of profiles
### if input is a list, user must provide second list with names of samples
def kratky_plotter(sample, sample_names=None):
    if type(sample) is list:
        plt.figure(num=1, figsize=(15,10))
        n=0
        for item in sample:
            x = item.q
            y = item.SA*item.q**2
            data_mask = np.array(x, dtype=bool)
            data_mask[x>0.3]=False
            x = x[data_mask]
            y = y[data_mask]
            plt.plot(x,y, label=sample_names[n])
            n+=1
        plt.xlabel("$q$")
        plt.ylabel("$I*q^2$")
        plt.title("Kratky Analysis")
        plt.legend()
        plt.show()
    elif type(sample) is table:
        plt.figure(num=1, figsize=(15,10))
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
    elif type(sample) is pd.core.frame.DataFrame:
        plt.figure(num=1, figsize=(15,10))
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
    else:
        print("profile format not supported")
    return



### for use with pkl or dat files
### input can be single profile or list of profiles
### if input is a list, user must provide second list with names of samples
def guinier_plotter(sample, sample_names=None):
    if type(sample) is list:
        plt.figure(num=1, figsize=(15,10))
        n=0
        for item in sample:
            x = item.q**2
            y = np.log(item.SA)
            data_mask = np.array(x, dtype=bool)
            data_mask[x>0.008]=False
            x = x[data_mask]
            y = y[data_mask]
            fit = np.polyfit(x,y,1)
            fit_fxn = np.poly1d(fit)
            plt.scatter(x,y, label=sample_names[n])
            plt.plot(x,fit_fxn(x), label=sample_names[n]+"_fitted")
            n+=1
        plt.xlabel("$q^2$")
        plt.ylabel("$\ln(I)$")
        plt.xlim(0.0,0.008)
        plt.title("Guinier Analysis")
        plt.legend()
        plt.show()
    elif type(sample) is table:
        plt.figure(num=1, figsize=(15,10))
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
        plt.show()
    elif type(sample) is pd.core.frame.DataFrame:
        plt.figure(num=1, figsize=(15,10))
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
        plt.show()
    else:
        print("profile format not supported")
    return


def guinier_residual_plotter(sample, sample_names=None):
    if type(sample) is list:
        plt.figure(num=1, figsize=(15,10))
        n=0
        for item in sample:
            x = item.q**2
            y = np.log(item.SA)
            data_mask = np.array(x, dtype=bool)
            data_mask[x>0.008]=False
            x = x[data_mask]
            y = y[data_mask]
            fit = np.polyfit(x,y,1)
            fit_fxn = np.poly1d(fit)
            plt.scatter(x,y-fit_fxn(x), label=sample_names[n])
#             plt.plot(x,fit_fxn(x), label=sample_names[n]+"_fitted")
            n+=1
        plt.xlabel("$q^2$")
        plt.ylabel("$\ln(I)$")
        plt.xlim(0.0,0.008)
        plt.title("Guinier Residuals")
        plt.legend()
        plt.show()
    elif type(sample) is table:
        plt.figure(num=1, figsize=(15,10))
        x = sample.q**2
        y = np.log(sample.SA)
        data_mask = np.array(x, dtype=bool)
        data_mask[x>0.008]=False
        x = x[data_mask]
        y = y[data_mask]
        fit = np.polyfit(x,y,1)
        fit_fxn = np.poly1d(fit)
        plt.scatter(x,y-fit_fxn(x))
#         plt.plot(x,fit_fxn(x))
        plt.xlabel("$q^2$")
        plt.ylabel("$\ln(I)$")
        plt.xlim(0.0,0.008)
        plt.title("Guinier Residuals")
        plt.show()
    elif type(sample) is pd.core.frame.DataFrame:
        plt.figure(num=1, figsize=(15,10))
        x = sample.q**2
        y = np.log(sample.SA)
        data_mask = np.array(x, dtype=bool)
        data_mask[x>0.008]=False
        x = x[data_mask]
        y = y[data_mask]
        fit = np.polyfit(x,y,1)
        fit_fxn = np.poly1d(fit)
        plt.scatter(x,y-fit_fxn(x))
#         plt.plot(x,fit_fxn(x))
        plt.xlabel("$q^2$")
        plt.ylabel("$\ln(I)$")
        plt.xlim(0.0,0.008)
        plt.title("Guinier Residuals")
        plt.show()
    else:
        print("profile format not supported")
    return



def second_virial(sample, sample_conc=None):
    from scipy.stats import linregress
    if type(sample) is list:
        plt.figure(num=1, figsize=(15,10))
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


