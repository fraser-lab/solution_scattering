### use python 3 error & statement
import datetime
import math

now = datetime.datetime.now()
release = datetime.datetime(2008, 12, 3)
deltaT = now-release
years  = math.floor((deltaT.days) / 365)
days = deltaT.days - years*365

import sys
if sys.version_info[0] < 3:
    raise Exception("Please use Python 3, it has been out for {} years {} days, and 2.7 may soon be unsupported (http://www.python3statement.org/)".format(years,days))

### import needed modules & functions
from time import clock
import parse
import pathlib


script, directory = sys.argv
reference = parse.parse("/Volumes/beryllium/saxs_waxs_tjump/Trypsin/Trypsin-BA-Buffer-1/xray_images/Trypsin-BA-Buffer-1_26_-10us-10.tpkl")





def experiment_map(direct):
    data_dir = pathlib.Path(direct)
    data_directories = [item for item in data_dir.iterdir() if item.is_dir()]
    buffer_datasets = []
    protein_datasets = []
    for item in data_directories:
        if 'static' in item.name:
            pass
        elif "Buffer" in item.name:
            buffer_datasets.append(item)
        else:
            protein_datasets.append(item)

    buffer_d = [item.name for item in buffer_datasets]
    buffer_dd = [item.replace("-Buffer-","-") for item in buffer_d]
    protein_d = [item.name for item in protein_datasets]

    ind_dict = dict((k,i) for i,k in enumerate(buffer_dd))
    matched_samples = set(protein_d).intersection(set(buffer_dd))
    indices = [ ind_dict[x] for x in matched_samples ]
    indices = sorted(indices)

    good_protein = [item for item in data_directories if item.name in matched_samples]
    good_buffers = [buffer_datasets[i] for i in indices]

    good_map = dict(zip(good_protein,good_buffers))

    return good_map




def sample_map(samp_dir):
    samp_files = list(samp_dir.glob(pattern='**/*.tpkl'))
    # buffer_files = list(buff.glob(pattern='**/*.tpkl'))
    t0 = clock()
    REPS = []
    TIMES = []
    for file in samp_files:
        # on = parse.parse(str(file))
        # on.alg_scale(reference)
        name = file.name
        parent = file.parent
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
    on_off_map = (dict(zip(clean_ons,sorted(OFFS, key=parse.alphanum_key))))

    return parent, samp, REPS, on_off_map


t_shortlist = ["-10.1us", "1us", "10us", "100us", "1ms"]
# t_shortlist = ["562ns"]
# t_shortlist = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us"]
# vectors = []
# print(sample_map(directory))
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from collections import OrderedDict
import numpy as np
from numpy.linalg import svd
from trace import Trace

### all that is needed for mpl to plotly conversions
# from plotly.offline import init_notebook_mode, plot_mpl


for protein, buff in experiment_map(directory).items():

    t0 = clock()

    parent, samp, reps, on_off_map = sample_map(protein)
    b_parent, b_samp, b_reps, b_on_off_map = sample_map(buff)
    plt.figure(figsize=(10,5),dpi=300)

    for n in reps:
        i = -1
        for on, off in on_off_map.items():

            if on in t_shortlist:

                on_data = parse.parse("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, on))
                on_data.alg_scale(reference)
                on_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, on))

                off_data = parse.parse("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, off))
                off_data.alg_scale(reference)
                off_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, off))

                buf_on_data = parse.parse("{0}/{1}_{2}_{3}.tpkl".format(b_parent, b_samp, n, on))
                buf_on_data.alg_scale(reference)
                buf_on_string = "{0}/{1}_{2}_{3}.tpkl".format(b_parent, b_samp, n, on)

                buf_off_data = parse.parse("{0}/{1}_{2}_{3}.tpkl".format(b_parent, b_samp, n, off))
                buf_off_data.alg_scale(reference)
                buf_off_string = "{0}/{1}_{2}_{3}.tpkl".format(b_parent, b_samp, n, off)

                # sub = on_data.scaled_SA - off_data.scaled_SA
                subSA = (on_data.scaled_SA - off_data.scaled_SA)
                subSAsig = np.sqrt((on_data.scale_factor*on_data.sigSA)**2+(off_data.scale_factor*off_data.sigSA)**2-(2*on_data.scale_factor*off_data.scale_factor*np.cov(on_data.sigSA,off_data.sigSA)[0][1]))
                sub_scaled = Trace(on_data.q, None, None, subSAsig, subSA, None)

                buf_subSA = (buf_on_data.scaled_SA - buf_off_data.scaled_SA)
                buf_subSAsig = np.sqrt((buf_on_data.scale_factor*buf_on_data.sigSA)**2+(buf_off_data.scale_factor*buf_off_data.sigSA)**2-(2*buf_on_data.scale_factor*buf_off_data.scale_factor*np.cov(buf_on_data.sigSA,buf_off_data.sigSA)[0][1]))
                buf_sub_scaled = Trace(buf_on_data.q, None, None, buf_subSAsig, buf_subSA, None)
                # sub.scaled_SA = on_data.scaled_SA - off_data.scaled_SA

                # vectors.append((off_data.as_vector(), off_string, on_string))
                # vectors.append((on_data.as_vector(), off_string, on_string))
                if i < 0:
                    i=0
                    ni = i
                else:
                    # i += int(255/len(on_off_map))
                    i += int(255/len(t_shortlist))
                    ni = 255-i

                plt.subplot(2, 1, 1)
                plt.ylabel('protein')
                plt.xscale("log")
                plt.plot(sub_scaled.q, sub_scaled.SA, label="{}".format(on), color=plt.cm.inferno(ni))
                plt.tight_layout()
                plt.subplot(2, 1, 2)
                plt.ylabel('Buffer')
                plt.xscale("log")
                plt.plot(buf_sub_scaled.q, buf_sub_scaled.SA, label="{}".format(on), color=plt.cm.inferno(ni))
                plt.tight_layout()


                # plt.errorbar(sub_scaled.q, sub_scaled.SA, label="{}".format(on), color=plt.cm.inferno(ni), yerr=sub_scaled.sigSA)

            else:
                pass


    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.savefig("{}_scaled_TR".format(samp))

    ### this will work to convert mpl to plotly, quite slow though
    # fig = plt.gcf()
    # plot_mpl(fig, filename="{}_scaled_TR.html".format(samp))

    t1 = clock()
    print("TR signal scaled & plotted for {:20s} took {:4.2f} seconds".format(samp, t1-t0))

#####
def chisquared(var, ref):
    nu = len(var)-1.0
    I, sigma = zip(*var) 
    Iref,sigmaref = zip(*ref)
    chi_squared = 1/nu*sum([(I[i]-Iref[i])**2/(sigmaref[i]**2)  for i in range(len(var))])
    # print chi_squared
    return chi_squared

def chi_outliers(vectors, reference_vector):
    list = [chisquared(i, reference_vector) for i in vectors]
    print(np.mean(list), np.std(list))
    outlier_list = []
    for i, val in enumerate(list):
        if val > CHI_OUTLIER:
            outlier_list.append(i)
            # print i, val
    return outlier_list

def combine_vectors_outliers(vectors):
    averaged_vector = []
    for i in range(len(vectors[0])):
        value_list = [v[i] for v in vectors]
        # averaged_vector.append(np.mean(value_list))
        #### Ensemble Weighting - mean = sum(value/sigma^2)/sum(1/sigma^2)
        ####                                            std = sqrt(1/sum(1/sigma^2))
        means, stds = zip(*value_list)
        avg_mean = np.mean(means)
        std_mean = np.std(means)
        std_prop = np.sqrt(sum([j**2 for j in stds]))/(len(stds)-1)
        std_tot = np.sqrt(std_mean**2+std_prop**2)
        if i == 0:
            print("Standard deviation of means: ", std_mean)
            print("Propogated standard deviation from data points: ", std_prop)
            print("Total propogated standard deviation: ", std_tot)
        # avg_mean = sum([means[i]/(stds[i]**2) for i,_ in enumerate(means)])/sum([1/stds[i]**2 for i,_ in enumerate(means)])
        # std_tot = np.sqrt(1/sum([1/stds[i]**2 for i,_ in enumerate(means)]))
        # print on.q[i], avg_mean, std_mean, std_prop, std_tot
        averaged_vector.append((avg_mean, std_tot))
    outlier_list = chi_outliers(vectors, averaged_vector)
    if len(outlier_list) > 0:
        new_vectors = [vector for i,vector in enumerate(vectors) if i not in outlier_list]
        print(len(new_vectors))
        averaged_vector = combine_vectors_outliers(new_vectors)
    return averaged_vector



####
    # matrix = np.matrix([i[0][6:] for i in vectors]).transpose()
    # # matrix = np.matrix(subtracted_vectors).transpose()

    # u,s,v = svd(matrix, full_matrices=False)

    # # print u.shape
    # # print s.shape
    # # print v[0]
    # # print v[1]
    # # for vector in v[0:7]:
    # #   print vector[0]

    # # print v.slice(0)

    # fig, ax = plt.subplots()
    # i = 0
    # for vector in v.tolist()[0:5]:
    #     # print vector
    #     # ax.plot(range(len(vectors)), [value+i for value in vector], "-")
    #     ax.plot([value+i*0.3 for value in vector], "-")
    #     i+=1
    # fig.savefig("timepoints.png")

    # fig2, ax2 = plt.subplots()
    # j = 0

    # for vector in u.transpose().tolist()[0:5]:
    #     # print vector
    #     # ax.plot(range(len(vectors)), [value+i for value in vector], "-")
    #     x = on_data.q[6:]    
    #     ax2.plot(x, [value for value in vector], "-")
    #     j+=0.3
    # fig2.savefig("Singular_vectors_WT_HD_Unscaled.png")

    # zero_vector_values_on = v.tolist()[0][1::2]
    # zero_vector_values_off = v.tolist()[0][0::2]
    # average_on = np.average(zero_vector_values_on)
    # std_on = np.std(zero_vector_values_on)
    # average_off =  np.average(zero_vector_values_off)
    # std_off = np.std(zero_vector_values_off)
    # print(average_on, average_off, std_on, std_off)
    # for index, value in enumerate(v.tolist()[0]):
    #     if index % 2 == 0:
    #         if np.fabs((value - average_off) / std_off) > 2.5:
    #             print(vectors[index][1])
    #             print(vectors[index][2])
    #     elif index % 2 == 1:
    #         if np.fabs((value - average_on) / std_on) > 2.5:
    #             print(vectors[index][1])
    #             print(vectors[index][2])


        
    # fig3, ax3= plt.subplots()
    # ax3.plot([np.log(i) for i in s][0:10], "-")
    # fig3.savefig("Singular_values_WT_HD_Unscaled.png")

    # fig4, ax4 = plt.subplots()
    # ax4.axhline(average_on+2.5*std_on)
    # ax4.axhline(average_on-2.5*std_on)
    # ax4.axhline(average_on, color="0.5")
    # ax4.plot(zero_vector_values_on)
    # fig4.savefig("On_vectors.png")
    # ax4.cla()
    # ax4.axhline(average_off+2.5*std_off)
    # ax4.axhline(average_off-2.5*std_off)
    # ax4.axhline(average_off, color="0.5")
    # ax4.plot(zero_vector_values_off)
    # fig4.savefig("Off_vectors.png")
    # ax4.cla()
    # ax4.axhline(average_off+2.5*std_off)
    # ax4.axhline(average_off-2.5*std_off)
    # ax4.axhline(average_off, color="0.5")
    # ax4.plot(v.tolist()[0])
    # fig4.savefig("All_vectors_off_lines.png")

    # plt.show()

########

    # print(protein)
    # print(buff)
    # print("\n")
    # print("*"*80)
    # print("\n")
    # print("REPS = {}\n".format(REPS))
    # print("OFFS = {}\n".format(sorted(OFFS, key=parse.alphanum_key)))
    # print("ONS = {}\n".format(clean_ons))
    # print("\nONS & OFFS Match = {}".format(len(OFFS)==len(clean_ons)))
    # print(dict(zip(clean_ons,sorted(OFFS, key=parse.alphanum_key))))


    # t1 = clock()
    # print("parsing & scaling for {:30s} took {:4.2f} seconds".format(protein.name, t1-t0))
    # for file in buffer_files:
    #     on = parse.parse(str(file))
    #     on.alg_scale(reference)
    # t2 = clock()
    # print("parsing & scaling for {:30s} took {:4.2f} seconds".format(buff.name, t2-t1))
    




