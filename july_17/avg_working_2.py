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
    print("\nException: Please use Python 3, it has been out for {} years {} days, and 2.7 may soon be unsupported (http://www.python3statement.org/)".format(years,days))
    sys.exit(1)

### import needed modules & functions
from time import clock
import parse
import pathlib
from scipy.stats import chisquare
import scipy.integrate as integrate
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
plt.style.use('ggplot')
from collections import OrderedDict
import numpy as np
from numpy.linalg import svd
from trace import Trace
import pandas as pd


# script, directory = sys.argv
# reference = parse.parse("/Volumes/beryllium/saxs_waxs_tjump/Trypsin/Trypsin-BA-Buffer-1/xray_images/Trypsin-BA-Buffer-1_26_-10us-10.tpkl")
reference = parse.parse("/Volumes/beryllium/saxs_waxs_tjump/cypa/APS_20170302/CypA-WT-1/xray_images/CypA-WT-1_9_4.22us.tpkl")
CHI_OUTLIER = 1.5
t_shortlist = ["-10.1us", "1us", "10us", "100us", "1ms"]
# t_shortlist = ["562ns"]
# t_shortlist = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us"]


def sample_map(samp_dir):
    samp_dir = pathlib.Path(samp_dir)
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

def alg_scale(self, ref):
    """Scale by projection of a vector onto a reference vector and determining the magnitude difference."""
    SA_ref = ref.SA
    SA_var = self.SA
    q = ref.q
    q_SA_ref = SA_ref*q
    q_SA_var = SA_var*q
    top = np.dot(q_SA_var,q_SA_ref)
    bottom = np.dot(q_SA_var,q_SA_var)
    scalar = top/bottom
    self.scale_factor = scalar
    self.scaled_SA = SA_var * scalar
    self.scaled_sigSA = self.sigSA * scalar
    return

def chi_stat(var, varsig, ref):
    n = len(var)-1.0
    I = var
    sigma = varsig
    Iref = ref[0]
    sigmaref = ref[1]
    chi_squared = np.sum((I-Iref)**2/(sigmaref**2))/n
    return chi_squared


def chi_outliers(vectors, reference_vector):
    chi_scores = np.apply_along_axis(chi_stat,1,vectors[:,0],vectors[:,1],reference_vector)
    print("({}, {})".format(np.mean(chi_scores), np.std(chi_scores)))
    inliers = chi_scores <= CHI_OUTLIER
    # outliers = [i for i, chi in enumerate(chi_scores) if chi>CHI_OUTLIER]
    # print("outliers= {}".format(outliers))
    return inliers


def combine_vectors_outliers(vectors):
    avg_mean = np.mean(vectors[:,0],axis=0)
    std_mean = np.std(vectors[:,0],axis=0)
    std_prop = np.sqrt(np.sum(vectors[:,1]**2 ,axis=0))/(len(vectors[:,1])-1)
    std_tot = np.sqrt(std_mean**2+std_prop**2)
    averaged_vector=(avg_mean, std_tot)
    inliers = chi_outliers(vectors, averaged_vector)
    if False not in inliers:
        pass
    else:
        clean_vectors = vectors[inliers]
        print(len(clean_vectors))
        averaged_vector = combine_vectors_outliers(clean_vectors)
    return averaged_vector


def subtract_scaled_traces(trace_one,trace_two):
    err_one = (trace_one.scale_factor*trace_one.sigSA)**2
    err_two = (trace_two.scale_factor*trace_two.sigSA)**2
    err_cov = (2*trace_one.scale_factor*trace_two.scale_factor*np.cov(trace_one.sigSA,trace_two.sigSA)[0][1])
    total_err = np.sqrt(np.abs(err_one+err_two-err_cov))
    output_SA = (trace_one.scaled_SA - trace_two.scaled_SA)
    # output = (trace_one.q, np.empty_like(output_SA, total_err)
    output = Trace(trace_one.q, np.empty_like(trace_one.q), np.empty_like(trace_one.q), total_err, output_SA, np.empty_like(trace_one.q))

    return output


def time_resolved_traces(parent, samp, reps, on_off_map):
    
    subtracted_vectors = {i: [] for i in on_off_map.keys()}
    # vectors = []

    for n in reps:
        for on, off in on_off_map.items():

            on_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, on))
            on_data = parse.parse(on_string)
            alg_scale(on_data,reference)

            off_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, off))
            off_data = parse.parse(off_string)
            alg_scale(off_data,reference)

            sub_scaled = subtract_scaled_traces(on_data,off_data)
            subtracted_vectors[on].append(sub_scaled)

            # vectors.append((off_data.SA, off_string, on_string))
            # vectors.append((on_data.SA, off_string, on_string))

    return pd.DataFrame(subtracted_vectors)


def unpack_traces(sampling):
    return[(thing.SA,thing.sigSA) for thing in sampling]

def unpack(packed_trace):
    return(packed_trace.SA,packed_trace.sigSA)

t0 = clock()
# parent, samp, reps, on_off_map = sample_map("/Volumes/beryllium/saxs_waxs_tjump/cypa/APS_20170302/CypA-WT-1/xray_images/")
# parent, samp, reps, on_off_map = sample_map(directory)
data_dir = "/Volumes/beryllium/saxs_waxs_tjump/cypa/APS_20170302/CypA-WT-1/xray_images/"

parent, samp, reps, on_off_map = sample_map(data_dir)

subtracted_vectors = time_resolved_traces(parent, samp, reps, on_off_map)


print(subtracted_vectors["-10.1us"][10].q)
# fig1, ax1 = plt.subplots(figsize=(6, 4),dpi=300)
# fig2, ax2 = plt.subplots(figsize=(6, 4),dpi=300)
# ii = -1
# for key in subtracted_vectors.keys():
#     if ii < 0:
#         ii=0
#         nii = ii
#     else:
#         N = plt.cm.inferno.N
#         ii += int(N/len(subtracted_vectors.keys()))
#         nii = N-ii
#     print("====")
#     print(key)

#     scrunch = np.stack([(item.SA,item.sigSA) for item in subtracted_vectors[key]],axis=0)
#     # for ind,row in enumerate(scrunch):
#     #     plt.scatter(ind,integrate.simps(row,sub_scaled.q),color="green")
#     biggie = [list(zip(item.q,item.SA)) for item in subtracted_vectors[key]]
#     sub_scaled = subtracted_vectors[key][0]
#     clean = combine_vectors_outliers(scrunch)
#     ### faster plotting
#     # ax = plt.axes()

#     line_segments = LineCollection(biggie,color="green",alpha=0.3,label="raw_prefilter")
#     ax2.add_collection(line_segments)
#     ax2.fill_between(sub_scaled.q, clean[0]+3*clean[1], clean[0]-3*clean[1], facecolor='red', alpha=0.3,label=r"3 $\sigma$")
#     ax2.plot(sub_scaled.q,clean[0],color="blue",alpha=1.0,label="mean_postfilter")
#     ax2.plot(sub_scaled.q,np.median(scrunch[:,0], axis=0),color="orange",label="median_prefilter")
#     ax2.set_ylabel(r'$\Delta I$ (A.U.)')
#     ax2.set_xlabel(r'q ($\AA^{-1}$)')
#     ax2.set_xscale("log")
#     ax2.set_title("{}_{}".format(samp,key))
#     ax2.legend()
#     fig2.tight_layout()
#     fig2.savefig("{}_{}_cleaned_scaled_TR.png".format(samp,key))
#     ax2.cla()

#     ### faster plotting
#     # ax.set_xscale("log")
#     # plt.show()
#     ax1.plot(sub_scaled.q, clean[0], label="{}".format(key), color=plt.cm.inferno(nii))
#     ax1.fill_between(sub_scaled.q, clean[0]+3*clean[1], clean[0]-3*clean[1], facecolor=plt.cm.inferno(nii), alpha=0.3)


# handles, labels = fig1.gca().get_legend_handles_labels()
# by_label = OrderedDict(zip(labels, handles))
# ax1.set_title('{} - Time Resolved Signal'.format(samp))
# ax1.set_ylabel(r'$\Delta I$ (A.U.)')
# ax1.set_xlabel(r'q ($\AA^{-1}$)')
# ax1.set_xscale("log")
# fig1.legend(by_label.values(), by_label.keys())
# fig1.tight_layout()
# fig1.savefig("{}_scaled_TR".format(samp))



# t1 = clock()
# print("TR signal scaled & plotted for {:20s} took {:4.2f} seconds".format(samp, t1-t0))


# matrix = np.matrix([i[0][6:] for i in vectors]).transpose()
# # matrix = np.matrix(subbed_vectors).transpose()

# u,s,v = svd(matrix, full_matrices=False)


# fig3, ax3 = plt.subplots()
# i = 0
# for vector in v.tolist()[0:5]:
#     ax3.plot([value+i*0.3 for value in vector], "-")
#     i+=1
# fig3.savefig("timepoints.png")

# fig4, ax4 = plt.subplots()
# j = 0

# for vector in u.transpose().tolist()[0:5]:
#     x = on_data.q[6:]    
#     ax4.plot(x, [value for value in vector], "-")
#     j+=0.3
# ax4.set_xscale("log")
# fig4.savefig("Singular_vectors_WT_HD_Unscaled.png")

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


    
# fig5, ax5= plt.subplots()
# ax5.plot([np.log(i) for i in s][0:10], "-")
# fig5.savefig("Singular_values_WT_HD_Unscaled.png")

# fig6, ax6 = plt.subplots()
# ax6.axhline(average_on+2.5*std_on)
# ax6.axhline(average_on-2.5*std_on)
# ax6.axhline(average_on, color="0.5")
# ax6.plot(zero_vector_values_on)
# fig6.savefig("On_vectors.png")
# ax6.cla()
# ax6.axhline(average_off+2.5*std_off)
# ax6.axhline(average_off-2.5*std_off)
# ax6.axhline(average_off, color="0.5")
# ax6.plot(zero_vector_values_off)
# fig6.savefig("Off_vectors.png")
# ax6.cla()
# ax6.axhline(average_off+2.5*std_off)
# ax6.axhline(average_off-2.5*std_off)
# ax6.axhline(average_off, color="0.5")
# ax6.plot(v.tolist()[0])
# fig6.savefig("All_vectors_off_lines.png")


    # scrunch_mean = np.mean(scrunch[:,0], axis=0)
    # scrunch_mean_arr = np.array([scrunch_mean for i in range(len(scrunch[:,0]))])
    # chistat, p_val = chisquare(scrunch[:,0],f_exp=scrunch_mean_arr,axis=1)
    # chistat = chistat/49
    # badx = np.where(chistat>1.5)
    # print("chibad = {}".format(set(badx[0])))

    # zscrunch = stats.zscore(scrunch)
    # zcum = np.sum(zscrunch,axis=1)/len(zscrunch)
    # baddx = np.where(zcum>(3*49))
    # badddx = np.where(zscrunch>5)
    # print("zbad = {}".format(set(badddx[0])))
    # print("zbad_cum = {}".format(baddx[0]))
    # scrunch_clean = scrunch[chi_stat<1.5]


