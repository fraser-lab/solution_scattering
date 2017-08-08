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
from new_analysis import real_space_plotter


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
    on_off_map = {k:v for k,v in (zip(clean_ons,sorted(OFFS, key=parse.alphanum_key)))}

    return parent, samp, REPS, on_off_map


def chi_stat(var, ref):
    n = len(var.q)-1.0
    chi_squared = np.sum((var.SA-ref.SA)**2/(ref.sigSA**2))/n
    return chi_squared


def chi_outliers(vectors, reference_vector):
    chi_scores = [chi_stat(i, reference_vector) for i in vectors]
    print("({}, {})".format(np.mean(chi_scores), np.std(chi_scores)))
    inliers = np.asarray(chi_scores) <= CHI_OUTLIER
    # outliers = [i for i, chi in enumerate(chi_scores) if chi>CHI_OUTLIER]
    # print("outliers= {}".format(outliers))
    return inliers

def average_traces(traces):
    one_curve = traces[0]
    mean_SA = np.mean([trace.SA for trace in traces],axis=0)
    std_err = np.std([trace.SA for trace in traces],axis=0)
    prop_err = np.sqrt(np.sum([trace.sigSA**2 for trace in traces],axis=0))/(len(traces)-1)
    tot_err = np.sqrt(std_err**2+prop_err**2)
    averaged_vector=Trace(one_curve.q, np.empty_like(one_curve.q), np.empty_like(one_curve.q), tot_err, mean_SA, np.empty_like(one_curve.q))
    return averaged_vector

def iterative_chi_filter(vectors):
    averaged_vector=average_traces(vectors)
    inliers = chi_outliers(vectors, averaged_vector)
    if False not in inliers:
        clean_vectors = vectors
    else:
        clean_vectors = [vectors[i] for i, x in enumerate(inliers) if x]
        print(len(clean_vectors))
        iterative_chi_filter(clean_vectors)
    return clean_vectors


def subtract_scaled_traces(trace_one,trace_two):
    err_one = (trace_one.scale_factor*trace_one.sigSA)**2
    err_two = (trace_two.scale_factor*trace_two.sigSA)**2
    err_cov = (2*trace_one.scale_factor*trace_two.scale_factor*np.cov(trace_one.sigSA,trace_two.sigSA)[0][1])
    total_err = np.sqrt(np.abs(err_one+err_two-err_cov))
    output_SA = (trace_one.scaled_SA - trace_two.scaled_SA)
    output = Trace(trace_one.q, np.empty_like(trace_one.q), np.empty_like(trace_one.q), total_err, output_SA, np.empty_like(trace_one.q))
    return output


def add_unscaled_traces(trace_one,trace_two):
    err_one = trace_one.sigSA**2
    err_two = trace_two.sigSA**2
    err_cov = (2*np.cov(trace_one.sigSA,trace_two.sigSA)[0][1])
    total_err = np.sqrt(np.abs(err_one+err_two+err_cov))
    output_SA = (trace_one.SA + trace_two.SA)
    output = Trace(trace_one.q, np.empty_like(trace_one.q), np.empty_like(trace_one.q), total_err, output_SA, np.empty_like(trace_one.q))
    return output

def subtract_unscaled_traces(trace_one,trace_two):
    err_one = trace_one.sigSA**2
    err_two = trace_two.sigSA**2
    err_cov = (2*np.cov(trace_one.sigSA,trace_two.sigSA)[0][1])
    total_err = np.sqrt(np.abs(err_one+err_two-err_cov))
    output_SA = (trace_one.SA - trace_two.SA)
    output = Trace(trace_one.q, np.empty_like(trace_one.q), np.empty_like(trace_one.q), total_err, output_SA, np.empty_like(trace_one.q))
    return output


def time_resolved_traces(parent, samp, reps, on_off_map):
    
    subtracted_vectors = {i: [] for i in on_off_map.keys()}

    for n in reps:
        for on, off in on_off_map.items():

            on_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, on))
            on_data = parse.parse(on_string)
            on_data.alg_scale(reference)

            off_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, off))
            off_data = parse.parse(off_string)
            off_data.alg_scale(reference)

            sub_scaled = subtract_scaled_traces(on_data,off_data)
            subtracted_vectors[on].append(sub_scaled)

    return subtracted_vectors

def all_off_traces(parent, samp, reps, on_off_map):
    
    off_vectors = []

    for n in reps:
        for off in on_off_map.values():

            off_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, off))
            off_data = parse.parse(off_string)
            off_data.alg_scale(reference)
            off_scaled = Trace(off_data.q, np.empty_like(off_data.q), np.empty_like(off_data.q), off_data.scaled_sigSA, off_data.scaled_SA, off_data.Nj)
            off_vectors.append(off_scaled)

    return off_vectors


def plot_all_traces(samp, trace_lib):
    fig, ax = plt.subplots(figsize=(6, 4),dpi=300)
    for key in trace_lib.keys():
        print("====")
        print("plotting all traces for {}".format(key))
        all_curves = [list(zip(item.q,item.SA)) for item in trace_lib[key]]
        mean_trace = average_traces(trace_lib[key])
        median_SA = np.median([trace.SA for trace in trace_lib[key]],axis=0)

        ### faster plotting
        # ax = plt.axes()
        line_segments = LineCollection(all_curves,color="green",alpha=0.3,label="raw_prefilter")
        ax.add_collection(line_segments)
        ax.plot(mean_trace.q,mean_trace.SA,color="blue",alpha=1.0,label="mean_postfilter")
        ax.fill_between(mean_trace.q, mean_trace.SA+3*mean_trace.sigSA, mean_trace.SA-3*mean_trace.sigSA, facecolor='red', alpha=0.3,label=r"3 $\sigma$")
        # ax.plot(one_curve.q,clean[0],color="blue",alpha=1.0,label="mean_postfilter")
        ax.plot(mean_trace.q,median_SA,color="orange",label="median_prefilter")
        ax.set_ylabel(r'$\Delta I$ (A.U.)')
        ax.set_xlabel(r'q ($\AA^{-1}$)')
        ax.set_xscale("log")
        ax.set_title("{}_{}".format(samp,key))
        ax.legend()
        ### faster plotting
        # ax.set_xscale("log")
        # plt.show()
        fig.tight_layout()
        fig.savefig("{}_{}_cleaned_scaled_TR.png".format(samp,key))
        ax.cla()

    return


def plot_mean_TR_traces(samp, trace_lib):

    fig, ax = plt.subplots(figsize=(6, 4),dpi=300)
    ii = -1
    for key in trace_lib.keys():
        if ii < 0:
            ii=0
            nii = ii
        else:
            N = plt.cm.inferno.N
            ii += int(N/len(trace_lib.keys()))
            nii = N-ii
        print("====")
        print("plotting mean trace for {}".format(key))
        mean_trace = average_traces(trace_lib[key])
        ax.plot(mean_trace.q, mean_trace.SA, label="{}".format(key), color=plt.cm.inferno(nii))
        ax.fill_between(mean_trace.q, mean_trace.SA+3*mean_trace.sigSA, mean_trace.SA-3*mean_trace.sigSA, facecolor=plt.cm.inferno(nii), alpha=0.3)
    handles, labels = fig.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.set_title('{} - Time Resolved Signal'.format(samp))
    ax.set_ylabel(r'$\Delta I$ (A.U.)')
    ax.set_xlabel(r'q ($\AA^{-1}$)')
    ax.set_xscale("log")
    fig.legend(by_label.values(), by_label.keys())
    fig.tight_layout()
    fig.savefig("{}_scaled_TR".format(samp))

    return


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
filtered_vectors = {key:iterative_chi_filter(subtracted_vectors[key]) for key in subtracted_vectors.keys()}

all_off_vectors = all_off_traces(parent, samp, reps, on_off_map)
filtered_off_vectors = iterative_chi_filter(all_off_vectors)

buffer_dir = "/Volumes/beryllium/saxs_waxs_tjump/cypa/APS_20170302/CypA-WT-Buffer-1/xray_images/"
parent2, samp2, reps2, on_off_map2 = sample_map(buffer_dir)
buffer_all_off_vectors = all_off_traces(parent2, samp2, reps2, on_off_map2)
buffer_filtered_off_vectors = iterative_chi_filter(buffer_all_off_vectors)


avg_filt_off = average_traces(filtered_off_vectors)
avg_filt_off.set_name("testing_testing_one_two")
print("mic check:\n{}".format(avg_filt_off.name))
buff_avg_filt_off = average_traces(buffer_filtered_off_vectors)
protein_only_avg_filt_off = subtract_unscaled_traces(avg_filt_off,buff_avg_filt_off)
mean_TR = [average_traces(filtered_vectors[key]) for key in filtered_vectors.keys()]
showme = [add_unscaled_traces(protein_only_avg_filt_off,item) for item in mean_TR]
showme[0].write_dat("first_output.dat")

real_space_plotter(showme)


# plot_mean_TR_traces(samp,filtered_vectors)



# plot_all_traces(samp,subtracted_vectors)
# for key in filtered_vectors.keys():
#     print("{} has {} vectors".format(key, len(filtered_vectors[key])))




# t1 = clock()
# print("TR signal scaled & plotted for {:20s} took {:4.2f} seconds".format(samp, t1-t0))


# # matrix = np.matrix([i[0][6:] for i in vectors]).transpose()
# # # matrix = np.matrix(subbed_vectors).transpose()

# # u,s,v = svd(matrix, full_matrices=False)


# # fig3, ax3 = plt.subplots()
# # i = 0
# # for vector in v.tolist()[0:5]:
# #     ax3.plot([value+i*0.3 for value in vector], "-")
# #     i+=1
# # fig3.savefig("timepoints.png")

# # fig4, ax4 = plt.subplots()
# # j = 0

# # for vector in u.transpose().tolist()[0:5]:
# #     x = on_data.q[6:]    
# #     ax4.plot(x, [value for value in vector], "-")
# #     j+=0.3
# # ax4.set_xscale("log")
# # fig4.savefig("Singular_vectors_WT_HD_Unscaled.png")

# # zero_vector_values_on = v.tolist()[0][1::2]
# # zero_vector_values_off = v.tolist()[0][0::2]
# # average_on = np.average(zero_vector_values_on)
# # std_on = np.std(zero_vector_values_on)
# # average_off =  np.average(zero_vector_values_off)
# # std_off = np.std(zero_vector_values_off)
# # print(average_on, average_off, std_on, std_off)
# # for index, value in enumerate(v.tolist()[0]):
# #     if index % 2 == 0:
# #         if np.fabs((value - average_off) / std_off) > 2.5:
# #             print(vectors[index][1])
# #             print(vectors[index][2])
# #     elif index % 2 == 1:
# #         if np.fabs((value - average_on) / std_on) > 2.5:
# #             print(vectors[index][1])
# #             print(vectors[index][2])


    
# # fig5, ax5= plt.subplots()
# # ax5.plot([np.log(i) for i in s][0:10], "-")
# # fig5.savefig("Singular_values_WT_HD_Unscaled.png")

# # fig6, ax6 = plt.subplots()
# # ax6.axhline(average_on+2.5*std_on)
# # ax6.axhline(average_on-2.5*std_on)
# # ax6.axhline(average_on, color="0.5")
# # ax6.plot(zero_vector_values_on)
# # fig6.savefig("On_vectors.png")
# # ax6.cla()
# # ax6.axhline(average_off+2.5*std_off)
# # ax6.axhline(average_off-2.5*std_off)
# # ax6.axhline(average_off, color="0.5")
# # ax6.plot(zero_vector_values_off)
# # fig6.savefig("Off_vectors.png")
# # ax6.cla()
# # ax6.axhline(average_off+2.5*std_off)
# # ax6.axhline(average_off-2.5*std_off)
# # ax6.axhline(average_off, color="0.5")
# # ax6.plot(v.tolist()[0])
# # fig6.savefig("All_vectors_off_lines.png")


# from new_analysis import real_space_plotter
# # mean_traces = [average_traces(filtered_vectors[key]) for key in filtered_vectors.keys()]
# mean_trace = [average_traces(filtered_vectors['422us'])]
# real_space_plotter(mean_trace)


