### use python 3 error

import sys
if sys.version_info[0] < 3:
    print("\nException: Please use Python 3")
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
from saxs_plots import real_space_plotter
import pickle
from argparse import ArgumentParser, RawDescriptionHelpFormatter





def sample_map_multitemp(samp_dir, multitemp=None, low_cutoff=0, high_cutoff=1000):
    samp_dir = pathlib.Path(samp_dir)
    samp_files = list(samp_dir.glob(pattern='**/*.tpkl'))
    # buffer_files = list(buff.glob(pattern='**/*.tpkl'))
    t0 = clock()
    REPS = []
    TIMES = []
    ITERATIONS = []
    TEMPS = []
    for file in samp_files:
        name = file.name
        parent = file.parent
        if multitemp:
            samp, iteration, temp, rep, time = name.split('_')
            time = time.replace('.tpkl','')
            ITERATIONS.append(iteration)
            TEMPS.append(temp)
        else:
            try:
                samp, rep, time, onoff = name.split('_')
            except:
                samp, rep, time = name.split('_')
                time = time.replace('.tpkl','')
        REPS.append(rep)
        TIMES.append(time)
        
    if multitemp:
        ITERATIONS = sorted(list(set(ITERATIONS)), key=float)
        TEMPS = list(set(TEMPS))
    else:
        pass
    REPS = sorted(list(set(REPS)), key=float)
    REPS = [x for x in REPS if int(x) > low_cutoff]
    REPS = [x for x in REPS if int(x) < high_cutoff]
    TIMES = list(set(TIMES))
    OFFS = [item for item in TIMES if "-10us" in item]
    ONS = [item for item in TIMES if "-10us" not in item]
    # OFFS = [item for item in TIMES if "off" in item]
    # ONS = [item for item in TIMES if "on" not in item]

    tup =  [parse.unit_sort(item) for item in ONS]
    tup_sort = sorted(tup, key=lambda item: (item[0],parse.natural_keys(item[1])))
    clean_ons = [item[1] for item in tup_sort]
    on_off_map = {k:v for k,v in (zip(clean_ons,sorted(OFFS, key=parse.alphanum_key)))}
    # on_off_map = {k:v for k,v in (zip(clean_ons,sorted(OFFS, key=parse.alphanum_key)))}
    return parent, samp, ITERATIONS, TEMPS, REPS, on_off_map

def sample_map(file_dir, low_cutoff=0, high_cutoff=1000):
    file_dir = pathlib.Path(file_dir)
    files = list(file_dir.glob(pattern='*/*.tpkl'))
    # buffer_files = list(buff.glob(pattern='**/*.tpkl'))
    t0 = clock()
    REPS = []
    TIMES = []
    for file in files:
        name = file.name
        parent = file.parent
        try:
            samp, rep, time, onoff = name.split('_')
        except:
            samp, rep, time = name.split('_')
            time = time.replace('.tpkl','')
        # samp, rep, time = name.split('_')
        # time = time.replace('.tpkl','')
        REPS.append(rep)
        TIMES.append(time)

    REPS = sorted(list(set(REPS)), key=float)
    REPS = [x for x in REPS if int(x) > low_cutoff]
    REPS = [x for x in REPS if int(x) < high_cutoff]
    TIMES = list(set(TIMES))
    OFFS = [item for item in TIMES if "-10us" in item]
    ONS = [item for item in TIMES if "-10us" not in item]

    tup =  [parse.unit_sort(item) for item in ONS]
    tup_sort = sorted(tup, key=lambda item: (item[0],parse.natural_keys(item[1])))
    clean_ons = [item[1] for item in tup_sort]
    on_off_map = {k:v for k,v in (zip(clean_ons,sorted(OFFS, key=parse.alphanum_key)))}

    return parent, samp, REPS, on_off_map

def static_map(samp_dir, buffer_d=None, low_cutoff=0, high_cutoff=1000):
    samp_dir = pathlib.Path(samp_dir)
    samp_files = list(samp_dir.glob(pattern='**/*.tpkl'))
    # buffer_files = list(buff.glob(pattern='**/*.tpkl'))
    # t0 = clock()
    # sample_map = []
    REPS = []
    TEMPS = []
    SERIES = []
    # TIMES = []
    if buffer_d:
        samp_files = [file for file in samp_files if "BT" in file.name]
    else:
        samp_files = [file for file in samp_files if "BT" not in file.name]
    for file in samp_files:
        name = file.name
        parent = file.parent
        samp, info, rep = name.split('_')
        if "BT" in info:
            dilution = info[3:4]
            temp = info[5:]
            # print(temp)
        elif "PC" in info:
            dilution = info[3:6]
            temp = info[7:]
            # print(temp)
        else:
            print(info+"failed")
        rep = rep.replace(".tpkl","")
        # on_data = parse.parse(str(file))
        # on_data.scale(reference)
        # sample_map.append((samp,dilution,temp,on_data))

        # time = time.replace('.tpkl','')
        REPS.append(rep)
        SERIES.append(dilution)
        TEMPS.append(temp)
        # TIMES.append(time)

    REPS = sorted(list(set(REPS)), key=float)
    REPS = [x for x in REPS if int(x) > low_cutoff]
    REPS = [x for x in REPS if int(x) < high_cutoff]
    TEMPS = sorted(list(set(TEMPS)), key=float)
    SERIES = sorted(list(set(SERIES)))
    # TIMES = list(set(TIMES))
    # OFFS = [item for item in TIMES if "-10us" in item]
    # ONS = [item for item in TIMES if "-10us" not in item]
    # OFFS = [item for item in TIMES if "off" in item]
    # ONS = [item for item in TIMES if "on" not in item]

    # tup =  [parse.unit_sort(item) for item in ONS]
    # tup_sort = sorted(tup, key=lambda item: (item[0],parse.natural_keys(item[1])))
    # clean_ons = [item[1] for item in tup_sort]
    # on_off_map = {k:v for k,v in (zip(clean_ons,clean_ons))}

    return parent, samp, REPS, TEMPS, SERIES

    # return sample_map


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


def time_resolved_traces(parent, samp, reps, on_off_map, option=None, multitemp=None, iterations=None, temp=None):
    
    subtracted_vectors = {i: [] for i in on_off_map.keys()}
    
    if multitemp:
        for iteration in iterations:
            for n in reps:
                for on, off in on_off_map.items():

                    # on_string = ("{0}/{1}_{2}_{3}_on.tpkl".format(parent, samp, n, on))
                    on_string = ("{0}/{1}_{2}_{3}_{4}_{5}.tpkl".format(parent, samp, iteration, temp, n, on))
                    # if option:
                    #     on_string = ("{0}/{1}_{2}_{3}_on.tpkl".format(parent, samp, n, on))
                    # else:
                    #     on_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, on))

                    # if multitemp:
                    #     on_string = ("{0}/{1}_{2}_{3}_{4}_{5}.tpkl".format(parent, samp, iteration, temp, n, on))
                    try:
                        on_data = parse.parse(on_string)
                        on_data.scale(reference, qmin=QMIN, qmax=QMAX)
                    except:
                        print(on_string+"\tfailed")
                        pass

                    # off_string = ("{0}/{1}_{2}_{3}_off.tpkl".format(parent, samp, n, off))
                    off_string = ("{0}/{1}_{2}_{3}_{4}_{5}.tpkl".format(parent, samp, iteration, temp, n, off))
                    # if option:
                    #     off_string = ("{0}/{1}_{2}_{3}_on.tpkl".format(parent, samp, n, off))
                    # else:
                    #     off_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, off))
                    
                    # if multitemp:
                    #     off_string = ("{0}/{1}_{2}_{3}_{4}_{5}.tpkl".format(parent, samp, iteration, temp, n, off))
                    
                    try:
                        off_data = parse.parse(off_string)
                        off_data.scale(reference, qmin=QMIN, qmax=QMAX)
                    except:
                        print(off_string+"\tfailed")
                        pass

                    if isinstance(on_data,type(None)) or isinstance(off_data,type(None)):
                        print(on_string+"\tfailed")
                        pass
                    else:
                        # sub_scaled = subtract_scaled_traces(on_data,off_data)
                        sub_scaled = on_data.subtract(off_data, scaled=True)
                        subtracted_vectors[on].append(sub_scaled)

    else:
        for n in reps:
            for on, off in on_off_map.items():
                on_data = None
                off_data = None

                if option:
                    on_string = ("{0}/{1}_{2}_{3}_on.tpkl".format(parent, samp, n, on))
                    off_string = ("{0}/{1}_{2}_{3}_on.tpkl".format(parent, samp, n, off))
                else:
                    on_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, on))
                    off_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, off))

                try:
                    on_data = parse.parse(on_string)
                    on_data.scale(reference, qmin=QMIN, qmax=QMAX)
                except:
                    print(on_string+"\tfailed")
                    pass
                
                try:
                    off_data = parse.parse(off_string)
                    off_data.scale(reference, qmin=QMIN, qmax=QMAX)
                except:
                    print(off_string+"\tfailed")
                    pass

                if isinstance(on_data,type(None)) or isinstance(off_data,type(None)):
                    print(on_string+"\tfailed")
                    pass
                else:
                    # sub_scaled = subtract_scaled_traces(on_data,off_data)
                    sub_scaled = on_data.subtract(off_data, scaled=True)
                    subtracted_vectors[on].append(sub_scaled)

    return subtracted_vectors


def static_traces(parent, samp, reps, temps, series, option=None):
    
    static_vectors = {i: {j: [] for j in series} for i in temps}

    for temp in temps:
        
        for dilution in series:
            static = []
            for n in reps:
                static_string = ("{0}/{1}_off{2}T{3}_{4}.tpkl".format(parent, samp, dilution, temp, n))
                # print(static_string)

                try:
                    static_data = parse.parse(static_string)
                    # print("test1")
                    # static_data.scale(reference, qmin=QMIN, qmax=QMAX)
                    static_data.scale(reference, qmin=QMIN, qmax=QMAX, approach='algebraic')
                    # print("test2")
                    static_scaled = Trace(static_data.q, np.empty_like(static_data.q), np.empty_like(static_data.q), static_data.scaled_sigSA, static_data.scaled_SA, static_data.Nj)
                    static.append(static_scaled)
                except:
                    # print(buff_string+"\tfailed")
                    print("{} failed to parse or scale".format(static_string))
                    pass
            try:
                static_filtered = iterative_chi_filter(static)
                static_filt_avg = average_traces(static_filtered)
                static_vectors[temp][dilution].append(static_filt_avg)
            except:
                print("temp {}C failed for {}".format(temp,dilution))
                
    return static_vectors

def all_off_traces(parent, samp, reps, on_off_map, option=None, multitemp=None, iterations=None, temp=None):
    
    off_vectors = []

    if multitemp:
        for iteration in iterations:

            for n in reps:
                for off in on_off_map.values():

                    off_string = ("{0}/{1}_{2}_{3}_{4}_{5}.tpkl".format(parent, samp, iteration, temp, n, off))
                    try:
                        off_data = parse.parse(off_string)
                        off_data.scale(reference, qmin=QMIN, qmax=QMAX)
                        off_scaled = Trace(off_data.q, np.empty_like(off_data.q), np.empty_like(off_data.q), off_data.scaled_sigSA, off_data.scaled_SA, off_data.Nj)
                        off_vectors.append(off_scaled)
                    except:
                        print(off_string+"\tfailed")

    else:
        for n in reps:
                    for off in on_off_map.values():

                        # off_string = ("{0}/{1}_{2}_{3}_off.tpkl".format(parent, samp, n, off))
                        if option:
                            off_string = ("{0}/{1}_{2}_{3}_on.tpkl".format(parent, samp, n, off))
                        else:
                            off_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, off))
                        try:
                            off_data = parse.parse(off_string)
                            off_data.scale(reference, qmin=QMIN, qmax=QMAX)
                            off_scaled = Trace(off_data.q, np.empty_like(off_data.q), np.empty_like(off_data.q), off_data.scaled_sigSA, off_data.scaled_SA, off_data.Nj)
                            off_vectors.append(off_scaled)
                        except:
                            print(off_string+"\tfailed")

    return off_vectors

def all_vectors(parent, samp, reps, on_off_map, option=None, multitemp=None, iterations=None, temp=None):
    
    all_vectors = []
    all_labels = []
    tr_vectors_labels = []


    if multitemp:
        for iteration in iterations:

            for n in reps:
                # for off in on_off_map.values():
                for on, off in on_off_map.items():
                    on_string = ("{0}/{1}_{2}_{3}_{4}_{5}.tpkl".format(parent, samp, iteration, temp, n, on))
                    off_string = ("{0}/{1}_{2}_{3}_{4}_{5}.tpkl".format(parent, samp, iteration, temp, n, off))
                    try:
                        on_data = parse.parse(on_string)
                        on_data = Trace(on_data.q, np.empty_like(on_data.q), np.empty_like(on_data.q), on_data.sigSA, on_data.SA, on_data.Nj)
                        on_data.scale(reference, qmin=QMIN, qmax=QMAX)

                        off_data = parse.parse(off_string)
                        off_data = Trace(off_data.q, np.empty_like(off_data.q), np.empty_like(off_data.q), off_data.sigSA, off_data.SA, off_data.Nj)
                        off_data.scale(reference, qmin=QMIN, qmax=QMAX)
                        # off_scaled = Trace(off_data.q, np.empty_like(off_data.q), np.empty_like(off_data.q), off_data.scaled_sigSA, off_data.scaled_SA, off_data.Nj)
                        # off_vectors.append(off_scaled)
                        if on_data:
                            if off_data:
                                all_vectors.append(off_data.SA[reference.q>0.03])
                                all_labels.append(off)
                                all_vectors.append(on_data.SA[reference.q>0.03])
                                all_labels.append(on)
                                # sub_scaled = subtract_scaled_traces(on_data,off_data)
                                sub_scaled = on_data.subtract(off_data, scaled=True)
                                tr_vectors_labels.append((sub_scaled.SA[reference.q>0.03], on))


                    except:
                        print(off_string+"\tfailed")

    # else:
    #     for n in reps:
    #                 for off in on_off_map.values():

    #                     # off_string = ("{0}/{1}_{2}_{3}_off.tpkl".format(parent, samp, n, off))
    #                     if option:
    #                         off_string = ("{0}/{1}_{2}_{3}_on.tpkl".format(parent, samp, n, off))
    #                     else:
    #                         off_string = ("{0}/{1}_{2}_{3}.tpkl".format(parent, samp, n, off))
    #                     try:
    #                         off_data = parse.parse(off_string)
    #                         off_data.scale(reference)
    #                         off_scaled = Trace(off_data.q, np.empty_like(off_data.q), np.empty_like(off_data.q), off_data.scaled_sigSA, off_data.scaled_SA, off_data.Nj)
    #                         off_vectors.append(off_scaled)
    #                     except:
    #                         print(off_string+"\tfailed")

    return all_vectors, all_labels, tr_vectors_labels


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

################################################################################

# if args.sample_directory:
#     data_dir = args.sample_directory
# elif args.buffer_directory:
#     data_dir = args.buffer_directory
# else:
#     print("\nException: No data provided")
#     sys.exit(1)



### Setup command-line input flags and help messages

# def main():
parser = ArgumentParser(usage='python3 reduce_data.py [options]', description='Analyze time-resolved solution scattering data.', formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-tr', '--time_resolved_differences', help='create average traces for each time point, with matching buffer subtraction', action='store_true')
parser.add_argument('-svd', help='singular value decomposition of a series of traces', action='store_true')
parser.add_argument('-mt', '--multitemp', help='Multiple temperatures within each sample directory', action='store_true')
parser.add_argument('-s', '--sample_directory', help='location of files to be analyzed')
parser.add_argument('-b', '--buffer_directory', help='location of matching buffer files for TR analysis')
parser.add_argument('-st', '--static_directory', help='location of static files')
parser.add_argument('-r', '--reference', help='Provide a reference image for scaling. If no reference is provided an arbitrary image will be chosen to scale all others against.')
parser.add_argument('-el', '--exclude_repeats_low', help='Choose repeats at the start of a dataset to exclude based on fast radiation damage effects', type=int, default=0)
parser.add_argument('-eh', '--exclude_repeats_high', help='Choose repeats at the end of a dataset to exclude based on slow radiation damage effects', type=int, default=10000)
args = parser.parse_args()


reference = parse.parse("/Volumes/beryllium/saxs_waxs_tjump/cypa/APS_20170302/CypA-WT-1/xray_images/CypA-WT-1_9_4.22us.tpkl")
if args.reference:
    reference = parse.parse(args.reference)
else:
    pass


QMIN = 0.0025
# QMAX = 5.1925
QMAX = 4.28
CHI_OUTLIER = 1.5
t_shortlist = ["-10.1us", "1us", "10us", "100us", "1ms"]
# t_shortlist = ["562ns"]
# t_shortlist = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us"]
TIMES = [ "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us"]
temp='3C'

################################################################################
###Begin time-resolved differences on scaled frames###


if args.time_resolved_differences:



    if args.multitemp:
        parent, samp, iterations, temps, reps, on_off_map = sample_map_multitemp(args.sample_directory, multitemp=args.multitemp, low_cutoff=args.exclude_repeats_low, high_cutoff=args.exclude_repeats_high)
        filtered_off_vectors = {}
        filtered_vectors = {}
        for temp in temps:
            subtracted_vectors = time_resolved_traces(parent, samp, reps, on_off_map, multitemp=args.multitemp, iterations=iterations, temp=temp)
            filtered_vectors[temp] = {key:iterative_chi_filter(subtracted_vectors[key]) for key in subtracted_vectors.keys()}
            all_off_vectors = all_off_traces(parent, samp, reps, on_off_map, multitemp=args.multitemp, iterations=iterations, temp=temp)
            filtered_off_vectors[temp] = iterative_chi_filter(all_off_vectors)

        with open("filtered_off_vectors_dict.pkl", "wb") as pkl:
            pickle.dump(filtered_off_vectors, pkl)
        with open("filtered_vectors_dict.pkl", "wb") as pkl:
            pickle.dump(filtered_vectors, pkl)

        parent2, samp2, iterations2, temps2, reps2, on_off_map2 = sample_map_multitemp(args.buffer_directory, multitemp=args.multitemp, low_cutoff=args.exclude_repeats_low, high_cutoff=args.exclude_repeats_high)
        buffer_filtered_off_vectors = {}
        buffer_filtered_vectors = {}

        for temp2 in temps2:
            buffer_TR_subtracted_vectors = time_resolved_traces(parent2, samp2, reps2, on_off_map2, multitemp=args.multitemp, iterations=iterations2, temp=temp2)
            buffer_filtered_vectors[temp2] = {key:iterative_chi_filter(buffer_TR_subtracted_vectors[key]) for key in buffer_TR_subtracted_vectors.keys()}
            buffer_all_off_vectors = all_off_traces(parent2, samp2, reps2, on_off_map2, multitemp=args.multitemp, iterations=iterations2, temp=temp2)


        with open("buffer_filtered_off_vectors_dict.pkl", "wb") as pkl:
            pickle.dump(buffer_filtered_off_vectors, pkl)
        with open("buffer_filtered_vectors_dict.pkl", "wb") as pkl:
            pickle.dump(buffer_filtered_vectors, pkl)



        avg_filt_off = {temp:average_traces(filtered_off_vectors[temp]) for temp in filtered_off_vectors.keys()}
        for key in avg_filt_off.keys():
            avg_filt_off[key].buffer_scale(avg_filt_off[key])
        with open("avg_filt_off_dict.pkl", "wb") as pkl:
            pickle.dump(avg_filt_off, pkl)

        buff_avg_filt_off = {temp:average_traces(buffer_filtered_off_vectors[temp]) for temp in buffer_filtered_off_vectors.keys()}
        for key in buff_avg_filt_off.keys():
            buff_avg_filt_off[key].buffer_scale(avg_filt_off[key])
        with open("buff_avg_filt_off_dict.pkl", "wb") as pkl:
            pickle.dump(buff_avg_filt_off, pkl)


        # protein_only_avg_filt_off = {temp:buffer_subtract_scaled_traces(avg_filt_off[temp],buff_avg_filt_off[temp]) for temp in filtered_off_vectors.keys()}
        protein_only_avg_filt_off = {temp:avg_filt_off[temp].subtract(buff_avg_filt_off[temp], buffer_scaled=True) for temp in filtered_off_vectors.keys()}

        with open("protein_only_avg_filt_off_dict.pkl", "wb") as pkl:
            pickle.dump(protein_only_avg_filt_off, pkl)


        # mean_TR = {temp:{key: subtract_unscaled_traces(average_traces(filtered_vectors[temp][key]),average_traces(buffer_filtered_vectors[temp][key])) for key in filtered_vectors[temp].keys()} for temp in filtered_off_vectors.keys()}
        mean_TR = {temp:{key: average_traces(filtered_vectors[temp][key]).subtract(average_traces(buffer_filtered_vectors[temp][key])) for key in filtered_vectors[temp].keys()} for temp in filtered_off_vectors.keys()}
        for temp in mean_TR.keys():
            for  diff in mean_TR[temp].keys():
                mean_TR[temp][diff].write_dat(samp+"_"+temp+"_diff_"+diff+".dat")

        with open("mean_TR_dict.pkl", "wb") as pkl:
            pickle.dump(mean_TR, pkl)

        # showme = {temp:{key: add_unscaled_traces(protein_only_avg_filt_off[temp],mean_TR[temp][key]) for key in mean_TR[temp].keys()} for temp in filtered_off_vectors.keys()}
        showme = {temp:{key: protein_only_avg_filt_off[temp].add(mean_TR[temp][key]) for key in mean_TR[temp].keys()} for temp in filtered_off_vectors.keys()}

        for temp in showme.keys():
            for itm in showme[temp].keys():
                showme[temp][itm].write_dat(samp+"_"+temp+"_"+itm+".dat")

        with open("TR-plus-avg_dict.pkl", "wb") as pkl:
            pickle.dump(showme, pkl)

        for temp in showme.keys():
            real_space_plotter(showme[temp],temp)


    else:

        parent, samp, reps, on_off_map = sample_map(args.sample_directory, low_cutoff=args.exclude_repeats_low, high_cutoff=args.exclude_repeats_high)
        subtracted_vectors = time_resolved_traces(parent, samp, reps, on_off_map, option=False)
        filtered_vectors = {key:iterative_chi_filter(subtracted_vectors[key]) for key in subtracted_vectors.keys()}

        all_off_vectors = all_off_traces(parent, samp, reps, on_off_map, option=False)
        filtered_off_vectors = iterative_chi_filter(all_off_vectors)

        parent2, samp2, reps2, on_off_map2 = sample_map(args.buffer_directory, low_cutoff=args.exclude_repeats_low, high_cutoff=args.exclude_repeats_high)
        buffer_TR_subtracted_vectors = time_resolved_traces(parent2, samp2, reps2, on_off_map2, option=False)
        buffer_filtered_vectors = {key:iterative_chi_filter(buffer_TR_subtracted_vectors[key]) for key in buffer_TR_subtracted_vectors.keys()}
        buffer_all_off_vectors = all_off_traces(parent2, samp2, reps2, on_off_map2, option=False)
        buffer_filtered_off_vectors = iterative_chi_filter(buffer_all_off_vectors)



        avg_filt_off = average_traces(filtered_off_vectors)
        avg_filt_off.buffer_scale(avg_filt_off)
        avg_filt_off.write_dat(samp+"_average_off_filtered_"+".dat")
        buff_avg_filt_off = average_traces(buffer_filtered_off_vectors)
        buff_avg_filt_off.buffer_scale(avg_filt_off)
        buff_avg_filt_off.write_dat(samp+"_buffer_average_off_filtered_"+".dat")
        protein_only_avg_filt_off = avg_filt_off.subtract(buff_avg_filt_off, buffer_scaled=True)
        protein_only_avg_filt_off.write_dat(samp+"_protein_only_average_off_filtered_"+".dat")
        # protein_only_avg_filt_off.buffer_scale(protein_only_avg_filt_off)
        # mean_TR = {key: subtract_unscaled_traces(average_traces(filtered_vectors[key]),average_traces(buffer_filtered_vectors[key])) for key in filtered_vectors.keys()}
        mean_TR = {key: average_traces(filtered_vectors[key]).subtract(average_traces(buffer_filtered_vectors[key])) for key in filtered_vectors.keys()}
        for diff in mean_TR.keys():
            mean_TR[diff].write_dat(samp+"_diff_"+diff+".dat")
        #     value.buffer_scale(protein_only_avg_filt_off)
        showme = {key: protein_only_avg_filt_off.add(mean_TR[key]) for key in mean_TR.keys()}
        # showme["750ns"].write_dat("first_output.dat")
        # print(showme.keys())
        for itm in showme.keys():
            showme[itm].write_dat(samp+"_sum_"+itm+".dat")

        # real_space_plotter(showme, name=samp)




###End time-resolved differences on scaled frames###
################################################################################
###Begin SVD on scaled frames###

elif args.svd:

    parent, samp, iterations, temps, reps, on_off_map = sample_map_multitemp(args.sample_directory, multitemp=args.multitemp)
    multi_all_vectors = {}
    multi_all_labels = {}
    tr_vectors_labels = {}
    # filtered_vectors = {}
    if args.multitemp:
        for temp in temps:
            multi_all_vectors[temp], multi_all_labels[temp], tr_vectors_labels[temp] = all_vectors(parent, samp, reps, on_off_map, multitemp=args.multitemp, iterations=iterations, temp=temp)
        full_list = [item[0] for item in tr_vectors_labels['3C']]
        full_labels = [item[1] for item in tr_vectors_labels['3C']]
    else:
        multi_all_vectors, multi_all_labels, tr_vectors_labels = all_vectors(parent, samp, reps, on_off_map, multitemp=args.multitemp)

        full_list = [item[0] for item in tr_vectors_labels]
        full_labels = [item[1] for item in tr_vectors_labels]

    # print(full_list)
    fig0, ax0 = plt.subplots()
    all_curves = [list(zip(reference.q[reference.q>0.03],item)) for item in full_list]
    # print(all_curves[0])
    line_segments = LineCollection(all_curves,color="green",lw=0.5)
    ax0.add_collection(line_segments)
    ax0.plot()
    # ax0.scatter(reference.q,full_list)
    fig0.savefig("scaled_radavgs.png", dpi=300)
    plt.show(block=False)


    u,s,v = svd(full_list, full_matrices=False)
    with open("svd_u.pkl", "wb") as pkl:
        pickle.dump(u, pkl)
    with open("svd_s.pkl", "wb") as pkl:
        pickle.dump(s, pkl)
    with open("svd_v.pkl", "wb") as pkl:
        pickle.dump(v, pkl)

    fig, ax = plt.subplots()
    i = 0
    #print("xx shape = {}".format(xx.shape))
    for vector in v[0:8]:
        # print vector
        #print("vector shape = {}".format(vector.shape))
        # ax.plot(range(len(vectors)), [value+i for value in vector], "-")
        ax.plot(reference.q[reference.q>0.03], vector+i*0.1, "-", label = "v{}".format(i), lw=0.7)
        i+=1
    plt.legend()
    ax.set_xscale('log')
    #fig.savefig("{}_svd.png".format(run_numb))
    fig.savefig("singular_vectors.png", dpi=300)
    plt.show(block=False)

    #np.save("time_dep_vector", v[2])

    fig2, ax2 = plt.subplots()
    i = 0
    for vector in u.T[0:4]:
        # print vector
        # ax.plot(range(len(vectors)), [value+i for value in vector], "-")
        # x = [i*0.025 for i in range(len(vector))] 
        ax2.plot(vector, label = "v{}".format(i), lw=0.5)
        # ax2.scatter(light, vector[light]+i*.3, marker='+', color='red', edgecolors="none", s=2, label = "v{} light".format(i))
        i+=1
    #plt.legend()
        
    #fig2.savefig("{}_result.png".format(run_numb))
    fig2.savefig("vector_per_image.png", dpi=300)
    fig2.set_figwidth(15)
    fig3, ax3= plt.subplots()
    ax3.plot([np.log(i) for i in s][0:8], "-")
    #fig3.savefig("{}_singular_values.png".format(run_numb))
    fig3.savefig("singular_values.png")
    plt.show(block=False)
    #plt.show(figsize(10,10),dpi=300)
        # print i
        # print ordered_keylist
    #fig2.show()
    fig4, ax4 = plt.subplots()
    i = 0

    time_resolved_vectors = {}
    for vector in u.T[0:8]:
        for time in TIMES:
            time_resolved_vectors[time] = np.mean(np.array(vector[full_labels.index(time)]))
            # ax4.scatter(parse.times_numeric(time),np.mean(np.subtract(np.array(vector[full_labels.index(time)]), np.array(vector[full_labels.index(time)-1]))), color="black", label=time)

            ax4.scatter(parse.times_numeric(time),np.mean(np.array(vector[full_labels.index(time)])), color="black", label=time)

        # print vector
        # ax.plot(range(len(vectors)), [value+i for value in vector], "-")
        # x = [i*0.025 for i in range(len(vector))] 
            # ax4.hist(vector[full_labels.index(time)], 100, color='blue', alpha=0.5, label = "v{}".format(str(i)+'_'+str(time)))
        # ax4.hist(vector[light], 500, color='red', alpha=0.5, label = "v{} light".format(i))
        ax4.set_xscale('log')
        ax4.set_xlabel('Time (ns)')
        # plt.legend()
        fig4.savefig("v{}_time_dependence.png".format(i), dpi=300)
        # plt.show()
        
        i+=1
        ax4.cla()
    plt.show()






################################################################################
###End SVD on scaled frames###

elif args.static_directory:
    protein = static_map(args.static_directory, low_cutoff=args.exclude_repeats_low, high_cutoff=args.exclude_repeats_high)
    parent, samp, reps, temps, series = protein
    buff = static_map(args.buffer_directory, buffer_d=True, low_cutoff=args.exclude_repeats_low, high_cutoff=args.exclude_repeats_high)
    parent2, samp2, reps2, temps2, series2 = buff
    statics = static_traces(parent, samp, reps, temps, series)
    buffs = static_traces(parent2, samp2, reps2, temps2, series2)
    for temp in statics.keys():
        for dilution in statics[temp].keys():
            kicker = statics[temp][dilution][0].subtract(buffs[temp]['B'][0])
            kicker.write_dat(samp+"_static_"+temp+"C_"+dilution+"_filtered_buffersubtracted_average"+".dat")
            # statics[temp][dilution][0].write_dat(samp+"_static_"+temp+"C_"+dilution+"_filtered_average"+".dat")
            # print("{} and {} worked".format(temp,dilution))

else:
    pass

#     return



# if __name__ == "__main__":
#     main()

