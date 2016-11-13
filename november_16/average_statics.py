import csv

import numpy as np

from parse import parse_tpkl, alg_scale



DIRECTORY_PREFIX = "/Users/benjaminbarad/Desktop/July_Trip_Data/CypA-WT/CypA-WT-static-1/xray_images/CypA-WT-static-1_off"
MUTANT = "WT"
VARS = ["B"]
REPS = range(33,65)
TEMPS = [21, 28]
CHI_OUTLIER = 1.5
reference = parse_tpkl("/Users/benjaminbarad/Desktop/July_Trip_Data/CypA-S99T/CypA-S99T-Buffer-2/xray_images/CypA-S99T-Buffer-2_17_-10us-14_on.tpkl")

def chisquared(var, ref):
  nu = len(var)-1.0
  I, sigma = zip(*var) 
  Iref,sigmaref = zip(*ref)
  chi_squared = 1/nu*sum([(I[i]-Iref[i])**2/(sigmaref[i]**2)  for i in range(len(var))])
  # print chi_squared
  return chi_squared

def chi_outliers(vectors, reference_vector):
  list = [chisquared(i, reference_vector) for i in vectors]
  print np.mean(list), np.std(list)
  outlier_list = []
  for i, val in enumerate(list):
    if val > CHI_OUTLIER:
      outlier_list.append(i)
      # print i, val
  return outlier_list

def combine_vectors(vectors):
  averaged_vector = []
  for i in range(len(vectors[0])):
    value_list = [v[i] for v in vectors]
    # averaged_vector.append(np.mean(value_list))
    #### Ensemble Weighting - mean = sum(value/sigma^2)/sum(1/sigma^2)
    ####                      std = sqrt(1/sum(1/sigma^2))
    means, stds = zip(*value_list)
    avg_mean = np.mean(means)
    std_mean = np.std(means)
    std_prop = np.sqrt(sum([j**2 for j in stds]))/(len(stds)-1)
    std_tot = np.sqrt(std_mean**2+std_prop**2)
    # avg_mean = sum([means[i]/(stds[i]**2) for i,_ in enumerate(means)])/sum([1/stds[i]**2 for i,_ in enumerate(means)])
    # std_tot = np.sqrt(1/sum([1/stds[i]**2 for i,_ in enumerate(means)]))
    # print on.q[i], avg_mean, std_mean, std_prop, std_tot
    averaged_vector.append((avg_mean, std_tot))
  outlier_list = chi_outliers(vectors, averaged_vector)
  if len(outlier_list) > 0:
    new_vectors = [vector for i,vector in enumerate(vectors) if i not in outlier_list]
    print len(new_vectors)
    averaged_vector = combine_vectors(new_vectors)
  return averaged_vector

def write_dat(averaged_data, filename):
  with open("static_dats/{}".format(filename), "w") as csvfile:
    writer = csv.writer(csvfile, delimiter=" ")
    for i, val in enumerate(averaged_data):
      writer.writerow([reference.q[i], val[0], val[1]])



for var in VARS:
  for temp in TEMPS:
    vectors = []
    for rep in REPS:
      vector = parse_tpkl("{}{}T{}_{}.tpkl".format(DIRECTORY_PREFIX, var, temp, rep))
      vector_scaled = alg_scale(reference, vector)
      vectors.append(zip(*vector_scaled))
    final_vector = combine_vectors(vectors)
    write_dat(final_vector, "{}_{}_{}.dat".format(MUTANT, var, temp))

