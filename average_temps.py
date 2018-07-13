import os
import sys
import numpy as np

input_filename = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])


cwd = os.getcwd()
input_file_fullpath = os.path.join(cwd, input_filename)
input_file = open(input_file_fullpath, 'r')

all_temps = []

for temp in input_file.readlines():
  all_temps.append(float(temp))

start_index = start-1
end_index = end

temp_array = np.asarray(all_temps[start_index:end_index])
avg = np.average(temp_array)
std = np.std(temp_array)

print "%8.4f" % (avg)
print "%8.4f" % (std)

