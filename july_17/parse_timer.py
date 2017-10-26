from sys import argv
from time import clock
from parse import parse


"""Statics"""
CHI_OUTLIER = 1.5
TEMPS = ["20C"]
TIMES = ["-10.1us", "562ns", "750ns", "1us", "1.33us", "1.78us", "2.37us", "3.16us", "4.22us", "5.62us", "7.5us", "10us", "13.3us", "17.8us", "23.7us", "31.6us", "42.2us", "56.2us", "75us", "100us", "133us", "178us", "237us", "316us", "422us", "562us", "750us", "1ms", '3.16ms', '10ms', '31.6ms', '100ms']
REPS = range(5,50)
PREFIX = "Trypsin-APMSF-1"
PKL_FILENAME = "Trypsin-APMSF-2_first-ten_full_algebraic_offs.pkl"
DATFILE_PREFIX = "Trypsin-APMSF-1"


directories = argv[1:]
reference = parse("/Volumes/beryllium/saxs_waxs_tjump/Trypsin/Trypsin-BA-Buffer-1/xray_images/Trypsin-BA-Buffer-1_26_-10us-10.tpkl")


# for directory in directories:
#     # files = listdir(directory)
#     # for index, _ in enumerate()
#     # print length
#     # for megarep in range(MEGAREPS):
#     t0 = clock()
#     for i in REPS: 
#         for temp in TEMPS:
#             for index, time in enumerate(TIMES):
#                 # try: 
#                     # onstring = subprocess.check_output("grep {0}_{1}_{2}_{3}_on beamstop-1.log".format(PREFIX, temp, i+1, time), shell=True)
#                     # onscale = int(onstring.split()[3])
#                     if index > 0:
#                         off_count = "-{}".format(index+1)
#                     else:
#                         off_count = ""

#                     on = parse("{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time))
                    
#                     on_scaled = on.alg_scale(reference)
#                     # ax1.plot(on.q, on_scaled[0],c='r')
#                     # on_scaled = lin_regress_scale(reference, on)
#                     # on_scaled = on.scale_isosbestic()[0]
#                     # on = parse_tpkl("{0}/{1}_{2}_{3}_{4}_on.tpkl".format(directory, PREFIX, temp, i+1, time)).as_vector()[80:]
                    
#                 # except:
#                 #     pass
#                 #     print "one or both of the on/off pairs was tossed:"
#                 #     print "{0}/{1}_{2}_{3}.tpkl".format(directory, PREFIX, i+1, time)
#     t1 = clock()
    # print("*"*80)

t2 = clock()
import pathlib
data_dir = pathlib.Path(directories[0])
files = list(data_dir.glob(pattern='**/*.tpkl'))
for file in files:
    on = parse(str(file))
    on.alg_scale(reference)

t3 = clock()
# print("\nparsing & scaling predefined took {} seconds".format(t1-t0))
print("parsing & scaling took {} seconds".format(t3-t2))




