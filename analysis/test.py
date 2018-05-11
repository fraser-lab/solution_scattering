from numpy import linalg as LA
import numpy as np
from trace import Trace
import parse
import matplotlib.pyplot as plt
plt.style.use('ggplot')



reference = parse.parse("/Volumes/beryllium/saxs_waxs_tjump/cypa/APS_20170302/CypA-WT-1/xray_images/CypA-WT-1_9_4.22us.tpkl")
test_one = parse.parse("/Volumes/beryllium/saxs_waxs_tjump/cypa/APS_20170302/CypA-WT-1/xray_images/CypA-WT-1_18_2.37us.tpkl")
test_two = parse.parse("/Volumes/beryllium/saxs_waxs_tjump/cypa/APS_20170302/CypA-WT-1/xray_images/CypA-WT-1_18_2.37us.tpkl")
test_three = parse.parse("/Volumes/beryllium/saxs_waxs_tjump/cypa/APS_20170302/CypA-WT-1/xray_images/CypA-WT-1_18_2.37us.tpkl")


print("test_one_SA equals test_two_SA : {}".format(test_one.SA == test_two.SA))


test_one.alg_scale(reference)
test_two.projection_scale(reference)
test_three.integration_scale(reference)


print("test_one_scaled_SA equals test_two_scaled_SA : {}".format(test_one.scaled_SA == test_two.scaled_SA))
print("test_one_scaled_SA equals test_three_scaled_SA : {}".format(test_one.scaled_SA == test_three.scaled_SA))
print("test_two_scaled_SA equals test_three_scaled_SA : {}".format(test_two.scaled_SA == test_three.scaled_SA))

print("alg_scale_sf : {}".format(test_one.scale_factor))
print("projection_scale_sf : {}".format(test_two.scale_factor))
print("integration_scale_sf : {}".format(test_three.scale_factor))



plt.plot(test_one.q,test_one.scaled_SA, lw=3, color='blue', label="alg_scaled")
plt.plot(test_two.q,test_two.scaled_SA, lw=2, color='red', label="projection_scaled")
plt.plot(test_three.q,test_three.scaled_SA, lw=1, color='green', label="integration_scaled")
plt.xscale('log')
plt.legend()
plt.show()








