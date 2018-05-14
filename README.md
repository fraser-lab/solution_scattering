# Solution Scattering
------------------

## Summary:

Code for analyzing solution scattering data, particularly focusing on 
time-resolved experiments.

## Authors:
Ben Barad (bbarad)  
Michael Thompson (miket928)  
Alex Wolff (LifeHasOrder)  


### Contents

```
parse.py							#library of functions for parsing dat files, tpkl files, etc
quickplots.py				   		#input parameters for diffuse data processing
saxs_plots.py 						#plotting library for SAXS data
reduce_data.py					    #reduce replicates, carry out buffer subtraction, and much more
structure_packing_calc.py			#calculate a SPF vector for each temperature
structure_packing_correction.py		#divide scattering vector by corresponding SPF vector to correct for packing effects
trace.py							#trace objects are scattering curves, with internal methods for scaling and arithmetic
```

### Usage

```
python3    reduce_data.py    -st directory_with_static_files    -b directory_with_static_buffer_files     #example of how to reduce data


python3    structure_packing_calc.py    directory_with_static_dat_files		#calculate a structure packing factor and write as a dat file


python3    structure_packing_correction.py    directory_with_files_to_be_corrected    SPF_dat_file    #divide all dats by SPF

```


