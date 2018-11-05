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
difference_dat_kinetics.py         #tools to calculate kinetics from changes between time-resolved dat files
parse.py                           #library for reading SAXS data
quickplots.py                      #visual dat files from reduce_data
saxs_plots.py                      #plotting library for SAXS data
reduce_data.py                     #reduce replicates
structure_packing_calc.py          #calculate a SPF(q)
structure_packing_correction.py    #divide I(q) by SPF(q)
trace.py                           #object and methods for handling I(q)
```

### Usage

```
python3  reduce_data.py  -st dir_static_files  -b dir_static_buffer_files    #example of how to reduce data


python3  structure_packing_calc.py  dir_static_dat_files    #calculate SPF, write as dat


python3  structure_packing_correction.py  dir_uncorrected  SPF_dat_file    #divide I(q) by SPF(q)

```


