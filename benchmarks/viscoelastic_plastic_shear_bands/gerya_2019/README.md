The files in this directory reproduce or analyze the results of
a benchmark for the shortening visco(elastic)-plastic block in
the absence of gravity. 

The benchmark is a modified version of Exercise 13.2 from   
Introduction to Numerical Geodynamic Modeling (Taras Gerya,
2019, doi:https://doi.org/10.1017/9781316534243).

Detailed file description:

- gerya_2019_vp.prm: Parameter file for the visco-plastic
version of the shortening benchmark.

- gerya_2019_vep.prm: Parameter file that is derived from 
gera_2019_vep.prm and modifies the rheological model to
include viscoelasticity (i.e., viscoelastic-plastic rheology).

- gerya_2019_analysis.py: Python script that loads in  
solution data and compares strain rate profiles.


 
