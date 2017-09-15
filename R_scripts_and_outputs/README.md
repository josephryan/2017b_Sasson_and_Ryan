# R commands used in Sasson and Ryan, 2017

#### determine best fitting model w/likelihood ratio test
##### (to run: adjust mydir variable)
```Rscript likelihood_ratio_test.R > likelihood_ratio_test.R.output```

#### stochastic character mapping analysis
##### (produces *simmap_ctenosis.pdf* & *simmap_spongesis.pdf*)
##### (to run: adjust mydir variable)
```Rscript simmap_sym.R > simmap_sym.R.output &```

# Explanation of files

#### charmat_names.dat
list of taxa corresponding to charmat_vals.dat used by simmap_sym.R.

#### charmat_nomissindata.dat
character matrix with missing data removed. Used in likelihood_ratio_test.R. 

#### charmat_vals.dat
comma-separated state probabilities corresponding to taxa in charmat_names.dat.
column 1: separate sexes, column 2: hermaphroditic, column 3: asexual
missing data is coded as 0.333,0.333,0.333; all others have a single 1.000 entry

#### ctenosis.nex
composite tree in Nexus format used by likelihood_ratio_test.R and simmap_sym.R

#### spongesis.nex
composite tree in Nexus format with sponges as sister lineage.
used by likelihood_ratio_test.R and simmap_sym.R

