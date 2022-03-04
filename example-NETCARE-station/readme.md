# input files for setting up a FLEXPART-WRF for a NETCARE case

Also includes plotting scripts for vertical cross section of FLEXPART-WRF output and  plotting scripts for partial and total columnFLEXPART-WRF output

## files included

set_flex31_run_netcare.m  -> matlab code for writing the flexwrf.input file and AVAILABLE for FLEXPART-WRF v3.1 based on a given source location, output grid, and FLEXPART-WRF run options

AVAILABLE -> list of WRF ouptut files used, written by set_flex31_run_netcare.m based on the path to your WRF ouptut 								

flexwrf.input -> the generated flexwrf.input file based on set_fex31_run_netcare.m

plot_flexpart_netcare.m -> plotting script for netcdf output from FLEXPAERT-WRF

FLEXPART_COLUMN_wrf_2014071900_flex_20140720143015.png -> example total column plot

CROSS_SECTION_wrf_2014071900_flex_20140720143015.png -> example vertical cross section plot	

code written by louis.marelle@latmos.ipsl.fr and modified by jennie.thomas@univ-grneoble-alpes.fr

