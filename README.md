# Self-organizing-maps
R code to run Kohonen Self Organizing Maps on geopotential height data. Code steps through how to pre-process data, run the self-organising map algorithm and evaluate the output using atmospheric data from reanalysis climate models. 

Adapted from code developed for https://journals.ametsoc.org/view/journals/clim/34/3/JCLI-D-20-0297.1.xml

Thanks to Peter Gibson for providing the initial SOM code in R to get me started. You can find his paper on using SOMs for studying climate extremes here: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JD026256
Also thanks to Marc Mallet for helping me improve the structure and flow of the code.

_What is a Self Organising Map?_
A self-organising map (SOM) is a type of artificial neural network that is trained using unsupervised learning to produce a reduced dimensional, discretised representation of the input dataset. In climate science, SOMs are often used to identify key regional weather patterns (e.g. High pressure/Low pressure systems) in an area of interest (without having to manually look at each daily weather map)

_Why use a Self Organising Map?_
The key purpose for using a SOM is to simplify the data by finding patterns. When grouping/classifying weather types a key assumption is that daily weather patterns can be split into a specified number of "types". For this reason the SOM algorithm is often preferred for classifying weather over traditional discrete clustering methods (e.g. k-means). SOMs are able to account for continuity and nonlinearity (Jiang et al. 2012) more than clustering which provides a more realistic representation of the continuous movement of weather patterns.

_What package to use?_
The 'Kohonen' R package is very popular in climate science. The documentation is available here: https://cran.r-project.org/web/packages/kohonen/kohonen.pdf

SOM parameters to consider: 
1.	Number of nodes and structure of grid (e.g. 3 x 3 grid = 9 output nodes)
The output results are most sensitive to this parameter. As you reduce the number of output nodes the results are more generalised & become closer to the mean of the field. As you increase the number of output nodes the results become more specific to the input data. Finding the balance of the number of nodes is subjective and depends on what type of features/questions you are trying to answer. 
2.	Number of iterations (rlen)
number of times the complete dataset will be presented to the network
3.	Neighbourhood Radius (rad)
Number of surrounding nodes from the winning node activated each time the winning node is selected. This parameter is a key difference between k-means and SOM clustering. This parameter is generally set to start high (e.g. 4) and reduce to either 1 (100% SOM for all iterations) or 0 (75% of iterations SOM & 25% discrete clustering). In the R Kohonen package, a radius value below 1 = no surrounding nodes activated = k-means clustering. 
4.	Learning rate (alp)
Magnitude that each winning node is influenced by the input data. For the first few iterations this is set to be a fast learning rate and reduces to a slower learning rate near the end of the interations to improve accuracy

DATA PRE-PROCESSING

Prior to loading the netcdf data into R the following steps are generally required (depending on your research questions):
1.	Download climate reanalysis data/variables of interest - this code uses ERA-Interim 500hpa geopotential height. 
2.	Calculate daily anomalies to remove seasonality influence on the SOM output. The seasonal signal will dominate the variability detected by the SOM algorithm - so unless you are interested in seasonality remove prior to running the SOM. Climate Data Operators tools are faster than using R or python to do this (https://code.mpimet.mpg.de/projects/cdo/)

Some useful cdo functions for processing netcdf data prior to loading into R
•	to merge multiple netcdf files together over time:
cdo -b F64 mergetime infiles*.nc mergetime_file.nc
•	to select a region (e.g Southern Hemisphere):
cdo -sellonlatbox,0,360,-90,0 mergetime_file.nc mergetime_regional_file.nc
•	to calculate daily anomalies:
•	a) convert 6hr dataset to daily using cdo daymean
•	b) create a file of daily climatology using cdo ydaymean
•	c) calculate the daily anomaly using cdo ydaysub infile1.nc climatologyfile.nc outfile.nc
e.g. cdo -r -f nc4 -b 64 ydaysub z_500_dailymean_ERAI_historical_an-pl_1979-2018_SH.nc z_500_dailymeanstat_ERAI_historical_an-pl_1979-2018_SH.nc daily_anomaly_file.nc

