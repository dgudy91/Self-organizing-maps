## R code to run Self Organising Maps on geopotential height daily anomaly
## Input data: ERA - Interim 500hpa data - daily anomalies
## one netcdf file over analysis period

#working example using google colab available here: https://github.com/dgudy91/2021_AAPP_ML_workshop/blob/main/AAPP_ML_workshop_SOM_tutorial_data_subset.ipynb

#### Load library packages ####
library("tidyverse")
install.packages("tidync")
library("tidync") #A Tidy Approach to 'NetCDF' Data Exploration and Extraction
install.packages("kohonen")
library("kohonen") #SOM package
install.packages("viridis")
library("viridis") #Colorblind-Friendly Color Maps for R
install.packages("ggcorrplot") 
library("ggcorrplot") #Visualization of a correlation matrix using ggplot2
library("ggplot2") #plot results
install.packages("yardstick")
library("yardstick") #package to estimate how well models are working using tidy data principles


# folder and file locations

data_dir <- "File_path/Data" #update to your filepath to your data folder
setwd(data_dir) #set working directory to folder containing the netcdf file you want to calculate SOM from
f_name  <- "z_500_dailyanomaly_ERAI_1979-2018_AAT.nc" #update to your file name
#use cdo mergetime function to create nc file of analysis period

#### Load netcdf file ####

#Make netcdf connection, filter by time, read into tibble, create date column

z500_daily <- tidync::tidync(f_name) %>%
  hyper_filter(time = lubridate::as_datetime(time * 3600, origin = "1900-01-01 00:00:00") < 
                 lubridate::as_datetime("1983-12-31 23:00:00")) %>%
  hyper_tibble() %>%
  mutate(date = lubridate::as_datetime(time * 3600, origin = "1900-01-01 00:00:00")) %>%
  dplyr::select(-c(lev, time)) #get rid of unwanted columns

#Set up output directory
rundate = 'YYYY-MM-DD' #update for each run
run_num = 'testrun_01' #update for each run


output <-"File_path/Results" #update to your filepath to results folder
dir.create(file.path(output,run_num), showWarnings = FALSE) #create new folder based on run number
out_dir <- file.path(output,run_num) 

data_var = 'daily_z500_anomaly'

setwd(out_dir) #set working directory to results folder (sub folder by run number)


#check data
head(z500_daily)
summary(z500_daily)

#convert geopotential pressure anomaly to geopotential height anomaly by dividing by gravity
#overwrite the z column with the calculated value
g = 9.80665
z500_daily$z = z500_daily$z/g
head(z500_daily)
summary(z500_daily)

#### Transform dataset for SOM algorithm ####

#Step 1 - convert the dataset from long format to wide format
#Step 2 - Remove the date column

#Currently the dataset consists of 4 columns: z, lon, lat and time - this is called long form
#The SOM algorithm input requires the dataset to have time as rows, and each combination of lat/lon as columns - this is called wide form

#Restructure the dataset using the pivot function to transform the long form dataset into wide from. 

#Pivots data from long to wide format. 
#Each row is a date, each column is the geopotential at an individual coordinate
z500_daily_wide <- z500_daily %>%
  pivot_wider(., id_cols = date, names_from = c(lat, lon), values_from = `z`)

#check transform is correct
head(z500_daily_wide)

#Remove the date column for the SOM calculation - make the dataset one dimensional
z500_daily_1d <- z500_daily_wide %>%
  dplyr::select(-date)

#### SET UP SOM ALGORITHM PARAMETERS ####

#parameters used in https://journals.ametsoc.org/view/journals/clim/34/3/JCLI-D-20-0297.1.xml

#SOM sequential mode: training is goverened by: 
# 1. learning rate parameter (alpha)
# 2. neighbourhood radius (radius)
# 3. number of iterations (rlen)

#default parameters: for training
#rlen = 100
#alp = c(0.05,0.01)
#rad = c(5,1)

# number of iterations - rlen
rlen = 1000             # number of times the complete dataset will be presented to the network
#rlen -> use 100 for initial testing -> increase to 1000+ when actually running

# learning rate parameter - alpha
alp = c(0.05,0.01)     # learning rate -> magnitude each node pattern is updated -> not used in batch learning mode
# Requires vector of two numbers indicating the amount of change
# Default (decrease from 0.05 to 0.01 linearly over rlen updates)

# neighbourhood radius - radius
rad = c(4,0)             # Number of surrounding nodes activated - # 75% SOM, 25% clustering.
# value < 1 no surrounding nodes activated = k-means clustering
#100% clustering example: rad = c(4,1)

#grid structure of SOM ("rect", "hexagonal)
gr = "rect"  #update manually in SOM algorithm.  This is for filename setup. 

#data type - "raw","detrended","daily_anom"
dt = "daily_anomaly"

## SOM parameters - more information in documentation - kohonen_R_2018.pdf ##
#https://cran.r-project.org/web/packages/kohonen/index.html

# parameter sensitivity to test -> from Gibson et al 2017 - https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JD026256 
# starting point of alpha (learning rate) and radius generally have limited influence on patterns produced and their associated realism  (pattern correlation to actual daily field)
# common practice - reduce radius parameter linearly to 1, then stop training at this point -> ensures winning and surrounding node updated and topological ordering preserved -> SOMs for projections

# optimise realism of SOM patterns -> SOMs for clustering -> reduce generalisation of SOM patterns
# Gibson et al 2017 test how end point of radius parameter influences SOM realisation

#SOM p (projection) parameters
# radius start point = 5, radius end 

#### SET UP FILENAMES FOR OUTPUT FILES/PLOTS ####

filename_end = paste("n",Nnodes,"_",dt,"_rlen",rlen,"_alp",alp[1],"to",alp[2],"_rad",rad[1],"to",rad[2],".png", sep="")  #sep removes spaces in filename
print(filename_end) #check filename

#write txt file with parameters for run - save in run_num subfolder
sink(paste0(run_num,"_Parameters.txt"))
cat(paste0("Run Date ",rundate))
cat("\n")
cat(paste0("Parameters: data = ",dt,",dist.fcts = euclidean, grid = ",gr,"Nodes = ",Nnodes, "rlen = ",rlen,", alp",alp[1],"to ",alp[2],", rad",rad[1],"to ",rad[2]))
sink()

#### RUN THE SOM ####

set.seed(5) #this is important for consistent results between runs

#input data 
data = z500_daily_1d

#Play around with number of nodes
Nnodes = 9
nx = 3
ny = 9

#SOM algorithm
print("running SOM algorithm.......")
K_SOM = som(X = data.matrix(data),
            grid = somgrid(nx,ny,"rect","gaussian"),
            rlen = rlen,
            alpha = alp,
            radius = rad,
            keep.data = T,
            dist.fcts = c("euclidean"),
            normalizeDataLayers = FALSE)

print("finished running SOM")


#Save the variables
K_SOM_d    <- data.frame(K_SOM$distances)        # get errors -> distances between grids
K_SOM_SOM  <- data.frame(K_SOM$codes)            # get clusters grid
tmp        <- data.frame(t(K_SOM_SOM))           # transpose
K_SOM_SOMc <- data.frame(unlist(tmp))           # concatenate codebook vectors to 1 row
K_SOM_win  <- data.frame(K_SOM$unit.classif)     # get clusters win
K_SOM_grid <- data.frame(K_SOM$grid[["pts"]])    # get grid structure

#column binds two vectors - one of the individual dates and one of the winning nodes
winning_nodes <- cbind(z500_daily_subset_wide$date,K_SOM_win)      #finds the winning nodes of each timestep - daily data      
colnames(winning_nodes) <- c("date","node")

#Attach the winning node data to original dataset - match by date column
z500_daily_winning_nodes <- dplyr::left_join(z500_daily, winning_nodes, by = "date")

#SAVE TXT FILE OF WINNING NODES 
write.table(winning_nodes, file=paste("SOM_tutorial_test_",data_var,"_",nx,"_",ny,".txt", sep=""),row.names=F, col.names=T, quote=F)
print("file written as ....", quote=F)
print(paste("SOM_tutorial_test_",data_var,"_",nx,"_",ny,".txt", sep=""))

#### Plot SOM composite maps ####

#Calculate the SOM composite for each node, for each coordinate
som_means <- z500_daily_subset_winning_nodes %>%
  group_by(lat, lon, node) %>%
  summarise(som_z500_mean = mean(`z`, na.rm = TRUE))


#Plot the SOM composites with coastlines.
install.packages("maptools")
library(maptools)
data(wrld_simpl)


ggplot() + 
  geom_tile(data = som_means, aes(x = lon, y = lat, fill = som_z500_mean)) + 
  scale_fill_viridis() +
  facet_wrap(~node) +
  theme_bw() +
  geom_polygon(data=wrld_simpl, 
               aes(x=long, y=lat, group=group), fill='grey', colour = 'black', alpha = 0.5)  +
  coord_cartesian(xlim=c(50,160), ylim=c(-70,-42)) +
  geom_contour(data = som_means, aes(x = lon, y = lat, z = som_z500_mean), colour = "grey30") +
  labs(title = "Composite z500 of each SOM - SOM input: daily anomalies 1979-1983",
       fill = "Geopotential height anomaly (m)")


#### EVALUATE SOM PERFORMANCE ####
#Evaluate the SOM performance for capturing 'real world' weather patterns
#Method: Pearson pattern correlation between the composite winning node pattern and each daily input assigned to the winning node.
#Refer to Gibson et al. 2017 ( https://doi.org/10.1002/2016JD026256); Udy et al. 2021 (https://doi.org/10.1175/JCLI-D-20-0297.1) for more details.


nodes <- seq(from = 1, to = Nnodes, by = 1)

#write a function to calculate correlation score for each day

calculate_som_temporal_correlation <- function(nodes){
  
  som_temporal_correlation <- left_join(z500_daily_winning_nodes %>%
                                          filter(node == nodes), 
                                        som_means %>% filter(node == nodes)) %>%
    group_by(date) %>%
    yardstick::rsq(., `z`, som_z500_mean) %>%
    mutate(node = nodes)
  
  som_temporal_correlation
  
}


# PLOT THE CORRELATION SCORES FOR EACH SOM NODE
#as the number of nodes increase, the correlation scores should improve
#as the number of nodes decrease, the correlation scores will likely go down. 
#Determining the right number of nodes is subjective, but you can be guided by how many patterns you roughly
#expect to find based on weather maps, and how sensitive the results are. 

som_temporal_correlations <- purrr::map_dfr(nodes, calculate_som_temporal_correlation)

ggplot(data = som_temporal_correlations, aes((.estimate ^ 0.5), colour = node, group = node)) +
  geom_density(size = 2) +
  scale_color_viridis() +
  #facet_wrap(~node) + #this line to plot them as separate facets
  theme_bw() +
  labs(x = "Pearson correlation")


#### OTHER PLOTS AVAILABLE IN KOHONEN PACKAGE ####

#### quality & count ####
#quality - shows the mean distance mapped to a unit to the codebook vector of that unit. 
# smaller distances = objects (e.g daily z500) represented better by the codebook vectors. 
# counts - number of objects mapped to SOM nodes - frequency?

counts <- plot(K_SOM, type = "counts", 
               palette = plasma, ncolors = 12, 
               shape = "straight")
quality <- plot(K_SOM, type = "quality", 
                palette = viridis, ncolors = 12,
                shape = "straight")

#U matrix/distance neighbours & hierachial clustering
plot(K_SOM,type = "dist.neighbours", 
     main = "SOM neighbour distances & clusters",
     palette = viridis,
     shape = "straight")
#hierarchial clustering to cluster the codebook vectors
som.hc <- cutree(hclust(object.distances(K_SOM,"codes")),6)  #random no. of clusters selected
add.cluster.boundaries(K_SOM, som.hc)

#Check the training progress (unlikely to reach optimal level in training)
plot(K_SOM, type="changes")


