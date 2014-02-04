####################################  DSS SSI DATA PROCESSING   #######################################
############################  Script downloads and import Reanalyses data  #######################################
#This script processes Maine data for DSS SSI to examined spatial patterns and compare them to daily devation 
#from CAI and FSS methods
#Figures and data for the contribution of covariate paper are also produced.                                                                     #
#AUTHOR: Benoit Parmentier                                                                      #
#DATE CREATED: 02/04/2014            
#DATE MODIFIED: 02/04/2014            
#Version: 1
#PROJECT: DSS-SSI project                                    #
#################################################################################################

###Loading R library and packages                                                      
library(gtools)                              # loading some useful tools 
library(mgcv)                                # GAM package by Simon Wood
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gstat)                               # Kriging and co-kriging by Pebesma et al.
library(fields)                              # NCAR Spatial Interpolation methods such as kriging, splines
library(raster)                              # Hijmans et al. package for raster processing
library(gdata)                               # various tools with xls reading
library(rasterVis)
library(parallel)
library(maptools)
library(maps)
library(reshape)
library(plotrix)
library(plyr)

#### FUNCTION USED IN SCRIPT

##############################
#### Parameters and constants  

script_path<-"/home/parmentier/Data/IPLANT_project/env_layers_scripts/" #path to script
#source(file.path(script_path,function_analyses_paper1)) #source all functions used in this script 1.

in_dir <- "/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/"
Maine_counties_file <- "county24.shp"
Maine_town_file <- "metwp24.shp"

CRS_reg <- "+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
file_format <- ".tif"
out_suffix <- "02042014"
setwd(in_dir)

out_dir <- in_dir

#out_path<-"/data/project/layers/commons/data_workflow/output_data"
out_dir <- file.path(out_dir,paste("output_data_",out_suffix,sep=""))

if (!file.exists(out_dir)){
  dir.create(out_dir)
  #} else{
  #  out_path <-paste(out_path..)
}

########## START SCRIPT #############

reg_counties <-readOGR(".",sub(".shp","",Maine_counties_file))                 #reading shapefile 
reg_town <-readOGR(".",sub(".shp","",Maine_town_file))                 #reading shapefile 

reg_extent <- bbox(reg_counties)
#<- "/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/ncar_ccsm3_0_sres_a1b_2050s_tmin_2_5min_no_tile_grd/tmin"
#r_fname <- "/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/ncar_ccsm3_0_sres_a1b_2050s_tmean_2_5min_no_tile_asc/tmean_1.asc"
#r <- raster(r_fname)


## Temp processing: 2020s,2030s,2040s,2050s for tmin, tmax and tmean (e.g. 4*3 directories with 12 files...)

#1.download...not automated...
#2.unzip asc or grd (arc grid)
#3.convert to tif (gdal_tranlslate)
#4.reproject and clip/subset for Maine region (clip using count24?)

#loop through files...
l_f <- list.files(pattern="*asc.zip")
l_dir <-lapply(l_f,function(i){sub(".zip","",i)})
unzip(l_f[[1]],exdir=l_dir[[1]])

# loop through files
f_list <- list.files(path=l_dir[[1]],pattern="*.asc",full.names=T)
src_dataset <- file.path(in_dir,f_list[[1]])
#in this case need to add folder name!!!
dst_dataset <- paste(l_dir[[1]],sub(".asc",file_format,basename(src_dataset)),sep="_")
#dst_dataset <- file.path(out_dir,sub(".asc",file_format,basename(src_dataset)))
dst_dataset <- file.path(out_dir,dst_dataset)
command_str <-paste("gdal_translate", src_dataset, dst_dataset ,sep=" ")
#Input file size is 161190, 104424
system(command_str)

src_dataset <- dst_dataset
#dst_dataset <- src_dataset
dst_dataset <- paste(sub(file_format,"_projected",basename(src_dataset)),file_format,sep="")
#dst_dataset <- file.path(out_dir,sub(".asc",file_format,basename(src_dataset)))
dst_dataset <- file.path(out_dir,dst_dataset)

command_str <- paste("gdalwarp", "-t_srs",paste("'",CRS_reg,"'",sep=""), src_dataset, dst_dataset,sep=" ")
#t_rsrs : target/output spatial ref system
system(command_str)


src_dataset <- dst_dataset
dst_dataset <- paste(sub(file_format,"_cliped",basename(src_dataset)),file_format,sep="")

reg_extent_str <- paste(as.character(as.vector(reg_extent)),collapse=" ")
command_str <-paste("gdal_translate", "-projwin", reg_extent_str, src_dataset, dst_dataset ,sep=" ")
#Input file size is 161190, 104424
system(command_str)

reg_extent_str <- paste(as.character(as.vector(reg_extent)),collapse=" ")
command_str <- paste("gdalwarp", "-t_srs",paste("'",CRS_reg,"'",sep=""), 
                     #"-te",paste("'",reg_extent_str,"'",sep=""),
                     "-te",reg_extent_str,
                     src_dataset, dst_dataset,sep=" ")

## NLCD processing: 1992, 2001, 2006

#1.download...not automated?
#2.unzip tif (arc grid)
#3.convert to tif (gdal_tranlslate)
#4.reproject and clip/subset for Maine region (clip using count24?)


r_nlcd2001 <- raster(f_list[[1]])
command_str <-paste("gdal_translate", "nlcd2006_landcover_4-20-11_se5.img","nlcd2006_landcover_4-20-11_se5.tif",sep=" ")
#Input file size is 161190, 104424
system(command_str)

gdalwarp -t_srs '+proj=utm +zone=11 +datum=WGS84' raw_spot.tif utm11.tif
#t_rsrs : target/output spatial ref system
in_r <- "" #input raster
out_r <- "" #output raster
paste("gdalwarp", "t_srs", paste(',CRS_reg,',sep=""), in_r, out_r, sep=" " )





