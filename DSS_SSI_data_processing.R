####################################  DSS SSI DATA PROCESSING   #######################################
############################  Script downloads and imports Reanalyses data  #######################################
#This script processes Maine data for Decision Support System for SSI.                                                                     #
#AUTHOR: Benoit Parmentier                                                                      #
#DATE CREATED: 02/04/2014            
#DATE MODIFIED: 03/19/2014            
#Version: 1
#PROJECT: DSS-SSI project
#TO DO:
#1) deal with NA
#2) reclassifify NLCD classes
#3) resample at 100 m (both temp and NLCD data)
#4) summarize trends of temp by cities and counties
#
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

#...insert here

##############################
#### Parameters and constants  

script_path<-"/home/parmentier/Data/IPLANT_project/env_layers_scripts/" #path to script
#source(file.path(script_path,function_analyses_paper1)) #source all functions used in this script 1.

in_dir <- "/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/"
Maine_counties_file <- "county24.shp"
Maine_town_file <- "metwp24.shp"

#EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
#+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
CRS_reg <- "+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
file_format <- ".tif"
out_suffix <- "03192014"
w_extent_str <- c("-72 48 -65 41") #minx,maxy (upper left), maxx,miny (lower right)
use_reg_extent<- TRUE
setwd(in_dir)

out_dir <- in_dir

#out_path<-"/data/project/layers/commons/data_workflow/output_data"
out_dir <- file.path(out_dir,paste("output_data_",out_suffix,sep=""))

if (!file.exists(out_dir)){
  dir.create(out_dir)
  #} else{
  #  out_path <-paste(out_path..)
}else{
  file.remove(list.files(path=out_dir))
}

########## START SCRIPT #############

reg_area_poly <-readOGR(".",sub(".shp","",Maine_counties_file))                 #reading shapefile 
reg_town <-readOGR(".",sub(".shp","",Maine_town_file))                 #reading shapefile 


#<- "/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/ncar_ccsm3_0_sres_a1b_2050s_tmin_2_5min_no_tile_grd/tmin"
#r_fname <- "/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/ncar_ccsm3_0_sres_a1b_2050s_tmean_2_5min_no_tile_asc/tmean_1.asc"
#r <- raster(r_fname)
if (use_reg_extent==TRUE){
  #get extent in WGS84 (lat long)
  reg_area_poly_WGS84 <- spTransform(reg_area_poly,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
  reg_extent <- bbox(reg_area_poly_WGS84)
  w_extent <- c(reg_extent[1,1],reg_extent[2,2],reg_extent[1,2],reg_extent[2,1])
  w_extent <- paste(as.character(w_extent),collapse=" ")
}else{
  w_extent <- w_extent_str #this is in WGS84
}

## Temp processing: 2020s,2030s,2040s,2050s for tmin, tmax and tmean (e.g. 4*3 directories with 12 files...)

#1.download...not automated...
#2.unzip asc or grd (arc grid)
#3.convert to tif (gdal_tranlslate)
#4.reproject and clip/subset for Maine region (clip using count24?)

#loop through files...
l_f <- list.files(pattern="*asc.zip")
l_dir <-lapply(l_f,function(i){sub(".zip","",i)})
#unzip(l_f[[1]],exdir=l_dir[[1]])

for (i in 1:length(l_f)){
  unzip(l_f[[i]],exdir=l_dir[[i]])
}

# loop through files
for (i in 1:length(l_f)){
  f_list <- list.files(path=l_dir[[i]],pattern="*.asc",full.names=T)
  
  for (j in 1:length(f_list)){
    
    src_dataset <- file.path(in_dir,f_list[[j]])
    #in this case need to add folder name!!!
    dst_dataset <- paste(l_dir[[i]],sub(".asc",file_format,basename(src_dataset)),sep="_")
    #dst_dataset <- paste(sub(file_format,"_clipped",basename(src_dataset)),file_format,sep="")

    #dst_dataset <- file.path(out_dir,sub(".asc",file_format,basename(src_dataset)))
    dst_dataset <- file.path(out_dir,dst_dataset)
    #command_str <-paste("gdal_translate", src_dataset, dst_dataset ,sep=" ") #without spatial subsetting
    command_str <-paste("gdal_translate", 
                    "-projwin", w_extent, 
                    "-a_srs", paste("'","+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs","'",sep=" "), 
                    src_dataset, dst_dataset ,sep=" ")
    #Input file size is 161190, 104424        
    system(command_str)

    ## Now reproject 
    src_dataset <- dst_dataset
    #dst_dataset <- src_dataset
    dst_dataset <- paste(sub(file_format,"_projected",basename(src_dataset)),file_format,sep="")
    #dst_dataset <- file.path(out_dir,sub(".asc",file_format,basename(src_dataset)))
    dst_dataset <- file.path(out_dir,dst_dataset)

    command_str <- paste("gdalwarp", "-t_srs",paste("'",CRS_reg,"'",sep=""), src_dataset, dst_dataset,sep=" ")
    #t_rsrs : target/output spatial ref system
    system(command_str)

    #reg_extent_str <- paste(as.character(as.vector(reg_extent)),collapse=" ")
    #command_str <- paste("gdalwarp", "-t_srs",paste("'",CRS_reg,"'",sep=""), 
    #                     #"-te",paste("'",reg_extent_str,"'",sep=""),
    #                     "-te",reg_extent_str,
    #                     src_dataset, dst_dataset,sep=" ")
  }

}

#Quick chech in R raster
#ncar_pred_files <- list.files(path=out_dir,pattern="*.projected",full.names=T)
#r_ncar_stack <- stack(ncar_pred_files)

#levelplot(r_ncar_stack,layer=1:3)
#plot(r_ncar_stack,y=1,col=matlab.like(25),colNA=c("white")) #deal with NA val!!!
## NLCD processing: 1992, 2001, 2006

#1.download...not automated?
#2.unzip tif (arc grid)
#3.convert to tif (gdal_tranlslate)
#4.reproject and clip/subset for Maine region (clip using count24?)

### Project the WGS84 extent to the target SRS
#Make this a function...!!!
#w_extent_num <- as.numeric(unlist(strsplit(w_extent," ")))
#mat_coord_extent<- matrix(c(-72,-65,48,41),2)
#mat_coord_extent<- t(matrix(w_extent_num,nrow=2,ncol=2,byrow=FALSE))
#dat<- as.data.frame(mat_coord_extent)
#coordinates(dat) <- mat_coord_extent
#proj4string(dat)<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#dat_proj <- spTransform(dat,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) 
dat_proj <- spTransform(reg_area_poly,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) 

#reg_area_poly
reg_extent<- bbox(dat_proj)
#reg_extent <- bbox(reg_area_poly_WGS84)
w_extent <- c(reg_extent[1,1],reg_extent[2,2],reg_extent[1,2],reg_extent[2,1])
w_extent <- paste(as.character(w_extent),collapse=" ")

#w_extent<-project(mat_coord_extent,
#    proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
#w_extent <- paste(as.character(w_extent),collapse=" ")
#w_extent <- paste(as.character(bbox(dat)),collapse=" ")

#as.character(as.vector(t(coordinates(dat))))
#w_extent <- paste(as.character(t(coordinates(dat))),collapse=" ")

### end of function

f_list <- list.files(path=in_dir,pattern="nlcd.*.img",full.names=T)

for(j in 1:length(f_list)){
    # ...make this part of function to match with earlier...
    src_dataset <- f_list[[j]]
    #in this case need to add folder name!!!
    out_prefix <- ""
    dst_dataset <- paste(out_prefix,sub(extension(src_dataset),file_format,basename(src_dataset)),sep="")
    #dst_dataset <- paste(sub(file_format,"_clipped",basename(src_dataset)),file_format,sep="")

    #dst_dataset <- file.path(out_dir,sub(".asc",file_format,basename(src_dataset)))
    dst_dataset <- file.path(out_dir,dst_dataset)
    #command_str <-paste("gdal_translate", src_dataset, dst_dataset ,sep=" ") #without spatial subsetting
    #clipped using extent...
    #w_extent <- as.numeric(unlist(strsplit(w_extent_str," ")))
    
    command_str <-paste("gdal_translate", 
                    "-projwin", w_extent, 
                    "-a_srs", paste("'","+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs","'",sep=" "), 
                    src_dataset, dst_dataset ,sep=" ")
    #Input file size is 161190, 104424        
    system(command_str)
    
    ## Now reproject 
    src_dataset <- dst_dataset
    #dst_dataset <- src_dataset
    dst_dataset <- paste(sub(file_format,"_projected",basename(src_dataset)),file_format,sep="")
    #dst_dataset <- file.path(out_dir,sub(".asc",file_format,basename(src_dataset)))
    dst_dataset <- file.path(out_dir,dst_dataset)

    command_str <- paste("gdalwarp", "-t_srs",paste("'",CRS_reg,"'",sep=""), src_dataset, dst_dataset,sep=" ")
    #t_rsrs : target/output spatial ref system
    system(command_str)


}

r_stack <- stack(list.files(pattern="^nlcd.*projected.tif$",path=out_dir,full.names=T))

### Now reclass the image

### Now adjust the resolution to be 100m

### Now summarize by table of towns/counties

list.files(pattern="tmean.*.projected.tif$")
list.files(path=out_dir,pattern="tmax.*.projected.tif$")
r_tmin <- stack(mixedsort(list.files(path=out_dir,pattern="tmin.*.projected.tif$",
                full.names=TRUE)))

r_reg_area <- rasterize(reg_area_poly,r_tmin,field="COUNTY") #this is the prediction from lm model
projection(r_tmin) <- projection(r_reg_area)
tb <- as.data.frame(zonal(r_tmin,r_reg_area))

col_pal <- c("blue","green","red","black")
y_range <- range(tb[1,])
plot(1:12,tb[1,2:13],type="b",col=col_pal[1],ylim=y_range) #2020s
lines(1:12,tb[1,14:25],type="b",col=col_pal[2]) #2030s
lines(1:12,tb[1,26:37],type="b",col=col_pal[3]) #2040s
lines(1:12,tb[1,38:49],type="b",col=col_pal[4]) #2050s


plot(1:4,tb[1,c(8,20,32,44)],type="b",col=col_pal[1]) #2020s

names(tb)<-c("zone",
  paste(1:12,c("y20"),sep=""),
  paste(1:12,c("y30"),sep=""),
  paste(1:12,c("y40"),sep=""),
  paste(1:12,c("y50"),sep=""))

tb$county <- as.character(unique((reg_area_poly$COUNTY)))
#change names and write out as point data.frame...
  
#gCentroid(reg_area_poly,byid=T,id="COUNTY")
write.table(tb,file="tmin_summary_by_county.txt",sep=",")

##Now try with extract

###### End of script ###########

# command_str <-paste("gdal_translate", "nlcd2006_landcover_4-20-11_se5.img","nlcd2006_landcover_4-20-11_se5.tif",sep=" ")
# #Input file size is 161190, 104424
# system(command_str)
# 
# gdalwarp -t_srs '+proj=utm +zone=11 +datum=WGS84' raw_spot.tif utm11.tif
# #t_rsrs : target/output spatial ref system
# in_r <- "" #input raster
# out_r <- "" #output raster
# paste("gdalwarp", "t_srs", paste(',CRS_reg,',sep=""), in_r, out_r, sep=" " )
# 
# command_str <- paste("gdalwarp", "-dstnodata 0", 
# "-cutline M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/county24.shp", 
# "-crop_to_cutline", 
# "-of GTiff", 
# "M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/nlcd2001_landcover_v2_2-13-11.img", 
# "M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/output_data_02042014/test55.tif",sep=" ")
# 
# dst_dataset <- file.path(out_dir,"test55.tif")
# command_str <- paste("gdalwarp", "-dstnodata", paste("'","0","'",sep=""), 
# "-cutline", paste("'","home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data//county24.shp","'",sep=""), 
# "-crop_to_cutline", 
# #"-of GTiff", 
# f_list[2], 
# dst_dataset,sep=" ")
# 
# command_str <- paste("gdalwarp", #"-dstnodata", paste("'","0","'",sep=""), 
# "-cutline home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data//county24.shp", 
# #"-of GTiff", 
# f_list[1], 
# dst_dataset,sep=" ")
# 
# system(command_str)
# 
# 
# gdalwarp -dstnodata 0 
# -q -cutline M:\Data\IPLANT_project\Maine_interpolation\DSS_SSI_data\county24.shp 
# -crop_to_cutline 
# -of GTiff 
# M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/nlcd2001_landcover_v2_2-13-11.img 
# M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/output_data_02042014/test02182014.tif
# 
# gdalwarp -s_srs EPSG:4326 
# -q -cutline M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/county24.shp 
# -dstalpha 
# -of GTiff 
# M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/ncar_ccsm3_0_sres_a1b_2020s_tmax_2_5min_no_tile_asc/tmax_4.asc 
# M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/output_data_02042014/test50.tif
# 
# gdalwarp -s_srs EPSG:4326 -q -cutline M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/county24.shp -of GTiff M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/ncar_ccsm3_0_sres_a1b_2020s_tmax_2_5min_no_tile_asc/tmax_4.asc M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/output_data_02042014/test52.tif
# gdalwarp -s_srs EPSG:4326 
# -t_srs "+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0" 
# -q -cutline M:\Data\IPLANT_project\Maine_interpolation\DSS_SSI_data\county24.shp 
# -dstalpha 
# -of GTiff M:\Data\IPLANT_project\Maine_interpolation\DSS_SSI_data\ncar_ccsm3_0_sres_a1b_2020s_tmax_2_5min_no_tile_asc\tmax_1.asc 
# M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/output_data_02042014/test45.tif
# 
# gdalwarp -s_srs EPSG:4326 
# -t_srs EPSG:2960 
# -q -cutline M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/county24.shp 
# -dstalpha -of GTiff 
# M:\Data\IPLANT_project\Maine_interpolation\DSS_SSI_data\ncar_ccsm3_0_sres_a1b_2020s_tmean_2_5min_no_tile_asc\tmean_3.asc 
# M:/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/output_data_02042014/test48.tif


#system("gdalwarp -t_srs EPSG:4326 /home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data//nlcd2006_landcover_4-20-11_se5.img test.55.tif")