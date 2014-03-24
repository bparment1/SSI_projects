#!/usr/bin/python
#
#Script to process data for the DSS-SSI website.
#Data are first clipped to match the Maine study area.
#
# TODO:
#  - functionalize to encapsulate high level procedural steps
#
# Authors: Benoit Parmentier 
# Created on: 03/24/2014
# Updated on: 03/24/2014

import os, glob
import subprocess
import re, zipfile
import datetime, calendar
import ftplib
#import grass.script as gs
import argparse
import shutil
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from osgeo import gdal_array
from osgeo import gdalconst

#------------------
# Functions used in the script 
#------------------


def create_dir_and_check_existence(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

#######################################################################
######################### BEGIN SCRIPT  ###############################
#--------------------------------------------
# Script run by arguments from the shell?
#--------------------------------------------
            
def main():
    #--------------------------------------------
    # Download and Calculate monthly climatology for daily LST time series for a specific area
    #--------------------------------------------

    # TODO: set up a (temporary?) GRASS database to use for processing? code
    # currently assumes it's being run within an existing GRASS session
    # using an appropriately configured LOCATION...
    #
 
    ########## START OF SCRIPT ##############
    #### Modified by Benoit on May 13, 2013 
    
    ### INPUT Parameters
    #Inputs from R?? there are 9 parameters
   #tiles = ['h11v08','h11v07','h12v07','h12v08','h10v07','h10v08'] #These tiles correspond to Venezuela.
    #tiles= ['h08v04','h09v04']    
    #start_year = 2001
    #end_year = 2010
    #start_month=1
    #end_month=12
    #hdfdir =  '/home/layers/commons/modis/MOD11A1_tiles' #destination file where hdf files are stored locally after download.
    #hdfdir =  '/home/parmentier/Data/IPLANT_project/MOD11A1_tiles' #destination file where hdf files are stored locally after download.
    #night=1    # if 1 then produce night climatology
    #out_suffix="_03192013"
    #download=0  # if 1 then download files
   
    #Passing arguments from the shell...using positional assignment
    parser = argparse.ArgumentParser()
    parser.add_argument("tiles", type=str, help="list of Modis tiles")
    parser.add_argument("start_year", type=int, help="start year")
    parser.add_argument("end_year", type=int, help="end year")
    parser.add_argument("start_month", type=int, help="start month")
    parser.add_argument("end_month", type=int, help="end month")
    parser.add_argument("hdfdir", type=str, help="destination/source directory for hdf file")
    parser.add_argument("night", type=int, help="night")
    parser.add_argument("download", type=int, help="out_suffix")
    parser.add_argument("out_suffix", type=str, help="out_suffix")

    myargs = parser.parse_args()
    
    #### parse input parameters
    tiles = myargs.tiles #These tiles correspond to Venezuela.
    start_year = myargs.start_year
    end_year = myargs.end_year 
    end_month = myargs.end_month #not used...to be removed
    start_month= myargs.start_month #not used...to be removed
    hdfdir =  myargs.hdfdir
    night= myargs.night    # if 1 then produce night climatology
    out_suffix= myargs.out_suffix #"_03192013"
    download=myargs.download# if 1 then download files
    
    tiles =tiles.split(",") #Create a list from string
    #need to add on the fly creation of folder for each tile!!
    
    in_dir = "/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/"
    Maine_counties_file = "county24.shp"
    Maine_town_file = "metwp24.shp"

    #EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
    #+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
    #CRS_reg = "+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    #http://spatialreference.org/ref/epsg/2037    
    CRS_reg = "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    file_format = ".tif"
    out_suffix = "03242014"
    w_extent_str = "-72 48 -65 41" #minx,maxy (upper left), maxx,miny (lower right)
    use_reg_extent = True
    os.chdir(in_dir)

    out_dir = in_dir

    #out_path<-"/data/project/layers/commons/data_workflow/output_data"
    out_dir = "output_data_"+out_suffix
    out_dir = os.path.join(in_dir,out_dir)
    create_dir_and_check_existence(out_dir)
            
    ########## START SCRIPT #############
    CRS_WGS84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    #Read in layers from data source,there is only one layer
    reg_area_poly = ogr.Open(Maine_counties_file).GetLayer()
    reg_town = ogr.Open(Maine_town_file).GetLayer()
    #reg_area_poly <-readOGR(".",sub(".shp","",Maine_counties_file))                 #reading shapefile 
    #reg_town <-readOGR(".",sub(".shp","",Maine_town_file))                 #reading shapefile 

    if use_reg_extent==True:
        #get extent in WGS84 (lat long)
        src_dataset = os.path.join(in_dir,Maine_counties_file)
        dst_dataset = os.path.splitext(src_dataset)[0]+"_WGS84"+os.path.splitext(src_dataset)[1]
        dst_dataset = os.path.join(out_dir,dst_dataset)

        #Write a quick wrapper function later on...
        cmd_str = "".join(["ogr2ogr",
                          " "+"-f "+"'ESRI Shapefile'",  
                          " "+"-t_srs"+" '"+CRS_WGS84+"'",
                          " "+"-overwrite",
                          " "+dst_dataset, 
                          " "+src_dataset])
        os.system(cmd_str)
        
        #reg_area_poly_WGS84 <- spTransform(reg_area_poly,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
        #reg_lyr= ogr.Open(dst_dataset).GetLayer(0)    
        #reg_extent = reg_lyr.GetExtent() #tuple object, not working right now
        #ogrinfo worldborder_sinusoidal.shp -so -al | grep Extent #givs xmin,ymin and xmax ymax
        cmd_str = "".join(["ogrinfo ",dst_dataset," -so -al"," | grep Extent"])
        os.system(cmd_str)
        
        reg_extent = subprocess.check_output(cmd_str, shell=True)
        reg_extent=reg_extent.split(" - ")
        coordpt1 = re.compile("\((.*)\)").search(reg_extent[0]).group(1)
        coordpt2 = re.compile("\((.*)\)").search(reg_extent[1]).group(1)
        w_extent = "".join([coordpt1.split(",")[0],
                           " "+coordpt2.split(",")[1],
                           " "+coordpt2.split(",")[0],
                           " "+coordpt1.split(",")[1]])


        #[1] "-71.083923738042 47.4598539782516 -66.8854440488051 42.9171281482886"
        
        reg_extent <- bbox(reg_area_poly_WGS84)
        w_extent <- c(reg_extent[1,1],reg_extent[2,2],reg_extent[1,2],reg_extent[2,1])
        w_extent <- paste(as.character(w_extent),collapse=" ")
    }else{
      w_extent <- w_extent_str #this is in WGS84
    }

    
    #in this case need to add folder name!!!
 
    
    ## Temp processing: 2020s,2030s,2040s,2050s for tmin, tmax and tmean (e.g. 4*3 directories with 12 files...)

    #1.download...not automated...
    #2.unzip asc or grd (arc grid)
    #3.convert to tif (gdal_tranlslate)
    #4.reproject and clip/subset for Maine region (clip using count24?)

    #loop through files...
    fileglob = "*asc.zip"
    pathglob = os.path.join(in_dir, fileglob)
    l_f = glob.glob(pathglob)
    l_dir = map(lambda x: os.path.splitext(x)[0],l_f) #similar to laplly in R

    #unzip(l_f[[1]],exdir=l_dir[[1]])
    for i in range(0,len(l_f)):
        with zipfile.ZipFile(l_f[i], "r") as z:
            z.extractall(l_dir[i])
            
 # loop through files
for i in 1:len(l_f)):
    fileglob="*.asc"
    pathglob = os.path.join(l_dir[i], fileglob)
    f_list = glob.glob(pathglob)
    
    for j in range(0,len(f_list)):
        src_dataset = os.path.join(in_dir,f_list[j])
        dst_dataset = os.path.splitext(os.path.basename(src_dataset))[0]+file_format
        dst_dataset = os.path.join(out_dir,dst_dataset)
        
        cmd_str = "".join(["gdal_translate",
                           " "+"-projwin"+" "+w_extent,
                           "-a_srs"+" '"+CRS_WGS84+"'",
                           " "+src_dataset, 
                           " "+dst_dataset])
                           
        os.system(cmd_str)
        

        #NOW REPROJECT
        
        src_dataset = dst_dataset
        #dst_dataset <- src_dataset
        dst_dataset = os.path.splitext(os.path.basename(src_dataset))[0]+"_projected"+file_format
        dst_dataset = os.path.join(out_dir,dst_dataset)

        #dst_dataset <- paste(sub(file_format,"_projected",basename(src_dataset)),file_format,sep="")
        dst_dataset <- file.path(out_dir,dst_dataset)
        cmd_str = "".join(["gdal_warp",
                           "-t_srs"+" '"+CRS_reg+"'",
                           " "+src_dataset, 
                           " "+dst_dataset])
                           
        os.system(cmd_str)

 
    
  f_list = list.files(path=l_dir[[i]],pattern="*.asc",full.names=T)
  
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


    ################## First step: download data ###################
    # Added tile loop 
    year_list=range(start_year,end_year+1) #list of year to loop through
   
    return None
    
#Need to add this to run
if __name__ == '__main__':
    main()


