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
# Updated on: 03/29/2014

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
    # Clip and reproject data given a shape file or extent
    #--------------------------------------------

    # TODO: set up a (temporary?) GRASS database to use for processing? code
    # currently assumes it's being run within an existing GRASS session
    # using an appropriately configured LOCATION...
    #
 
    ########## START OF SCRIPT ##############
    #### Modified by Benoit on May 13, 2013 
        
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
     
    #this can be a funciontion !!!! 
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

        #w_extent= "-71.083923738042 47.4598539782516 -66.8854440488051 42.9171281482886"
    elif use_reg_extent==False:
        w_extent <- w_extent_str #this is in WGS84
    #end if

    
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
    for i in range(0,len(l_f)):
        fileglob="*.asc"
        pathglob = os.path.join(l_dir[i], fileglob)
        f_list = glob.glob(pathglob) #sort it!!!
        
        #This part needs to be turned into a function!!!
        for j in range(0,len(f_list)):
            
            #Should be made to act only per file so this can be parallelized!!!
            ## FIrst for every file in hte list crop
            #if crop==True:
            src_dataset = os.path.join(in_dir,f_list[j])
            dst_dataset = os.path.splitext(os.path.basename(src_dataset))[0]+file_format
            #Need to attach origingal folder name to keep trackof all types of predictions
            dst_dataset = os.path.basename(l_dir[i])+"_"+dst_dataset
            dst_dataset = os.path.join(out_dir,dst_dataset)
            
            #src_dataset <- file.path(in_dir,f_list[[j]])
            #in this case need to add folder name!!!
            #dst_dataset <- paste(l_dir[[i]],sub(".asc",file_format,basename(src_dataset)),sep="_")
            #dst_dataset <- paste(sub(file_format,"_clipped",basename(src_dataset)),file_format,sep="")

            #dst_dataset <- file.path(out_dir,sub(".asc",file_format,basename(src_dataset)))
            #dst_dataset <- file.path(out_dir,dst_dataset)
            #command_str <-paste("gdal_translate", src_dataset, dst_dataset ,sep=" ") #without spatial subsetting
            
            #Add out_suffix!!! and turn this into a function
            cmd_str = "".join(["gdal_translate",
                               " "+"-projwin"+" "+w_extent,
                               " "+"-a_srs"+" '"+CRS_WGS84+"'",
                               " "+src_dataset, 
                               " "+dst_dataset])
                               
            os.system(cmd_str)
            
            #Seond REPROJECT
            #if reproject ==True
            
            src_dataset = dst_dataset
            #dst_dataset <- src_dataset
            dst_dataset = os.path.splitext(os.path.basename(src_dataset))[0]+"_projected"+file_format
            dst_dataset = os.path.join(out_dir,dst_dataset)
    
            #dst_dataset <- paste(sub(file_format,"_projected",basename(src_dataset)),file_format,sep="")
            #dst_dataset <- file.path(out_dir,dst_dataset)
            cmd_str = "".join(["gdalwarp",
                               " "+"-t_srs"+" '"+CRS_reg+"'",
                               " "+src_dataset, 
                               " "+dst_dataset])
                               
            os.system(cmd_str)

 

    ################## NEXTstep: process landcover data ###################
    # Added tile loop 
   
    #dat_proj <- spTransform(reg_area_poly,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) 

    #reg_area_poly
    #reg_extent<- bbox(dat_proj)
    #reg_extent <- bbox(reg_area_poly_WGS84)
    #w_extent <- c(reg_extent[1,1],reg_extent[2,2],reg_extent[1,2],reg_extent[2,1])
    #w_extent <- paste(as.character(w_extent),collapse=" ")

    #Write a quick wrapper function later on... should have prefix!!!
    src_dataset = os.path.join(in_dir,Maine_counties_file) 
    dst_dataset = os.path.splitext(src_dataset)[0]+"_US_AEA"+os.path.splitext(src_dataset)[1]
    dst_dataset = os.path.join(out_dir,dst_dataset)

    #Make this general provide a target projection sys ...via CRS_proj instead of CRS_WGS84
    CRS_aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
    CRS_dst = CRS_aea
    
    cmd_str = "".join(["ogr2ogr",
                       " "+"-f "+"'ESRI Shapefile'",  
                       " "+"-t_srs"+" '"+CRS_dst+"'",
                       " "+"-overwrite",
                       " "+dst_dataset, 
                       " "+src_dataset])
    os.system(cmd_str)
        
    #reg_area_poly_WGS84 <- spTransform(reg_area_poly,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
    #reg_lyr= ogr.Open(dst_dataset).GetLayer(0)    
    #reg_extent = reg_lyr.GetExtent() #tuple object, not working right now
    #Make this a function?
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

    #w_extent<-project(mat_coord_extent,
#   # proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    #w_extent <- paste(as.character(w_extent),collapse=" ")
    #w_extent <- paste(as.character(bbox(dat)),collapse=" ")

    #as.character(as.vector(t(coordinates(dat))))
    #w_extent <- paste(as.character(t(coordinates(dat))),collapse=" ")

    ### end of function

    f_list <- list.files(path=in_dir,pattern="nlcd.*.img",full.names=T)
    
    fileglob = "*nlcd*.img"
    pathglob = os.path.join(in_dir, fileglob)
    f_list = glob.glob(pathglob)
#    l_dir = map(lambda x: os.path.splitext(x)[0],l_f) #similar to laplly in R

        for j in range(0,len(f_list)):
            #### Function parameters
            #in_dir: input directory
            #out_dir: output dir
            #in_file: input raster file
            #out_suffix: string
            #CRS_dst: projecction for the output
            #CRS_src: projection for the input
            #file_format: output file format...(tiff as default)
            
            in_file=f_list[j]
            CRS_src = CRS_aea
            CRS_dst = CRS_reg
            #Should be made to act only per file so this can be parallelized!!!
            ## FIrst for every file in hte list crop
            #if crop==True:
            src_dataset = os.path.join(in_dir,in_file)
            dst_dataset = os.path.splitext(os.path.basename(src_dataset))[0]+"_clipped_"+out_suffix+file_format
            #Need to attach origingal folder name to keep trackof all types of predictions
            #dst_dataset = os.path.basename(l_dir[i])+"_"+dst_dataset
            dst_dataset = os.path.join(out_dir,dst_dataset)
            
            #src_dataset <- file.path(in_dir,f_list[[j]])
            #in this case need to add folder name!!!
            #dst_dataset <- paste(l_dir[[i]],sub(".asc",file_format,basename(src_dataset)),sep="_")
            #dst_dataset <- paste(sub(file_format,"_clipped",basename(src_dataset)),file_format,sep="")

            #dst_dataset <- file.path(out_dir,sub(".asc",file_format,basename(src_dataset)))
            #dst_dataset <- file.path(out_dir,dst_dataset)
            #command_str <-paste("gdal_translate", src_dataset, dst_dataset ,sep=" ") #without spatial subsetting
            
            #Add out_suffix!!! and turn this into a function
            
            cmd_str = "".join(["gdal_translate",
                               " "+"-projwin"+" "+w_extent,
                               " "+"-a_srs"+" '"+CRS_src+"'",
                               " "+src_dataset, 
                               " "+dst_dataset])
                               
            os.system(cmd_str)
            
            #Seond REPROJECT
            #if reproject ==True
            
            src_dataset = dst_dataset #change in case the first step is not carried out
            #dst_dataset <- src_dataset
            dst_dataset = os.path.splitext(os.path.basename(src_dataset))[0]+"_projected_"+out_suffix+file_format
            dst_dataset = os.path.join(out_dir,dst_dataset)
    
            #dst_dataset <- paste(sub(file_format,"_projected",basename(src_dataset)),file_format,sep="")
            #dst_dataset <- file.path(out_dir,dst_dataset)
            cmd_str = "".join(["gdalwarp",
                               " "+"-t_srs"+" '"+CRS_reg+"'",
                               " "+src_dataset, 
                               " "+dst_dataset])
                               
            os.system(cmd_str)
            
            
            #return
            #end of function

    return None
    
#Need to add this to run
if __name__ == '__main__':
    main()


#ogrinfo clipping_mask.shp -so -al | grep Extent

# which gives the extent as xMin,yMin, xMax, yMax:
#Extent: (268596, 5362330) - (278396, 5376592)
# which is (xMin,yMin) - (xMax,yMax)
#Then copy and paste that text to create your gdal_translate clipping command:

# -projwin's ulx uly lrx lry is equivalent to xMin, yMax, xMax, yMin so just switch the Y coordinates
# For the above Extent that would turn into:

#gdal_translate -projwin 268596 5376592 278396 5362330 src_dataset dst_dataset