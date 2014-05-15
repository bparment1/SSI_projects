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
# Updated on: 05/15/2014

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
import numpy as np

#------------------
# Functions used in the script 
#------------------

def create_dir_and_check_existence(path):
    #Create a new directory
    try:
        os.makedirs(path)
        
    except:
        print "directory already exists"

def calculate_region_extent(reg_outline,out_suffix_dst,CRS_dst,out_dir):
    
    #get extent in WGS84 (lat long)
    src_dataset = reg_outline
    dst_dataset = os.path.splitext(src_dataset)[0]+out_suffix_dst+os.path.splitext(src_dataset)[1]
    dst_dataset = os.path.join(out_dir,dst_dataset)

    #Write a quick wrapper function later on...
    
    cmd_str = "".join(["ogr2ogr",
                      " "+"-f "+"'ESRI Shapefile'",  
                      " "+"-t_srs"+" '"+CRS_dst+"'",
                      " "+"-overwrite",
                      " "+dst_dataset, 
                      " "+src_dataset])
    os.system(cmd_str) #Should try to use subprocess.call with list!!! no need of managing spaces
    
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
    #w_extent formated for use in gdal!!
    return(w_extent,dst_dataset)

def create_raster_region(j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,out_suffix,out_dir,
                         w_extent,clip_param=True,reproject_param=True):
    
    #### Function parameters
    #in_dir: input directory
    #out_dir: output dir
    #in_file: input raster file
    #out_suffix: string
    #CRS_dst: projecction for the output
    #CRS_src: projection for the input
    #file_format: output file format...(tiff as default)
    #w_extent: should be in coordinate system of src!!
    
    ####### START SCRIPT/FUNCTION  #########
    
    in_file= infile_l_f[j]
    #CRS_src = CRS_aea
    #CRS_dst = CRS_reg
    
    ## FIRST CLIP for every file in hte list crop
    if clip_param == True:
        
        #Add out_suffix!!! and turn this into a function
        
        src_dataset = in_file
        #build output name and remove the input path
        dst_dataset = os.path.splitext(os.path.basename(src_dataset))[0]+"_clipped_"+out_suffix+file_format
        dst_dataset = os.path.join(out_dir,dst_dataset)
                            
        cmd_str = "".join(["gdal_translate",
                           " "+"-projwin"+" "+w_extent,
                           " "+"-a_srs"+" '"+CRS_src+"'",
                           " "+src_dataset, 
                           " "+dst_dataset])                     
        os.system(cmd_str)
        
        ### END OF FUNCTION
        clipped_dataset = dst_dataset
        #Also add option to reproject to specificied image of specific resolution
     
    ##SECOND REPROJECT 
    if reproject_param == True:
        #Need to take into account first step
        #prepare input and output
        if clip_param == True:
            src_dataset = clipped_dataset
            dst_dataset = "".join([os.path.splitext(os.path.basename(in_file))[0],
                                   "_clipped_projected_",
                                   out_suffix,
                                   file_format])   
            dst_dataset = os.path.join(out_dir,dst_dataset)
                                   
        elif clip_param == False:
            src_dataset = in_file
            dst_dataset = "".join([os.path.splitext(os.path.basename(in_file))[0],
                                   "_projected_",
                                   out_suffix,
                                   file_format])   
            dst_dataset = os.path.join(out_dir,dst_dataset)                       
    
        cmd_str = "".join(["gdalwarp",
                           " "+"-t_srs"+" '"+CRS_reg+"'",
                           " "+src_dataset, 
                           " "+dst_dataset])               
        os.system(cmd_str)
    
    output_dataset = dst_dataset                
    
    return(output_dataset)
    #end of function


##Change resolution of images
#gdal.ReprojectImage or gdal_warp
#help(gdal.ReprojectImage)

#tto big?? can crash if not enough memory so read by row!!
def find_unique_val_raster(in_raster,band_nb=1):
    ds = gdal.Open(in_raster)
    ncol=ds.RasterYSize
    nrow=ds.RasterXSize
    band = ds.GetRasterBand(band_nb)
    data_type = gdal.GetDataTypeName(band.DataType)

    #ar = numpy.zeros((nrow,ncol),dtype=data_type) #data type 'Byte' not understood
    ar = np.zeros((nrow,ncol))


    ar = np.array(band.ReadAsArray()) #reading everythin in memory!!! needs to be changed!!!
    values = np.unique(ar)
    
    return values
    
#def breakout_raster(in_raster):
#    ...
    
    
def change_resolution_raster(file_format,out_suffix,out_dir,res_val,in_file,out_file=None,resamp_opt="near"):
    #
    src_dataset = in_file
    #resamp_opt = "average" #this will compute the average for given pixel when changing resolution...
    #create output name if out_file=None
    if out_file==None:
        dst_dataset = "".join([os.path.splitext(os.path.basename(in_file))[0],
                                       "_res_",res_val,
                                       out_suffix,
                                       file_format])   
        dst_dataset = os.path.join(out_dir,dst_dataset)                       
    
    cmd_str = "".join(["gdalwarp",
                           " "+"-tr"+" "+res_xy+"'",
                           " "+"-r" +" "+resamp_opt, 
                           " "+src_dataset, 
                           " "+dst_dataset])               
    os.system(cmd_str)
    
    #print "Resampling"
    #outFnLow = "%s/%s_LST_Day_1km_%s_%s_sinu_lowRes.tif" % (outDir,l,monthDict[int(m)],m) 
    #inputStr = ['gdalwarp','-srcnodata',"0","-tr","10000","10000","-dstnodata","0","-wm","500","-overwrite",outFn,outFnLow]
    #out = subprocess.call(inputStr)


    output_dataset = dst_dataset 
    
    return ouput_dataset             

#find number of each category in the new resolution...

#######################################################################
######################### BEGIN SCRIPT  ###############################
#--------------------------------------------
# Script run by arguments from the shell?
#--------------------------------------------
            
def main():
    #
    #--------------------------------------------
    # Clip and reproject data given a shape file or extent
    #--------------------------------------------
 
    ########## READ AND PARSE PARAMETERS AND ARGUMENTS ######### 

    #in_dir = "/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/"
    in_dir ="/ssi-dss-data/DSS_SSI_data" #DSS SSI Maine

    #Input shape file used to define the zonal regions: could be town or counties in this context
    shp_fname = os.path.join(in_dir,"county24.shp")
    #input shp defining study area: can be the same as shp_fname or different                                       
    shp_reg_outline = os.path.join(in_dir,"county24.shp")

    #EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
    #+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
    #CRS_reg = "+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
    #http://spatialreference.org/ref/epsg/2037    
    CRS_reg = "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    CRS_aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 

    file_format = ".tif"
    out_suffix = "05032014"
    w_extent_str = "-72 48 -65 41" #minx,maxy (upper left), maxx,miny (lower right)
    use_reg_extent = True
    os.chdir(in_dir)

    out_dir = in_dir

    #out_path<-"/data/project/layers/commons/data_workflow/output_data"
    out_dir = "output_data_"+out_suffix
    out_dir = os.path.join(in_dir,out_dir)
    create_dir_and_check_existence(out_dir)
            
    os.chdir(out_dir)        #set working directory
    
    ########## START SCRIPT #############
    
    CRS_WGS84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    CRS_dst = CRS_WGS84 #this is the projection system of climate ccsm predictions
    
    #this can be a funciontion !!!! 
    out_suffix_dst = "_WGS84"

    #Read in layers from data source,there is only one layer
    reg_area_poly = ogr.Open(shp_reg_outline).GetLayer()
    
    if use_reg_extent==True:
        w_extent, reg_area_poly_wgs84 = calculate_region_extent(shp_reg_outline,out_suffix_dst,CRS_dst,out_dir)
        #w_extent= "-71.083923738042 47.4598539782516 -66.8854440488051 42.9171281482886"
    elif use_reg_extent==False:
        w_extent = w_extent_str #this is in WGS84
        #end if
     
    ##### PART II: PROCESSING CLIMATE PREDICION  ##########
    
    ### Temp pocessing: 2020s,2030s,2040s,2050s for tmin, tmax and tmean (e.g. 4*3 directories with 12 files...)

    #1.download...not automated...
    #2.unzip asc or grd (arc grid)
    #3.convert to tif (gdal_tranlslate)
    #4.reproject and clip/subset for Maine region (clip using count24?)

    ### EXTRACT ZIPPED FILES ###
    fileglob = "*asc.zip"
    pathglob = os.path.join(in_dir, fileglob)
    l_f = glob.glob(pathglob)
    l_f.sort() #order input by decade
    l_dir = map(lambda x: os.path.splitext(x)[0],l_f) #remmove extension
    l_dir = map(lambda x: os.path.join(out_dir,os.path.basename(x)),l_dir) #set the directory output

    #unzip(l_f[[1]],exdir=l_dir[[1]])
    for i in range(0,len(l_f)):
        with zipfile.ZipFile(l_f[i], "r") as z:
            z.extractall(l_dir[i])
            
    ###  CLIP AND REPROJECT ###
    
    # loop through dir and files
    for i in range(0,len(l_dir)):
        fileglob="*.asc"
        pathglob = os.path.join(l_dir[i], fileglob)
        f_list = glob.glob(pathglob) #sort it!!!
        infile_l_f = f_list
        #This part needs to be turned into a function!!!
        for j in range(0,len(f_list)):
            CRS_src = CRS_WGS84
            CRS_dst = CRS_reg #CRS ofthe study area
            model_info = os.path.basename(os.path.dirname(infile_l_f[0])) #dirname contains the info
            out_suffix_reg = model_info+out_suffix
            
            #outfile = pdb.runcall(create_raster_region,j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,
            #                      out_suffix_reg,out_dir,w_extent,clip_param=True,reproject_param=True)
            outfile = create_raster_region(j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,
                                  out_suffix_reg,out_dir,w_extent,clip_param=True,reproject_param=True)

 
    ##### PART III: PROCESSING LAND COVER  ##########
        
    #TO avoid reproject the large layer (NLCD 30m) we find the extent of the region outile (e.g. ccounties)
    #in the CONUS region for clipping.
    #3x3 and 30x30
    
    out_suffix_dst = "_US_AEA"
    CRS_dst = CRS_aea
    #shp_reg_outline = os.path.join("/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/"
    #                                       ,"county24.shp")

    if use_reg_extent==True:
        w_extent_aea, reg_area_poly_us_aea = calculate_region_extent(shp_reg_outline,out_suffix_dst,CRS_dst,out_dir)
        #w_extent= "-71.083923738042 47.4598539782516 -66.8854440488051 42.9171281482886"
    elif use_reg_extent==False:
        w_extent = w_extent_str #this is in WGS84
        #end if

    f_list <- list.files(path=in_dir,pattern="nlcd.*.img",full.names=T)
    
    fileglob = "*nlcd*.img"
    pathglob = os.path.join(in_dir, fileglob)
    f_list = glob.glob(pathglob)
    #    l_dir = map(lambda x: os.path.splitext(x)[0],l_f) #similar to laplly in R
  
    for j in range(0,len(f_list)):
        CRS_src = CRS_aea #source projection syst
        CRS_dst = CRS_reg  #target projection syst
        w_extent = w_extent_aea
        infile_l_f = f_list
        #outfile = pdb.runcall(create_raster_region,j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,out_suffix,out_dir,w_extent,clip_param=True)
        outfile = create_raster_region(j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,out_suffix,out_dir,w_extent,clip_param=True)
        
            
    ### NOW RECLASSIFY VALUES: i.e. BREAKOUT IDRISI like for unique values (land cover type)
    unique_val = find_unique_val_raster(outfile)
    
    #gdal_calc.py -A crop.tif --outfile=level0100.tif --calc="100*(A>100)"     --NoDataValue=0
    #gdal_calc.py -A crop.tif --outfile=level0100.tif --calc="100*(100<A<200)"     --NoDataValue=0 #reclass between 100 and 200 to 100 else 0
    os.system("gdal_calc.py -A nlcd2001_landcover_v2_2-13-11_clipped_projected_05032014.tif --outfile=ncld_rec11.tif --calc=\"11*(A==11)\"     --NoDataValue=0")
    #create function break out...
    ## THEN use gdal wrap to calcuate the proportion of each land cover types at 100m
    
    ## LAST step use  function to create summary by polygon for counties (proprotion of each land cover)
    
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