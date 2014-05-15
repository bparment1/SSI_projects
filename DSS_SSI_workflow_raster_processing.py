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
    #ar = np.zeros((nrow,ncol))

    ar = np.array(band.ReadAsArray()) #reading everythin in memory!!! needs to be changed!!!
    values = np.unique(ar)
    
    return values
    
## Creating boolean images from categorical rasters...such as NLCD   
def breakout_raster_categories(in_file,out_dir,out_suffix_s,unique_val,file_format,NA_flag_val= -9999,boolean=True):
        #This functions creates an image per unique value input categorical image.
        #The output images can be boolean or have the unique values extracted.
        #This can be used to reclassify images...!!
        
        ### Add loop here
        l_out_file = []
        for i in range(0,len(unique_val)):
            val_in  = str(unique_val[i]) #value to be reclassified
            if boolean==True:
                val_out = str(1)
            if boolean==False:
                val_out = str(unique_val[i]) #value to assign
            #val_out = str(unique_val[0]) #value to assign
            #val = 11
            NA_flag_val_str = str(NA_flag_val)
            #in_file="nlcd2001_landcover_v2_2-13-11_clipped_projected_05152014.tif"
            #out_file = os.path.splitext(os.path.basename(in_file))[0]+"_rec_"+str(val)+file_format
            out_file = os.path.splitext(os.path.basename(in_file))[0]+"_rec_"+val_in+"_"+out_suffix_s+file_format
            out_file = os.path.join(out_dir,out_file)

            #Use this formatting!! this is more efficient than dealing with spaces...and can deal with 
            #additional options later on...
            #cmdStr = ['gdal_calc.py','-A',in_file,"--outfile="+"ncld_rec11.tif","--calc="+"11*(A==11)","--NoDataValue="+"0"]
            cmdStr = ['gdal_calc.py','-A',in_file,"--outfile="+out_file,"--calc="+val_out+"*(A=="+val_in+")","--NoDataValue="+NA_flag_val_str]
            out = subprocess.call(cmdStr)
            l_out_file.append(out_file)
        
        return l_out_file
    
def change_resolution_raster(in_file,res_xy_val,out_suffix_s,out_dir,file_format,output_type=None,NA_flag_val=-9999,out_file=None,resamp_opt="near"):
    #basic command: gdalwarp -tr 10 10 input.tif output.tif
    src_dataset = in_file
    #resamp_opt = "average" #this will compute the average for given pixel when changing resolution...
    #create output name if out_file=None
    if out_file==None:
        dst_dataset = "".join([os.path.splitext(os.path.basename(in_file))[0],
                                       "_res_",
                                       out_suffix_s,
                                       file_format])   
        dst_dataset = os.path.join(out_dir,dst_dataset)                       
    if out_file!=None:
        dst_dataset = out_file
    #resamp_opt = "average"
    NA_flag_val_str = str(NA_flag_val)

    if ouput_type==None:
            cmd_str = ["gdalwarp",
              "-tr",str(res_xy_val[0]),str(res_xy_val[1]), 
              "-r",resamp_opt, 
              "-dstnodata", NA_flag_val_str,
              "-overwrite",
              src_dataset, 
              dst_dataset]     

    if ouput_type!=None:
            cmd_str = ["gdalwarp",
              "-tr",str(res_xy_val[0]),str(res_xy_val[1]), 
              "-r",resamp_opt, 
              "-ot",output_type,
              "-dstnodata", NA_flag_val_str,
              "-overwrite",
              src_dataset, 
              dst_dataset]
            
    print "Resampling :changing resolution"
    
    out = subprocess.call(cmd_str)   
    #os.system(cmd_str)
    
    out_file = dst_dataset 
    
    return ouput_file             

def raster_calc_operation_on_list(l_rast,out_dir,out_suffix_s,file_format,operation="+",NA_flag_val= -9999,out_file=None):
            
    ### Add loop here
    l_out_file = []
    NA_flag_val_str = str(NA_flag_val)
    
    if out_file==None:
        #create an output filename first:
        file_combined = os.path.basename(l_rast[0]) #assuming this is the same file type? should change
        #To avoid having a suffix attached to a name several times remove and replace by new suffix...
        file_combined = file_combined.replace(out_suffix_s+file_format,"_calc_combined"+out_suffix_s+file_format)
        file_combined = os.path.join(out_dir,file_combined)
        in_file_1 = os.path.basename(l_rast[0])        
    if out_file!=None:
        #create an output filename first:
        file_combined = os.path.join(out_dir,out_file) #assuming this is the same file type? should change
        in_file_1 = os.path.basename(l_rast[0])        
  
    cmdStr = ['gdal_calc.py',
                  '-A',in_file_1,
                  "--outfile="+file_combined,
                  "--calc=1*(A*0)",
                  "--NoDataValue="+NA_flag_val_str]
    out = subprocess.call(cmdStr)
  
    #start at zero but length is reduced by one!!!
    for i in range(0,len(l_rast)-1):
        #val_out = str(unique_val[0]) #value to assign
        #val = 11
        in_file_1 = file_combined
        if i==0:
            in_file_1 = l_rast[0]
            
        out_file = file_combined
        in_file_2 = l_rast[i+1]
        #out_file = os.path.join(out_dir,"tmp_combined"+ file_format)
        #Use this formatting!! this is more efficient than dealing with spaces...and can deal with 
        #additional options later on...
        #cmdStr = ['gdal_calc.py','-A',in_file,"--outfile="+"ncld_rec11.tif","--calc="+"11*(A==11)","--NoDataValue="+"0"]
        cmdStr = ['gdal_calc.py',
                  '-A',in_file_1,
                  '-B',in_file_2,
                  "--outfile="+out_file,
                  "--calc=1*(A"+operation+"B)",
                  "--NoDataValue="+NA_flag_val_str]
        out = subprocess.call(cmdStr)
        file_combined = out_file
        #l_out_file.append(out_file)
    
    return file_combined


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
    CRS_WGS84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

    file_format = ".tif"
    NA_flag_val = -9999
    out_suffix = "05152014"
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
     
    ##### PART I: PROCESSING CLIMATE PREDICION  ##########
    
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

 
    ##### PART II: PROCESSING LAND COVER  ##########
        
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
   
    fileglob = "*nlcd*.img"
    pathglob = os.path.join(in_dir, fileglob)
    f_list = glob.glob(pathglob)
    #    l_dir = map(lambda x: os.path.splitext(x)[0],l_f) #similar to laplly in R
    
    outfile_list = []
    outfile_breakout_list = []
    
    for j in range(0,len(f_list)):
        CRS_src = CRS_aea #source projection syst
        CRS_dst = CRS_reg  #target projection syst
        w_extent = w_extent_aea
        infile_l_f = f_list
        #outfile = pdb.runcall(create_raster_region,j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,out_suffix,out_dir,w_extent,clip_param=True)
        outfile = create_raster_region(j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,out_suffix,out_dir,w_extent,clip_param=True)
        outfile_list.append(outfile)
        ### NOW RECLASSIFY VALUES: i.e. BREAKOUT IDRISI like for unique values (land cover type)
        unique_val = find_unique_val_raster(outfile) #This is inefficient! bad function!! rewrite!!
        outfile_breakout_nlcd = breakout_raster_categories(outfile,out_dir,out_suffix,unique_val,file_format,NA_flag_val= -9999,boolean=True)
        outfile_breakout_list.append(outfile_breakout_nlcd)
        
    #nlcd92mosaic_clipped_projected_05152014_rec_0_05152014    
     outfile_breakout_list = []
     outfile_breakout_list.append(glob.glob(os.path.join(out_dir,"*nlcd92mosaic*rec*.tif")))
     outfile_breakout_list.append(glob.glob(os.path.join(out_dir,"*nlcd2001*rec*.tif")))
     outfile_breakout_list.append(glob.glob(os.path.join(out_dir,"*nlcd2006*rec*.tif")))

     ## Need to reclassify ?
     #Combine all forest, urban, agri, wetland...?
     #NLCD 2001,2006
     20: urban , Developed (all 21,22,23,24)
     40: forest, Forested upland (all 41,41,43)
     80: agriculture, planted/cultivated (all,81,82)
     90: wetland, (all,91,92,93...)

     ## Clean this up to make it more general and shorter...
     l_f_nlcd92_subset_cat = []
     l_f_nlcd92_subset_cat.append(filter(lambda x: re.search(r'rec_2',x),outfile_breakout_list[0]))        
     l_f_nlcd92_subset_cat.append(filter(lambda x: re.search(r'rec_4',x),outfile_breakout_list[0]))        
     l_f_nlcd92_subset_cat.append(filter(lambda x: re.search(r'rec_8',x),outfile_breakout_list[0]))        
     l_f_nlcd92_subset_cat.append(filter(lambda x: re.search(r'rec_9',x),outfile_breakout_list[0]))        

     l_f_nlcd2001_subset_cat = []
     l_f_nlcd2001_subset_cat.append(filter(lambda x: re.search(r'rec_2',x),outfile_breakout_list[1]))        
     l_f_nlcd2001_subset_cat.append(filter(lambda x: re.search(r'rec_4',x),outfile_breakout_list[1]))        
     l_f_nlcd2001_subset_cat.append(filter(lambda x: re.search(r'rec_8',x),outfile_breakout_list[1]))        
     l_f_nlcd2001_subset_cat.append(filter(lambda x: re.search(r'rec_9',x),outfile_breakout_list[1]))        

     l_f_nlcd2006_subset_cat = []
     l_f_nlcd2006_subset_cat.append(filter(lambda x: re.search(r'rec_2',x),outfile_breakout_list[2]))        
     l_f_nlcd2006_subset_cat.append(filter(lambda x: re.search(r'rec_4',x),outfile_breakout_list[2]))        
     l_f_nlcd2006_subset_cat.append(filter(lambda x: re.search(r'rec_8',x),outfile_breakout_list[2]))        
     l_f_nlcd2006_subset_cat.append(filter(lambda x: re.search(r'rec_9',x),outfile_breakout_list[2]))        
        
     l_f_nlcd_subset_cat = {"nlcd92": l_f_nlcd1992_subset_cat, "nlcd2001": l_f_nlcd2001_subset_cat, "nlcd2006": l_f_nlcd2006_subset_cat }

     for i in range(0,len(l_f_nlcd_subset_cat)):
         l_f_nlcd =l_f_nlcd_subset_cat[i]
         
         for j in range(0,len(l_f_nlcd)):
             out_file = "test.tif" 
             nlcd2006_dss_cat = raster_calc_operation_on_list(l_f_nlcd2006_subset_cat[0],out_dir,out_suffix_s,file_format,operation="+",NA_flag_val= -9999,out_file)
               
        
    #This specific to nlcd    
    for j in range(0,len(outfile_breakout_list):
        outfile_breakout_nlcd = outfile_breakout_list[j]
        
        for i in range(0,len( outfile_breakout_nlcd)):
            nlcd_cat_file = outfile_breakout_nlcd[i]
            res_xy_val = [30,30] #this is for visualization...
            change_resolution_raster(nlcd_cat_file,res_xy_val,out_suffix_s,out_dir,file_format,output_type=None,NA_flag_val=-9999,out_file=None,resamp_opt="near"):
            res_xy_val = [300,300]
            change_resolution_raster(nlcd_cat_file,res_xy_val,out_suffix_s,out_dir,file_format,output_type=None,NA_flag_val=-9999,out_file=None,resamp_opt="near"):

    #generate new resolution...
    ## first at 3x30m for display (per land cover categories) use gdal wrap to calcuate the proportion of each land cover types at 90m
    ## second at 30x30m for calculations of proportions in POSTGIS

    ## LAST step use  function to create summary by polygon for counties (proprotion of each land cover)
    ### The last step takes outside this script
    #add look up table for legend? It wll be important for the name.
    
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