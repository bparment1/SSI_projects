#!/usr/bin/python
#
######## PROCESS RASTER IMAGES FOR STUDY AREA  ########
#
#Script to process data for the Decision Support System, DSS-SSI website.
#It is aimed at general purpose (i.e. script must strives for generality and automation)
#Data are first clipped and reprojected to match the Maine study area.
#NLCD data are reclassified into general classes (forest,agriculture,urban and wetland)
#NLCD general classes are resampled at higher resolution and proportion of land cover  calculated.
#
## TODO:
#  - Add paralelization...
#  - get statistics for layers in tif
#  - create and apply mask for ME to all outputs
#  - solve issues with nlcd_2001_urban_proportion_900_900_05152014.tif and 90m
#  - find out about palette storage
#  - change projection to 26919?
#  - divide by 10 the temp...add option in region function??
#
# Authors: Benoit Parmentier 
# Created on: 03/24/2014
# Updated on: 06/11/2014
# Project: DSS-SSI
#
####### LOAD LIBRARY/MODULES USED IN THE SCRIPT ###########
#    
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
from multiprocessing import Process, Manager, Pool #parallel processing
import pdb                   #for debugging
import pandas as pd          #DataFrame object and other R like features for data munging


################ NOW FUNCTIONS  ###################

##------------------
# Functions used in the script 
##------------------

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

def create_raster_region(j,in_dir,infile_l_f,CRS_dst,file_format,out_suffix,out_dir,
                         w_extent,CRS_src=None,NA_flag_val=None,output_type=None,clip_param=True,reproject_param=True):
    
    #### Function parameters
    #in_dir: input directory
    #out_dir: output dir
    #in_file: input raster file
    #out_suffix: string
    #CRS_dst: projecction for the output
    #CRS_src: projection for the input
    #file_format: output file format...(tiff as default)
    #w_extent: should be in coordinate system of src!!
    #NA_flag_val: value for the Nodata pixels
    #output_type: user defined such as "Float32", if None, the same datatype is kept
    
    ####### START SCRIPT/FUNCTION  #########
    
    in_file= infile_l_f[j]
    #CRS_src = CRS_aea
    #CRS_dst = CRS_reg
    
    ds = gdal.Open(in_file)
    band = ds.GetRasterBand(1)
    data_type = gdal.GetDataTypeName(band.DataType)
    no_data_val = band.GetNoDataValue()
    
    Projection = osr.SpatialReference() #create a projection object
    EPSG_code= Projection.GetAuthorityCode(None)
    Projection.ImportFromWkt(ds.GetProjectionRef())
    proj_str = Projection.ExportToProj4() #show the format of the CRS projection object

    if NA_flag_val==None:
        NA_flag_val_str = str(no_data_val)
    if NA_flag_val!=None:
        NA_flag_val_str = str(NA_flag_val)
    if CRS_src==None:
        CRS_src = proj_str
    
    ## FIRST CLIP for every file in hte list crop
    if clip_param == True:
        
        #Add out_suffix!!! and turn this into a function
        
        src_dataset = in_file
        #build output name and remove the input path
        dst_dataset = os.path.splitext(os.path.basename(src_dataset))[0]+"_clipped_"+out_suffix+file_format
        dst_dataset = os.path.join(out_dir,dst_dataset)
        
        if output_type==None:
            cmd_str = "".join(["gdal_translate",
                               " "+"-projwin"+" "+w_extent,
                               " "+"-a_srs"+" '"+CRS_src+"'",
                               " -a_nodata "+NA_flag_val_str,
                               " "+src_dataset, 
                               " "+dst_dataset])            
        if output_type!=None:
            cmd_str = "".join(["gdal_translate",
                               " -ot "+output_type,
                               " "+"-projwin"+" "+w_extent,
                               " "+"-a_srs"+" '"+CRS_src+"'",
                               " -a_nodata "+NA_flag_val_str,
                               " "+src_dataset, 
                               " "+dst_dataset])                     

        os.system(cmd_str)
        #-a_nodata value: no data value....
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
        
        if output_type==None:
            cmd_str = "".join(["gdalwarp",
                               " "+"-t_srs"+" '"+CRS_reg+"'",
                               " -srcnodata "+NA_flag_val_str,
                               " -dstnodata "+NA_flag_val_str,
                               " -overwrite",
                               " "+src_dataset, 
                               " "+dst_dataset])             

        if output_type!=None:
            cmd_str = "".join(["gdalwarp",
                               " "+"-t_srs"+" '"+CRS_reg+"'",
                               " -ot "+output_type,
                               " -srcnodata "+NA_flag_val_str,
                               " -dstnodata "+NA_flag_val_str,
                               " -overwrite",
                               " "+src_dataset, 
                               " "+dst_dataset])             
                           
                           
        os.system(cmd_str)
    #can also add srcnodata
    output_dataset = dst_dataset                
    band = None
    ds = None
    return(output_dataset)
    #end of function


##Change resolution of images
#gdal.ReprojectImage or gdal_warp
#help(gdal.ReprojectImage)

#tto big?? can crash if not enough memory so read by row!!
#Ok changed...reading rows by row binary values but quite  slow...maybe read by 2 rows???
#this currently taks aobut 5minutes for ths large 19913*15583 or about 310 millions of pixels in nlcd raster files
#At later stage test np.memmap and paralelization? Sure...
def find_unique_val_raster(in_raster,band_nb=1):
    #Function  uses the classic approach for large raster images...ie. read by chunk
    #and process rows by rows simiar to IDRISI delphi code.
    #in_raster = "county24reg_06042014_rast_06042014.tif"
    #in_raster = "nlcd2006_landcover_4-20-11_se5_clipped_projected_06042014.tif"
    ds = gdal.Open(in_raster)
    nrow = ds.RasterYSize
    ncol =ds.RasterXSize        

    band = ds.GetRasterBand(band_nb)
    data_type = gdal.GetDataTypeName(band.DataType)
    block_sizes = band.GetBlockSize()  #this is equal to number of rows in a Geo
    no_data_val = band.GetNoDataValue()


    list_unique_val = []
    for i in range (nrow):
        #data = band.ReadAsArray(0,i,nrow,1)
        #data = band.ReadAsArray(0,1,i,ncol)
        #data = band.ReadAsArray(1,0,1,ncol) #read first row  
        data = band.ReadAsArray(0,i,ncol,1) #read by row and place in a buffer array        
        val = np.unique(data)
        list_unique_val= list_unique_val+(val.tolist())
        unique_val = list(set(list_unique_val))
    #img = np.memmap(rast_fname, dtype=np.float32, shape=(19913, 15583))? could be an alternative using
    #memory mapping for large numpy array??
    #ar = numpy.zeros((nrow,ncol),dtype=data_type) #data type 'Byte' not understood
    #ar = np.zeros((nrow,ncol))
    #ar = np.array(band.ReadAsArray()) #reading everythin in memory!!! needs to be changed!!!
    #ar = band.ReadAsArray() #reading everythin in memory!!! needs to be changed!!!
    #values = np.unique(ar)
    #values = np.array(unique_val)
    #values = np.where(values!= no_data_val)
    band = None
    ds = None
    return unique_val
    
## Creating boolean images from categorical rasters...such as NLCD   
def breakout_raster_categories(in_file,out_dir,out_suffix_s,unique_val,file_format,NA_flag_val= -9999,boolean=True):
    
    #This functions creates an image per unique value input categorical image.
    #The output images can be boolean or have the unique values extracted.
    #This can be used to reclassify images...!!
        
    ### Add loop here
    #should only convert integer val to float if the inputfile is float!!!
    #
    ds = gdal.Open(in_file)
    band = ds.GetRasterBand(1)
    data_type = gdal.GetDataTypeName(band.DataType)
    no_data_val = band.GetNoDataValue()

    l_out_file = []
    for i in range(0,len(unique_val)):
        #d
        #val_in  = str(unique_val[0][i]) #value to be reclassified   
        if data_type == "Float32":
            val_in  = str(float(unique_val[i])) #value to be reclassified
        if data_type == "UInt32":
            val_in  = str(unique_val[i]) #value to be reclassified
        if data_type == "Byte":
            val_in  = str(unique_val[i]) #value to be reclassified
        if boolean==True:
            val_out = str(1)
        if boolean==False:
            val_out = str(val_in) #value to assign
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
        #cmd_str = "".join(["gdal_calc.py",
        #                  # " "+"-projwin"+" "+w_extent,
        #                " "+"-A "+in_file,
        #                " --outfile="+out_file,
        #                "--calc="+val_out+"*(A=="+val_in+")",
        #                "--NoDataValue="+NA_flag_val_str])       
        #os.system(cmd_str)
        
        l_out_file.append(out_file)
        
        band = None
        ds = None
    return l_out_file

## Creating boolean images from categorical rasters...such as NLCD   
def reclass_raster_categories(in_file,out_dir,out_suffix_s,unique_val,out_val,file_format,NA_flag_val= -9999):
        #This functions creates an image per unique value input categorical image.
        #The output images can be boolean or have the unique values extracted.
        #This can be used to reclassify images...!!
        
        ### Add loop here
        l_out_file = []
        for i in range(0,len(unique_val)):
            val_in  = str(unique_val[i]) #value to be reclassified
            val_out = str(out_val[i]) #value to assign
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
            cmd_str = "".join(["gdal_calc.py",
                               #" -at ", #all touched pixel are converted...
                               " "+"-a "+att_field,
                               " -l "+layer_name,
                               " -a_srs "+"'"+proj_str+"'",
                               " -a_nodata "+NA_flag_val_str,
                               " -tr "+r_res_str,
                               " -te "+r_extent_str,
                               " -ot "+output_type,
                               " "+src_dataset, 
                               " "+dst_dataset])        

            l_out_file.append(out_file)
        
        return l_out_file

## Creating mask image from rasterinput: this creates a boolean image from input layer   
def create_raster_mask(in_file,out_dir,out_suffix_s,file_format,CRS_reg,NA_flag_val= -9999,out_file=None):
    #This functions creates an mask from an input raster.
    #ALl values that are not NA will be reclassified as 1 and all other as NA.
    #this function will be improve later on...  
      
    NA_flag_val_str = str(NA_flag_val)

    NA_flag_val_in = getNoDataValue(in_file) #get no data val from file...
    if out_file==None:
        out_file = os.path.splitext(os.path.basename(in_file))[0]+"_masked_"+out_suffix_s+file_format
        out_file = os.path.join(out_dir,out_file)

    #additional options later on...
    val_out = str(1)
    val_in = str(NA_flag_val_in)
    tmp_file = "tmp_"+out_suffix_s+file_format #temporary file, should clean this up...add later...
    #cmdStr = ['gdal_calc.py','-A',in_file,"--outfile="+"ncld_rec11.tif","--calc="+"11*(A==11)","--NoDataValue="+"0"]
    cmdStr = ['gdal_calc.py','-A',in_file,"--outfile="+tmp_file,"--calc="+val_out+"*(A!="+val_in+")","--NoDataValue="+NA_flag_val_str,"--overwrite"]
    out = subprocess.call(cmdStr)
    
    src_dataset = tmp_file
    dst_dataset = out_file
    CRS_src = CRS_reg
    cmd_str = "".join(["gdal_translate",
                              # " "+"-projwin"+" "+w_extent,
                            " "+"-a_srs"+" '"+CRS_src+"'",
                            " -a_nodata "+NA_flag_val_str,
                            " "+src_dataset, 
                            " "+dst_dataset])            
    os.system(cmd_str)
    
    return out_file

def apply_raster_mask(in_file,mask_file,out_dir,out_suffix_s,file_format,NA_flag_val= -9999,EPSG_code=None,out_file=None):
    #This functions creates an mask from an input raster.
    #ALl values that are not NA will be reclassified as 1 and all other as NA.
    #this function will be improve later on...  
      
    NA_flag_val_str = str(NA_flag_val)
    #NA_flag_val_in = getNoDataValue(in_file) #get no data val from file...
    if out_file==None:
        out_file = os.path.splitext(os.path.basename(in_file))[0]+"_masked_"+out_suffix_s+file_format
        out_file = os.path.join(out_dir,out_file)
        
    #out_suffix_s = "_masked_"+out_suffix
    #out_file = lf_temp[i]
    #in_file = lf_temp[i]
    #mask_file = mask_rast_file
    out_file = out_file.replace(out_suffix,out_suffix_s) #remove the suffix if it is there in the file name
        
    if EPSG_code!=None:
        #d
        #d
        tmp_file = "tmp.tif"
        src_dataset = in_file
        dst_dataset = tmp_file
        CRS_src = "EPSG:"+EPSG_code
        cmd_str = "".join(["gdal_translate",
                              # " "+"-projwin"+" "+w_extent,
                            " "+"-a_srs"+" '"+CRS_src+"'",
                            " -a_nodata "+NA_flag_val_str,
                            " "+src_dataset, 
                            " "+dst_dataset])            
        os.system(cmd_str)
        in_file=tmp_file
       #
        tmp_mask_file = "tmp_mask.tif"
        src_dataset = mask_file
        dst_dataset = tmp_mask_file
        CRS_src = "EPSG:"+EPSG_code
        cmd_str = "".join(["gdal_translate",
                              # " "+"-projwin"+" "+w_extent,
                            " "+"-a_srs"+" '"+CRS_src+"'",
                            " -a_nodata "+NA_flag_val_str,
                            " "+src_dataset, 
                            " "+dst_dataset])            
        os.system(cmd_str)
        mask_file=tmp_mask_file
       
    #additional options later on...
    #cmdStr = ['gdal_calc.py','-A',in_file,"--outfile="+"ncld_rec11.tif","--calc="+"11*(A==11)","--NoDataValue="+"0"]   
    #mask_file =lf_temp[2]
    #out_file="test.tif
    #cmdStr = ['gdal_calc.py',
    #              '-A',mask_file,
    #              "--outfile="+out_file,
    #              "--calc=1*(A*1)",
    #              "--NoDataValue="+NA_flag_val_str]
    #out = subprocess.call(cmdStr)
    out_file = os.path.basename(out_file)
    #cmdStr = ['gdal_calc.py','-A',in_file,'-B',mask_file,"--outfile="+out_file,"--calc=1*(A*B)","--NoDataValue="+NA_flag_val_str]
    #cmdStr = ['gdal_calc.py','-A',in_file,'-B',mask_file,"--outfile="+out_file,"--calc="+"'(A*B)'","--NoDataValue="+NA_flag_val_str]
    #out = subprocess.call(cmdStr)
    #out = subprocess.call(cmdStr)
    cmd_str = "".join(["gdal_calc.py",
                              # " "+"-projwin"+" "+w_extent,
                            " "+"-A "+in_file,
                            " "+"-B "+"tmp_mask.tif",
                            " --outfile="+out_file,
                            " --calc="+"'(A*B)'"])       
    os.system(cmd_str)
    
    return out_file
    
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

    if output_type==None:
            cmd_str = ["gdalwarp",
              "-tr",str(res_xy_val[0]),str(res_xy_val[1]), 
              "-r",resamp_opt, 
              "-dstnodata", NA_flag_val_str,
              "-overwrite",
              src_dataset, 
              dst_dataset]     

    if output_type!=None:
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
    
    return out_file             

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

## ADD new functions...
#  - get statistics for layers in tif
#  - create and apply mask for ME to all outputs
#http://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
def getNoDataValue(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    return band.GetNoDataValue()
    
def raster_to_poly_operation_on_list(in_vect,in_rast,out_suffix_s,att_field,file_format,output_type,all_touched=False,NA_flag_val= -9999,CRS_reg=None,out_file=None):
#gdal_rasterize -a ROOF_H -where 'class="A"' -l footprints footprints.shp city_dem.tif

    NA_flag_val_str = str(NA_flag_val) #note that if integer NA is not set to -9999?
    layer_name = os.path.splitext(os.path.basename(in_vect))[0]
    
    #This should be an option....
    ds = gdal.Open(in_rast)
    ncol=ds.RasterYSize
    nrow=ds.RasterXSize
    #band = ds.GetRasterBand(band_nb)
    band = ds.GetRasterBand(1)
    data_type = gdal.GetDataTypeName(band.DataType)
   
    ##Make a general function to get info from  raster??
    # Get raster georeference info: the extent of the raster!
    if CRS_reg==None:        
        Projection = osr.SpatialReference()
        EPSG_code= Projection.GetAuthorityCode(None)
        Projection.ImportFromWkt(ds.GetProjectionRef())
        proj_str = Projection.ExportToProj4() #show the format of the CRS projection object
    #Problem even though the EPSG code is 26919 the proj4 exported does the not the input:
   #
    if CRS_reg!=None:
        proj_str =CRS_reg
        
    geoTransform = ds.GetGeoTransform()
    xmin = geoTransform[0] #Topleft x
    ymax = geoTransform[3] #Topleft y
    xmax = xmin + geoTransform[1]*ds.RasterXSize
    ymin = ymax + geoTransform[5]*ds.RasterYSize
    xres = geoTransform[1] # w-e pixel res
    yres = geoTransform[5] # n-s pixel res
    
    r_extent = [xmin, ymin, xmax, ymax]
    r_extent_str = ' '.join(map(str, r_extent))
    
    r_res = [abs(xres), abs(yres)] 
    r_res = [xres, yres] 
    
    r_res_str = ' '.join(map(str, r_res))    

    if out_file==None:
        #create an output filename first:
        #To avoid having a suffix attached to a name several times remove and replace by new suffix...
        out_file = layer_name+"_rast_"+out_suffix_s+file_format
        out_file = os.path.join(out_dir,out_file)
    if out_file!=None:
        #create an output filename first:
        out_file = os.path.join(out_dir,out_file) #assuming this is the same file type? should change

    #Add a option later to create a raster without having to give another  raster  image!!
    #att_field ="CNTYCODE" #could get get the type from the field?
    #CNTYCODE
    #output_type = "UInt32"
    #output_type = "Float32"
    #if all_touched == False:
    src_dataset = in_vect
    dst_dataset = out_file
    #dst_dataset ="test3.tif"
    #os.system("gdal_rasterize -a CNTYCODE -l county24 -te 326674.8647578625 4755744.830938448 671552.8095416913 5262918.279149961 -tr 4057.387585692103 -4057.387585692103 /ssi-dss-data/DSS_SSI_data/county24.shp  test2.tif") 
    #326674.8647578625, 4755744.830938448, 671552.8095416913, 5262918.279149961    
    #os.system("gdal_rasterize -a CNTYCODE -l county24 -te 326674.8647578625 4755744.830938448 671552.8095416913 5262918.279149961 -tr 4057.387585692103 -4057.387585692103 /ssi-dss-data/DSS_SSI_data/county24.shp  test2.tif") 
    if all_touched==True:
        cmd_str = "".join(["gdal_rasterize",
                               " -at ", #all touched pixel are converted...
                               " "+"-a "+att_field,
                               " -l "+layer_name,
                               " -a_srs "+"'"+proj_str+"'",
                               " -a_nodata "+NA_flag_val_str,
                               " -tr "+r_res_str,
                               " -te "+r_extent_str,
                               " -ot "+output_type,
                               " "+src_dataset, 
                               " "+dst_dataset])        
    if all_touched==False:                               
        cmd_str = "".join(["gdal_rasterize",
                               #" -at ", #all touched pixel are converted...
                               " "+"-a "+att_field,
                               " -l "+layer_name,
                               " -a_srs "+"'"+proj_str+"'",
                               " -a_nodata "+NA_flag_val_str,
                               " -tr "+r_res_str,
                               " -te "+r_extent_str,
                               " -ot "+output_type,
                               " "+src_dataset, 
                               " "+dst_dataset])        
                               
    os.system(cmd_str)
     
    #cmd_str = "".join(["gdal_translate",
    #                          # " "+"-projwin"+" "+w_extent,
    #                           " "+"-a_srs"+" '"+CRS_src+"'",
    #                           " -a_nodata "+NA_flag_val_str,
    #                           " "+src_dataset, 
    #                          " "+dst_dataset])            
   
    #not assigning proj_str properply so...
    
    #Problem with cmdStr and subprocess.. find out later...works for now
    #cmdStr = ['gdal_rasterize',
    #         #"-b", "1",
    #        "-a" ,att_field,
    #       "-l", layer_name,
     #       #"-of", file_format,
    #          "-a_srs", proj_str,
    #          "-a_no_data_value",NA_flag_val_str,
    #           "-te ", "xmin ymin xmax ymax", #output extent , use the one already defined earlier!! use tap?
    #          "-te ", r_extent_str, #output extent , use the one already defined earlier!! use tap?
    #          #"-tr ", "xres yres",
     #         "-tr ", r_res_str,
     #         "-ot", output_type,
     #         src_dataset, 
      #       dst_dataset]
     #out = subprocess.call(cmdStr)
    return dst_dataset
   
def calculate_raster_stat(rasterfn,generate_xml=False):
    ds = gdal.Open(rasterfn)
    band = ds.GetRasterBand(1)

    stats = band.ComputeStatistics(1)
    band = None
    ds = None
    #stats contains: min , max, mean ,std
    if generate_xml==True:
        cmd_str = "gdalinfo -stats "+rasterfn
        os.system(cmd_str)
        
    return stats

### Create a custom function to plot image...
#import matplotlib.image as mpimg
#img=mpimg.imread('MARBLES.TIF ')
#imgplot = plt.imshow(img)
#if 
#rast_fname, to create an array not read in memory use memmap (memory mapping)!!!!!
#img = np.memmap(rast_fname, dtype=np.float32, shape=(19913, 15583))

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
    #CRS_reg = "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    CRS_reg = "+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" #using EPSG 26919
    CRS_reg = "+proj=utm +zone=19 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 

    CRS_aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
    CRS_WGS84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

    file_format = ".tif"
    NA_flag_val = -9999
    output_type = "Float32"
    out_suffix = "06042014"
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

    ### EXTRACT THE ASC ZIPPED FILES ###
    
    fileglob = "*asc.zip"
    pathglob = os.path.join(in_dir, fileglob)
    l_f = glob.glob(pathglob)
    l_f.sort() #order input by decade
    l_dir = map(lambda x: os.path.splitext(x)[0],l_f) #remmove extension
    l_dir = map(lambda x: os.path.join(out_dir,os.path.basename(x)),l_dir) #set the directory output

    #unzip(l_f[[1]],exdir=l_dir[[1]])
    #This part takes a lot of time!!!
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
            #outfile = create_raster_region(j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,
            #                      out_suffix_reg,out_dir,w_extent,NA_flag_val,output_type,clip_param=True,reproject_param=True)
            outfile = create_raster_region(j,in_dir,infile_l_f,CRS_dst,file_format,
                                  out_suffix,out_dir,w_extent,CRS_src,NA_flag_val=None,output_type=None,clip_param=True,reproject_param=True)


    #lf_temp = glob.glob(*)
    #fileglob_pattern = "*projected*ncar*.tif"
    #pathglob = os.path.join(out_dir, fileglob_pattern)
    lf_temp = glob.glob(os.path.join(out_dir, "*projected*ncar*.tif")) #this contains the raster variable files that need to be summarized
    lf_temp.sort()
    
    #Now Apply mask :line
     
    out_suffix_reg = "reg_"+out_suffix #the region file projected in the defined projectio 
    w_extent_reg, reg_area_poly_projected = calculate_region_extent(shp_fname,out_suffix_reg,CRS_reg,out_dir)
    
    #first create mask from region definition:
    in_vect=reg_area_poly_projected
    in_rast=lf_temp[0]
    #out_suffix_s
    att_field="CNTYCODE"
    #file_format
    output_type="Float32"
    all_touched=True
    #NA_flag_val= -9999
    out_file=None
        
    #region_rast_fname = raster_to_poly_operation_on_list(in_vect,in_rast,out_suffix_s,att_field,file_format,output_type,all_touched,NA_flag_val,out_file)
    region_rast_fname = raster_to_poly_operation_on_list(in_vect,in_rast,out_suffix,att_field,file_format,output_type,all_touched,NA_flag_val,CRS_reg,out_file)

    # now either reclass or apply directly the mask
    #This functions creates an mask from an input raster.
    #ALl values that are not NA will be reclassified as 1 and all other as NA.
    #this function will be improve later on...  

    in_file = region_rast_fname
    out_file = "mask_regions_"+out_suffix+file_format
    #mask_rast_file = create_raster_mask(in_file,ouout_dir,out_suffix_s,file_format,NA_flag_val,out_file)
    mask_rast_file = create_raster_mask(in_file,out_dir,out_suffix,file_format,CRS_reg,NA_flag_val,out_file)
    
    #problem with the mask does not contain full
    lf_temp_masked = []
    for i in range(0,len(lf_temp)):
        out_suffix_s = "_masked_"+out_suffix
        out_file = lf_temp[i]
        in_file = lf_temp[i]
        mask_file = mask_rast_file
        out_file = out_file.replace(out_suffix,out_suffix_s) #remove the suffix if it is there in the file name
        f_masked = apply_raster_mask(in_file,mask_file,out_dir,out_suffix_s,file_format,NA_flag_val,out_file)
        lf_temp_masked.append(f_masked)
        #Now get stat
        
    #last step...recalculate stat?? maybe also put in the mask function...

    ##################################################    
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
        #outfile = create_raster_region(j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,out_suffix,out_dir,w_extent,clip_param=True)
        #output_type="UInt32"
        output_type = None   #does not work with real number  assignment...     
        #outfile = create_raster_region(j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,
        #                          out_suffix,out_dir,w_extent,NA_flag_val,output_type,clip_param=True,reproject_param=True)
        outfile = create_raster_region(j,in_dir,infile_l_f,CRS_dst,file_format,
                                  out_suffix,out_dir,w_extent,CRS_src,NA_flag_val=None,output_type=None,clip_param=True,reproject_param=True)

        outfile_list.append(outfile)
        ### NOW RECLASSIFY VALUES: i.e. BREAKOUT IDRISI like for unique values (land cover type)
        unique_val = find_unique_val_raster(outfile) #This has been modified processing by blocks (rows)
        #val=unique_val[0]
        
        outfile_breakout_nlcd = breakout_raster_categories(outfile,out_dir,out_suffix,unique_val,file_format,NA_flag_val= -9999,boolean=True)
        outfile_breakout_list.append(outfile_breakout_nlcd)
        



    #nlcd92mosaic_clipped_projected_05152014_rec_0_05152014    
    #outfile_breakout_list = []
    #outfile_breakout_list.append(glob.glob(os.path.join(out_dir,"*nlcd92mosaic*rec*.tif")))
    #outfile_breakout_list.append(glob.glob(os.path.join(out_dir,"*nlcd2001*rec*.tif")))
    #outfile_breakout_list.append(glob.glob(os.path.join(out_dir,"*nlcd2006*rec*.tif")))

    ## Need to reclassify ?
    #Combine all forest, urban, agri, wetland...?
    #NLCD 1992, 2001,2006
    #20: urban , Developed (all 21,22,23,24)
    #40: forest, Forested upland (all 41,41,43)
    #80: agriculture, planted/cultivated (all,81,82)
    #90: wetland, (all,91,92,93...)

    ## Clean this up to make it more general and shorter...
        
    l_f_nlcd1992_subset_cat_dict = {
    "nlcd_1992_urban":(filter(lambda x: re.search(r'rec_2',x),outfile_breakout_list[0])),
    "nlcd_1992_forest":(filter(lambda x: re.search(r'rec_4',x),outfile_breakout_list[0])),
    "nlcd_1992_agriculture":(filter(lambda x: re.search(r'rec_8',x),outfile_breakout_list[0])),
    "nlcd_1992_wetland":(filter(lambda x: re.search(r'rec_9',x),outfile_breakout_list[0]))}
        
    l_f_nlcd2001_subset_cat_dict = {
    "nlcd_2001_urban":(filter(lambda x: re.search(r'rec_2',x),outfile_breakout_list[1])),
    "nlcd_2001_forest":(filter(lambda x: re.search(r'rec_4',x),outfile_breakout_list[1])),
    "nlcd_2001_agriculture":(filter(lambda x: re.search(r'rec_8',x),outfile_breakout_list[1])),
    "nlcd_2001_wetland":(filter(lambda x: re.search(r'rec_9',x),outfile_breakout_list[1]))}

    l_f_nlcd2006_subset_cat_dict = {
    "nlcd_2006_urban":(filter(lambda x: re.search(r'rec_2',x),outfile_breakout_list[2])),
    "nlcd_2006_forest":(filter(lambda x: re.search(r'rec_4',x),outfile_breakout_list[2])),
    "nlcd_2006_agriculture":(filter(lambda x: re.search(r'rec_8',x),outfile_breakout_list[2])),
    "nlcd_2006_wetland":(filter(lambda x: re.search(r'rec_9',x),outfile_breakout_list[2]))}

    l_f_nlcd_subset_cat_dict = {"nlcd92": l_f_nlcd1992_subset_cat_dict, "nlcd2001": l_f_nlcd2001_subset_cat_dict, "nlcd2006": l_f_nlcd2006_subset_cat_dict }
    nlcd_dss_cat = {}
    ##Combine categories in general cat...     
    for i in range(0,len(l_f_nlcd_subset_cat_dict)):
        l_f_nlcd = {}
        l_f_nlcd.update(l_f_nlcd_subset_cat_dict.values()[i])
         
        for j in range(0,len(l_f_nlcd)):
            var_name = l_f_nlcd.keys()[j]
            list_var_rast = l_f_nlcd.values()[j]
            out_file = var_name+"_" + out_suffix+file_format
            operation = "+"
            nlcd_dss_rast = raster_calc_operation_on_list(list_var_rast,out_dir,out_suffix,file_format,operation,NA_flag_val,out_file)
            #def raster_calc_operation_on_list(l_rast,out_dir,out_suffix_s,file_format,operation="+",NA_flag_val= -9999,out_file=Non
            #nlcd_dss_cat.update(nlcd_dss_rast) 
            #Add to dictionary
            nlcd_dss_cat[var_name] = nlcd_dss_rast
             
    ## Now change resolution for processing and display for DSS
    #generate new resolution...
    ## first at 3x30m for display (per land cover categories) use gdal wrap to calcuate the proportion of each land cover types at 90m
    ## second at 30x30m for calculations of proportions in POSTGIS
    lf_nlcd_prop900 = []
    for i in range(0,len(nlcd_dss_cat)):
        var_name = nlcd_dss_cat.keys()[i]
        infile = nlcd_dss_cat.values()[i]
        res_xy_val = [90,90] #this is for visualization...
        out_file = var_name+"_proportion_"+str(res_xy_val[0])+"_"+str(res_xy_val[1]) + "_" + out_suffix+file_format 
        resamp_opt = "average"
        output_type = "Float32"
       #[-ot {Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/
        change_resolution_raster(infile,res_xy_val,out_suffix,out_dir,file_format,output_type,NA_flag_val,out_file,resamp_opt)
        res_xy_val = [900,900] #this is for computation with postgis
        out_file = var_name+"_proportion_"+str(res_xy_val[0])+"_"+str(res_xy_val[1]) + "_" + out_suffix+file_format 
        change_resolution_raster(infile,res_xy_val,out_suffix,out_dir,file_format,output_type,NA_flag_val,out_file,resamp_opt)
        #change_resolution_raster(infile,res_xy_val,out_suffix_s,out_dir,file_format,output_type=None,NA_flag_val=-9999,out_file,resamp_opt="average"):
        lf_nlcd_prop900.append(out_file)

    #### Start of future function
    #Now Apply mask :line, make this a function...
     
    out_suffix_reg = "reg900_"+out_suffix #the region file projected in the defined projectio 
    w_extent_reg, reg_area_poly_projected = calculate_region_extent(shp_fname,out_suffix_reg,CRS_reg,out_dir)
    
    #first create mask from region definition:
    in_vect=reg_area_poly_projected
    in_rast=lf_nlcd_prop900[0]
    #out_suffix_s
    att_field="CNTYCODE"
    #file_format
    output_type="Float32"
    all_touched=True
    #NA_flag_val= -9999
    out_file=None
    EPSG_code="26919"
    #region_rast_fname = raster_to_poly_operation_on_list(in_vect,in_rast,out_suffix_s,att_field,file_format,output_type,all_touched,NA_flag_val,out_file)
    region_rast_fname = raster_to_poly_operation_on_list(in_vect,in_rast,out_suffix,att_field,file_format,output_type,all_touched,NA_flag_val,CRS_reg,out_file)

    # now either reclass or apply directly the mask
    #This functions creates an mask from an input raster.
    #ALl values that are not NA will be reclassified as 1 and all other as NA.
    #this function will be improve later on...  

    in_file = region_rast_fname
    out_file = "mask_regions_"+out_suffix_reg+file_format
    #mask_rast_file = create_raster_mask(in_file,ouout_dir,out_suffix_s,file_format,NA_flag_val,out_file)
    mask_rast_file = create_raster_mask(in_file,out_dir,out_suffix,file_format,CRS_reg,NA_flag_val,out_file)
    
    #problem with the mask does not contain full
    lf_var_masked = []
    lf_var = lf_nlcd_prop900
    for i in range(0,len(lf_var)):
        out_suffix_s = "_masked_"+out_suffix
        out_file = lf_var[i]
        in_file = lf_var[i]
        mask_file = mask_rast_file
        out_file = out_file.replace(out_suffix,out_suffix_s) #remove the suffix if it is there in the file name
        #f_masked = apply_raster_mask(in_file,mask_file,out_dir,out_suffix_s,file_format,NA_flag_val,out_file,EPSG_code)
        f_masked = apply_raster_mask(in_file,mask_file,out_dir,out_suffix_s,file_format,NA_flag_val,EPSG_code,out_file)
        #f_masked = pdb.runcall(apply_raster_mask,in_file,mask_file,out_dir,out_suffix_s,file_format,NA_flag_val,out_file,EPSG_code)
        lf_var_masked.append(f_masked)
        #Now get stat
    
    ## End of future function
    
    ###########################################################
    ##### Now process house density ...
    #dir_path = os.path.join(in_dir,"VarsFeb2014")
    lf_var_pop = glob.glob(os.path.join(in_dir,"se*pop.tif")) #US AEA, use none
    lf_var_hou = glob.glob(os.path.join(in_dir,"se*housing.tif")) #US AEA, use none
    lf_var_hist_temp = glob.glob(os.path.join(in_dir,"me_hist*.tif")) #Lat long EPSG 4326
    lf_var = lf_var_pop + lf_var_hou + lf_var_hist_temp

    #l_CRS_dst[j]
    l_CRS_dst=  [CRS_reg]*len(lf_var)
    l_CRS_src = ['+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs ']*(len(lf_var)-len(lf_var_hist_temp))
    l_CRS_src = l_CRS_src + [CRS_WGS84]*len(lf_var_hist_temp)

    #out_dir
    #in_dir 
    #file_format
    #out_suffix #set earlier
    w_extent_str=None #if null then use_reg_extent set to True
    NA_flag_val=None #[could be useful to have lists...]
    output_type=None
    clip_param=True #clip before reprojection
    reproject_param=True
    
   
    #Start function here ... 
        
    f_list = lf_var
    outfile_list = []
    
    #this can be paralellized later on...
    for j in range(0,len(f_list)):
        #d
        CRS_src = l_CRS_src[j] #source projection syst
        CRS_dst = l_CRS_src[j]  #target projection syst in this case CRS reg...

        out_suffix_dst = "_projected_"+out_suffix
        #CRS_dst = CRS_aea
        #shp_reg_outline = os.path.join("/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/"
        #                                       ,"county24.shp")
        #Combine with the other execution!!!
        if w_extent_str==None:
            use_reg_extent=True
        if w_extent_str!=None:
            use_reg_extent=False
            
        if use_reg_extent==True:
            w_extent, reg_area_poly_projected = calculate_region_extent(shp_reg_outline,out_suffix_dst,CRS_dst,out_dir)
        #w_extent= "-71.083923738042 47.4598539782516 -66.8854440488051 42.9171281482886"
        if use_reg_extent==False:
            w_extent = w_extent_str #this is in WGS84
        #end if

        #w_extent = w_extent_str #can vary by input?
        #outfile = pdb.runcall(create_raster_region,j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,out_suffix,out_dir,w_extent,clip_param=True)
        #outfile = create_raster_region(j,in_dir,infile_l_f,CRS_src,CRS_dst,file_format,out_suffix,out_dir,w_extent,clip_param=True)
        #output_type="Float64"
        #output_type = None   #does not work with real number  assignment...     
        CRS_dst = l_CRS_dst[j]  #target projection syst in this case CRS reg...
 
        outfile = create_raster_region(j,in_dir,f_list,CRS_dst,file_format,
                                  out_suffix,out_dir,w_extent,CRS_src,NA_flag_val=None,output_type=None,clip_param=True,reproject_param=True)

        outfile_list.append(outfile)
        
        
    return None
    
#Need to add this to run
if __name__ == '__main__':
    main()

#gdal_translate /ssi-dss-data/DSS_SSI_data/ser90_housing.tif projwin 1930763.951542  3013034.834472 2263786.175566  2487923.505723
#dst_dataset="test44.tif"
#src_dataset = "/ssi-dss-data/DSS_SSI_data/ser90_housing.tif"
#cmd_str = "".join(["gdal_translate",
#                              " "+"-projwin"+" "+w_extent,
#                              #" "+"-a_srs"+" '"+CRS_src+"'",
#                              #" -a_nodata "+NA_flag_val_str,
#                             " "+src_dataset, 
#                              " "+dst_dataset])            
#os.system(cmd_str)
#os.system("gdalsrsinfo -o proj4 /ssi-dss-data/DSS_SSI_data/ser90_housing.tif")
#'+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs '
#This is not matching the aea description assigned!!!

################# End of script ################# 