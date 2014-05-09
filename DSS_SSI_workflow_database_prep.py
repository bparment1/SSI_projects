#!/usr/bin/python
#
######## SUMMARIZE BY ZONES (POLYGON) WITH RASTER FROM POSTGIS  ########
#Script to process data for the DSS-SSI website.General purpose.
#This script computes average means for a vector polygon layer
#using a variable stored in raster format.
#Data inserted in a Postgis database and tables are created. 
##
## Authors: Benoit Parmentier 
# Created on: 04/02/2014
# Updated on: 05/03/2014

# TODO:
#  - functionalize to encapsulate high level procedural steps
### Need to add more functions by breaking out code!!...
## add function for use in multiprocessing pool

######## LOAD LIBRARY USED IN THE SCRIPT ###########

import os, glob, sys   #System tools: OS, files, env var etc.
import subprocess      #thread and processes
import re, zipfile     #Regular expression and zip tools
import datetime, calendar #Date processing
import ftplib             #Downloading ftp library
import argparse       #Agumnt for terminal callable scripts
import shutil        #Shell utilities
from osgeo import gdal       #Raster processing tools for geographic ata
from osgeo import ogr        #Vector processing tools for geographic data
from osgeo import osr        #Georeferencing: Spatial refreferen system for geographic data 
from osgeo import gdal_array #gdal geographic library
from osgeo import gdalconst  #gdal geographic library
import psycopg2              # Postgres binding, SQL query etc
import numpy as np           #Array, matrices and scientific computing
import pickle                #Object serialization
from multiprocessing import Process, Manager, Pool #parallel processing
import pdb                   #for debugging
import pandas as pd

################ NOW FUNCTIONS  ###################
#------------------
# Functions used in the script 
#------------------

    
def create_dir_and_check_existence(path):
    #Create a new directory
    try:
        os.makedirs(path)        
    except:
        print "directory already exists"
            
def save_data(data,fname):
    #save python object from memory to disk
    import pickle
    with open(fname, "wb") as f:
        pickle.dump(data, f)
        
def load_data(fname):
    #load python object from disk to memory
    import pickle
    try:
        with open(fname) as f:
            data_obj = pickle.load(f)
    except:
        print "No data"
        data_obj=[]
    return data_obj
            
def atoi(text):
    #convert str number to int in a string chain
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def create_summary_table(list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix):
    #This works but really needs to be improved!!
    #This function creates a table from the list wiht averages by
    #polygons. We used a numpy array rather than a numpy record array.
    #This only takes the mean right now so could store other information
    #later!!!
    
    table = numpy.ones((nb_rows,nb_columns))
    #add ID first...?
    for i in range(0,nb_columns):
        rows_output = list_rows[i]
        for j in range(0,nb_rows):
            tu = rows_output[j]
            table[j,i] = tu[element_nb]
    fname = "table_"+out_suffix+".txt"       
    outfile1 = os.path.join(out_dir,fname)
    numpy.savetxt(outfile1, table, fmt='%-7.6f')
    return table
    ## NEED TO ADD and ID Column!!!
    #May be add function to also store other values...        

def create_var_names_from_files(f_list):
    
    #assume the following string structure
    #tmax_1_clipped_projected_ncar_ccsm3_0_sres_a1b_2030s_tmax_2_5min_no_tile_asc04122014
    list_var_info_dict = []
    for j in range(0,len(f_list)):
        
        filename= os.path.basename(f_list[j])
        list_str= filename.split("_")
        month = list_str[1]
        decade = list_str[9][0:4]
        var = list_str[0]
        var_info_dict = {}
        var_info_dict["var"] = var
        var_info_dict["month"] = month
        var_info_dict["decade"] = decade
        
        list_var_info_dict.append(var_info_dict)
    
    return list_var_info_dict
    
def shp_to_data_frame(shp_fname):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds =driver.Open(shp_fname,0)
    if ds is None:
        print "Could not open file"
        sys.exit(1)
    lyr= ds.GetLayer(0)    
    #feature= lyr.GetFeature(0)
    #Now get the attributes...
    #make this a function
    layerDefinition = lyr.GetLayerDefn()
    column_names = []
    #ncolumn = layerDefinition.GetFieldCount()
    for i in range(layerDefinition.GetFieldCount()):       
        column_names.append(layerDefinition.GetFieldDefn(i).GetName())
    #table = {}
    
    #list_values = []
    dict_values = {} #this will store the column values...
    for column in column_names:
        name=column
        value_att = []
        #value_att= {}
        for i in range(lyr.GetFeatureCount()):
            #f is a featuur 
            f =lyr.GetFeature(i)
            val = f.GetField(name)
            value_att.append(val)
            f.Destroy()
            
            #value_att.append(f.GetField(name))
            #value_att.update(f.GetField(name))
            #f.Destroy()
        #Maybe should make val_att a data.frame here or a series... 
        #pd.Series()
        d = {name:value_att}
        #list_values.append(value_att)
        dict_values.update(d)
    #now make it a data.frame
    df_table = pd.DataFrame(dict_values,columns=column_names)
    #df_table = pd.DataFrame(list_values,columns=column_names)
    return df_table

    #print "Name  -  Type  Width  Precision"
    

def add_field_to_shp(in_fname,out_fname,dict_col=None):
    #This functions add a field to a shape file. It can also be used to make a copy of the shapefil.
    #Add option for reprojection!!!
    #Also clean so it works for any vector type (ie other than ESRI Shapefile)
    
    ### START SCRIPT ####
    
    OGRTypes = {int: ogr.OFTInteger, str: ogr.OFTString, float: ogr.OFTReal} #note there are missing types!!

    #Set driver
    DriverName = "ESRI Shapefile"      # e.g.: GeoJSON, ESRI Shapefile
     
    #Read in data source 
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds =driver.Open(in_fname,1) #writeable for one
    if ds is None:
        print "Could not open file"
        sys.exit(1)
        
    #Get the input layer    
    lyr = ds.GetLayer(0)    
     
    #Create destination data file and source 
    driver = ogr.GetDriverByName(DriverName)
    if os.path.exists(out_fname):
        driver.DeleteDataSource(out_fname)

    dst_ds = driver.CreateDataSource( out_fname )
    
    #Create output laer
    
    #proj = osr.SpatialReference()  
    #proj.SetWellKnownGeogCS( "EPSG:4326" )  
    proj =lyr.GetSpatialRef()
    geometry_type = lyr.GetLayerDefn().GetGeomType()    
    #outLayer = outDataSource.CreateLayer("states_convexhull", geom_type=geometry_type)
    dst_layer = dst_ds.CreateLayer('region',
                            srs = proj,
                            geom_type=geometry_type)
    inFeature  = lyr.GetFeature(0)
    
    
    # add fields: use input layer definitions
    inLayerDefn = lyr.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        dst_layer.CreateField(fieldDefn)

    # get the output layer's feature definition
    dst_layer_Defn = dst_layer.GetLayerDefn()

    # loop through the input features
    #inFeature = inLayer.GetNextFeature()

    #Now copy the feature in the new layer
    #reproject if True 

    for j in range(lyr.GetFeatureCount()):
        #f is a featuur 
        # get the input geometry
        inFeature = lyr.GetFeature(j) #getting feature i    
        geom = inFeature.GetGeometryRef() #get geometry
        # reproject the geometry
        #geom.Transform(coordTrans)
        # create a lynew feature
        feature_new = ogr.Feature(dst_layer_Defn)
        # set the geometry and attribute
        feature_new.SetGeometry(geom)
        for i in range(0, dst_layer_Defn.GetFieldCount()):
            feature_new.SetField(dst_layer_Defn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))

            
        # add the feature to the shapefile
        dst_layer.CreateFeature(feature_new)
        # destroy the features and get the next input feature
        feature_new.Destroy()
        inFeature.Destroy()

    #CLOSE SHAPEFILES...
    dst_ds.Destroy()
    ds.Destroy()


    #If new field is True:
    if dict_col is None:
        print "No field added "
    #    dst_ds.Destroy()
    #    ds.Destroy()
        sys.exit(1)

    #Read in data source 
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds =driver.Open(out_fname,1) #writeable for one
    if ds is None:
        print "Could not open file"
        sys.exit(1)
    dst_layer = dst_ds.GetLayer(0)    
    # get the output layer's feature definition
    dst_layer_Defn = dst_layer.GetLayerDefn()

    field_name = dict_col.keys()[0]
    data_type =  type(dict_col.values()[0][0]) #get python data type from the first item in the dict 
    new_field =  ogr.FieldDefn(field_name, OGRTypes[data_type])
    #fldDef = ogr.FieldDefn('Name', ogr.OFTString)
    #fldDef2.SetWidth(16) #16 char string width
    dst_layer.CreateField(new_field)
    
    #inFeature = inLayer.GetNextFeature()
    #while inFeature:

    for j in range(dst_layer.GetFeatureCount()):
        #f is a featuur 
            # get the input geometry
        inFeature = dst_layer.GetFeature(j) #getting feature i    
        val = dict_col.values()[0][j]

        inFeature.SetField(field_name,val) # assign value to feature for new field
        #dst_layer.SetFeature(feat) 
        # add the feature to the shapefile
        #dst_layer.CreateFeature(inFeature)
        dst_layer.SetFeature(inFeature)  #NEED TO USE SET FEATURE NOT CREATE FEATURE!!!
        
        # destroy the features and get the next input feature
        #feature_new.Destroy()
        inFeature.Destroy()

    #CLOSE SHAPEFILES...
    dst_ds.Destroy()
    #ds.Destroy()
    #END OF FUNCTION

def get_vct_FieldName(shp_fname):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds =driver.Open(shp_fname,0)
    if ds is None:
        print "Could not open file"
        sys.exit(1)
    lyr= ds.GetLayer(0)    
    #make this a function
    layerDefinition = lyr.GetLayerDefn()
    list_names = []
    #print "Name  -  Type  Width  Precision"
    for i in range(layerDefinition.GetFieldCount()):  
        fieldName =  layerDefinition.GetFieldDefn(i).GetName()
        fieldTypeCode = layerDefinition.GetFieldDefn(i).GetType()
        fieldType = layerDefinition.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)
        fieldWidth = layerDefinition.GetFieldDefn(i).GetWidth()
        GetPrecision = layerDefinition.GetFieldDefn(i).GetPrecision()
        info = [fieldName,fieldType,str(fieldWidth),GetPrecision]
        list_names.append(info)
        #print fieldName + " - " + fieldType+ " " + str(fieldWidth) + " " + str(GetPrecision)
    nf = lyr.GetFeatureCount() #number of features
    return list_names, nf

def caculate_zonal_statistics(i,out_dir,out_suffix,rast_fname,shp_fname,SRS_EPSG,region_id_col,postgres_path_bin,db_name,user_name):
#def caculate_zonal_statistics(i,out_dir,out_suffix,rast_fname,shp_fname,SRS_EPSG,postgres_path_bin,db_name,user_name):
    
    #This function calculate summary statistic i.e. average per region in the
    #input shapefile.
    #Raster and polygon files are imported in postgis where caluclation takes place.

    #Write a function to excute the next lines in a loop (for paralelization)
    #Multiprocessing module requires individual db connections for each process 
    
    #INPUT PARAMETERS
    #i
    #l_f
    #out_dir
    #out_suffix
    #rast_fname
    #shp_fname
    #in_dir
    #SRS_EPSG
    #postgres_path_bin    
    #db_name ="test_ssi_db2"
    #user_name = "benoit"
    
    #polygon_input_shp_file = "county24.shp"
    #polygon_input_shp_file = shp_fname
    #polygon_input_shp_file = "metwp24.shp"
    #out_suffix = "03242014"
    #os.chdir(in_dir) #set current dir
    #os.getcwd() #check current dir
    #out_dir = in_dir
    #out_path<-"/data/project/layers/commons/data_workflow/output_data"
    #out_dir = in_dir
    #SRS_EPSG = "2037"             
    #Database information       
    #db_name ="test_ssi_db2"
    #user_name = "benoit"
    #user_name = "parmentier" 
    #postgres_path_bin = "/usr/lib/postgresql/9.1/bin/"
    #postgres_path_bin = ""
    
    polygon_input_shp_table = "counties"    #can be any name
                                            # polygons can be any regions
    #polygon_input_shp_table = "towns"
    
    rast_input = rast_fname[i]
    #output_table_rast = "mtemp"
    output_table_rast ="".join(["mtemp_",str(i)])
    poly_table = "".join(["mht_poly_use_",str(i)]) #This is the name of the table that will store the values of 
    #output_table_rast ="mtemp"
    #poly_table = "mht_poly_use" #This is the name of the table that will store the values of 
    
    #polygon_input_shp_file
    #shp_input = os.path.join(in_dir,polygon_input_shp_file)
    shp_input = shp_fname
    
    output_table_shp = "".join([polygon_input_shp_table,str(i)])
    #output_table_shp = polygon_input_shp_table

    ####
    ##OPEN A CONNECTION IN THE DATABASE
    
    ##Make this a function
    sys.path.append(postgres_path_bin)
    #sys.path.append(postgres_path_bin)
    try:
        cmd_db_str = "".join(["dbname=",db_name," ","user=",user_name])
        conn= psycopg2.connect(cmd_db_str)
    except:
        print "I am unable to connect to the database"
    
    conn.set_isolation_level(0) #change isolation level,so that we are not in transaction mode
    cur = conn.cursor()
    #end of function
    
    ## IMPORT SHAPEFILE IN TABLE
    SQL_str = "DROP TABLE IF EXISTS %s;" % (output_table_shp)
    cur.execute(SQL_str) #Should collect all commands executed in a file (for later)
    
    #Make this a function
    
    cmd_str = "".join([postgres_path_bin,
                      "shp2pgsql",
                      " ","-s ",SRS_EPSG,             #projection system add as input
                      " ","-c -I",               #additional options
                      " ",shp_input,
                      " ",output_table_shp,
                      " ","| psql",
                      " ","-U ",user_name,
                      " ","-d ",db_name])
    os.system(cmd_str)
    
    ## IMPORT RASTER IN POSTGIS
    SQL_str = "DROP TABLE IF EXISTS %s;" % (output_table_rast)
    cur.execute(SQL_str) #Should collect all commands executed in a file (for later)
    
    cmd_str = "".join([postgres_path_bin,
                      "raster2pgsql",
                      " ","-s ",SRS_EPSG,             #projection system add as input
                      " ","-I -t 1x1",               #additional options
                      " ",rast_input,
                      " ",output_table_rast,
                      " ","| psql",
                      " ","-U ",user_name,
                      " ","-d ",db_name,
                      " >rast.log"])
    os.system(cmd_str)
    
    #First convert raster image into a polygon...
    
    #cur.execute("DROP TABLE IF EXISTS mht_poly_use");
    SQL_str = "DROP TABLE IF EXISTS %s;" % (poly_table)
    cur.execute(SQL_str) #Should collect all commands executed in a file (for later)
    
    #Note that rast is column that contains values of tiles imported from the raster
    SQL_str = "SELECT ST_ConvexHull(rast) AS pixelgeom, ST_Area(ST_ConvexHull(rast)) AS pixelarea,(ST_SummaryStats(rast)).sum AS pixelval INTO %s FROM %s WHERE (ST_SummaryStats(rast)).sum IS NOT NULL;" % (poly_table,output_table_rast)
    cur.execute(SQL_str) 
    #NOW create an index
    #cur.execute("CREATE INDEX idx_pixelgeom ON mht_poly_use USING GIST(pixelgeom);")
    idx_pixelgeom = "idx_pixelgeom_%s" % (str(i))
    SQL_str = "CREATE INDEX %s ON %s USING GIST(pixelgeom);" % (idx_pixelgeom,poly_table)
    cur.execute(SQL_str)
    #Change cntycode to Name...also save the type column?--> make it a table...
    #region_id_col
    #SQL_str ="WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, cntycode FROM %s, %s WHERE ST_Intersects(pixelgeom,geom)) SELECT cntycode, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by cntycode;" % (poly_table,output_table_shp)
    #SQL_str ="CREATE TABLE zonal_stat_regions AS WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, cntycode FROM %s, %s WHERE ST_Intersects(pixelgeom,geom)) SELECT cntycode, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by cntycode;" % (poly_table,output_table_shp)
    SQL_str ="CREATE TABLE zonal_stat_regions AS WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, %s FROM %s, %s WHERE ST_Intersects(pixelgeom,geom)) SELECT %s, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by %s;" % (region_id_col,poly_table,output_table_shp,region_id_col,region_id_col)

    #cmd_str ="WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, cntycode FROM mht_poly_use, counties WHERE ST_Intersects(pixelgeom,geom)) SELECT cntycode, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by cntycode;"
    cur.execute(SQL_str)       
    SQL_str = "SELECT * FROM zonal_stat_regions"
    cur.execute(SQL_str)       

    rows_table = cur.fetchall()
    
    #list_rows.append(rows)
    #writeout list....    
    fname = "rows"+str(i)+"_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(rows_table,fname)
    
    #Close connection with postgres
    cur.close()
    conn.close()
    
    return rows_table
    
    
#######################################################################
######################### BEGIN SCRIPT  ###############################
#--------------------------------------------
# Script run by arguments from the shell?
#--------------------------------------------
            
def main():
    #
    #--------------------------------------------
    #Compute summary by zones/regions using POSTGIS
    #--------------------------------------------
 
    ########## READ AND PARSE PARAMETERS AND ARGUMENTS #########
        
    #in_dir_poly =     #in_dir = "/ssi-dss-data/DSS_SSI_data/"
    #in_dir ="/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/"
    in_dir ="/ssi-dss-data/DSS_SSI_data"
    
    shp_fname = os.path.join(in_dir,"county24.shp")
    region_name = "COUNTY" #name of the column containing the id for each entities...
    #if region_name is None then create an ID? Add this option later on.
    region_type = "C"    #type for the region entity: C for county, T for town  
    #valueType = temperature

    #polygon_input_shp_file = "metwp24.shp"

    #EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
    #+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
    out_suffix = "05032014"
    os.chdir(in_dir) #set current dir
    #os.getcwd() #check current dir
    out_dir = in_dir

    #out_path<-"/data/project/layers/commons/data_workflow/output_data"
    out_dir = "output_data_"+out_suffix
    out_dir = os.path.join(in_dir,out_dir)
    create_dir_and_check_existence(out_dir)
    
    os.chdir(out_dir) #set working directory        

    #EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
    #+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
    #file_format = ".tif"

    SRS_EPSG = "2037"             
    #Database information       
    db_name ="test_ssi_db3"
    user_name = "benoit"
    #user_name = "parmentier" 
    postgres_path_bin = "/usr/lib/postgresql/9.1/bin/"  #on SSI Maine
    #postgres_path_bin = ""  #on ATLAS NCEAS
    
    #polygon_input_shp_table = "counties"    
    #polygon_input_shp_table = "towns"
    
    fileglob_pattern = "*projected*ncar*.tif"
    pathglob = os.path.join(out_dir, fileglob_pattern)
    l_f_temp = glob.glob(pathglob)
    l_f_temp.sort(key=natural_keys) #mixed sorting of files

    ########## START SCRIPT #############
         
    
    ### FIRST ADD A COLUMN NAME, VALUE TYPE
    
    #get counties names from shapefile,type=c,Time=year,valueType=tmin_1,tmax_1 etc., value
    #test = shp_to_data_frame(shp_fname)
    ### NOW REORGANIZE THE DATA INTO...
    #Name,type,Time,valueType,Value    
    #../county24.shp
    
    layer_name = os.path.splitext(os.path.basename(shp_fname))[0]
    #os.system("ogrinfo "+shp_fname+" -sql \"SELECT COUNT(*) FROM "+layer_name+"\"")
    nf = get_vct_FieldName(shp_fname)[1] #extract number of feature using user defined function (see above)
    region_type_list = list(region_type)*nf
    dict_region_type = {"Type":region_type_list}
    
    df_region = shp_to_data_frame(shp_fname)
    region_name_list = list(df_region[region_name])
    dict_region_name = {"Name":region_name_list}
    #add_field_to_shp(in_fname,out_fname,dict_col=None):
    out_fname1 = layer_name+"_out_"+out_suffix+".shp"
    add_field_to_shp(shp_fname,out_fname1,dict_region_name)
    #df_region1 = shp_to_data_frame(out_fname1)
    out_fname2 = layer_name+"_out_"+out_suffix+".shp"
    add_field_to_shp(out_fname1,out_fname2,dict_region_type)
    df_region2 = shp_to_data_frame(out_fname1)
                 
    ## Temp processing: 2020s,2030s,2040s,2050s for tmin, tmax and tmean (e.g. 4*3 directories with 12 files...)

    #1.download...not automated...
    #2.unzip asc or grd (arc grid)
    #3.convert to tif (gdal_tranlslate)
    #4.reproject and clip/subset for Maine region (clip using count24?)

    #Get tmin files
    l_f_tmin = filter(lambda x: re.search(r'tmin',x),l_f_temp)
    #Get tmax files
    l_f_tmax = filter(lambda x: re.search(r'tmax',x),l_f_temp)
    #Get tmean files
    l_f_tmean = filter(lambda x: re.search(r'tmean',x),l_f_temp)
    
    dict_rast_fname = {"tmin": l_f_tmin, "tmax": l_f_tmax, "tmean": l_f_tmean }
    #debug
    #test = pdb.runcall(caculate_zonal_statistics,i,out_dir,out_suffix,rast_fname,shp_fname,SRS_EPSG,postgres_path_bin,db_name,user_name)
    
    ### IMPORTANT BEFORE RUNNING THE QUERY WE NEED TO MODIFY THE SHP TO HAVE A ID FIELD, NAME FIELD AND TYPE FIELD
    
    
    #This will be parallelized and looped through dict_rast_fname
    #can make one additional loop to reduce the repitition with tmin, tmax and tmean

    list_rows_tmin = [] # defined lenth right now    
    rast_fname = l_f_tmin #copy by refernce!!

    for i in range(0,len(rast_fname)):
    #for i in range(0,2):    
        out_suffix_s = "tmin_"+out_suffix
        rows = caculate_zonal_statistics(i,out_dir,out_suffix_s,rast_fname,shp_fname,SRS_EPSG,postgres_path_bin,db_name,user_name)
        list_rows_tmin.append(rows) #add object rows to list

    list_rows_tmax = [] # defined lenth right now
    rast_fname = l_f_tmax #copy by refernce!!
    
    for i in range(0,len(rast_fname)):
    #for i in range(0,2):      
        out_suffix_s = "tmax_"+out_suffix
        rows = caculate_zonal_statistics(i,out_dir,out_suffix_s,rast_fname,shp_fname,SRS_EPSG,postgres_path_bin,db_name,user_name)
        list_rows_tmax.append(rows) #add object rows to list
        
    list_rows_tmean = [] # defined lenth right now
    rast_fname = l_f_tmean #copy by refernce!!
    
    for i in range(0,len(rast_fname)):
    #for i in range(0,2):     
        out_suffix_s = "tmean_"+out_suffix
        rows = caculate_zonal_statistics(i,out_dir,out_suffix_s,rast_fname,shp_fname,SRS_EPSG,postgres_path_bin,db_name,user_name)
        list_rows_tmean.append(rows) #add object rows to list
        #Note the first columns from row contains: region ID
        #Other columns contain in hte order: region ID, sums, means, counts, maxes, mins, area

    #Write out results by combining
    #test1 =load_data(fname)
    fname = "list_rows_min_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_rows_tmin,fname)
    
    nb_columns = len(list_rows_tmin) # + 1
    nb_rows = len(list_rows_tmin[0])
    #table_test = pdb.runcall(create_summary_table,list_rows_tmin,nb_rows,nb_columns, out_dir,out_suffix)
    #table_tmin = create_summary_table(list_rows_tmin,nb_rows,nb_columns, out_dir,out_suffix)
    #create_summary_table(list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix)

    ##Write out tmax
    fname = "list_rows_max_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_rows_tmax,fname)
    
    list_rows_tmin = load_data("list_rows_min_05032014.dat")
    #list_rows_summary = [list_rows_tmin,list_rows_tmax,list_rows_tmean]
    dict_rows_summary = {"tmin":list_rows_tmin,"tmax":list_rows_tmax,"tmean":list_rows_tmean}
    dict_table = {}
    
    for i in range(0,len(dict_rows_summary)):
        list_rows = dict_rows_summary.values()[i]
        var_name = dict_rows_summary.keys()[i]
        nb_columns = len(list_rows) # + 1
        nb_rows = len(list_rows[0])
        element_nb = 2 #this is the column containing the mean for each entity (third column from postgis table)    
        #table_test = pdb.runcall(create_summary_table,list_rows_tmin,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        table = create_summary_table(list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        dict_table[var_name] = table
        #Now add id back...
        nb_columns = 1
        element_nb = 0 #this is the column containing the mean for each entity (third column from postgis table)    
        id_reg = create_summary_table(list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        region_id = ["%.0d" % number for number in id_reg]
        #"".join([["f_"]*len(region_id),region_id])    
        #NOW FORMAT TALBE TO PUT IN POSTGIS!!!
        l_f = dict_rast_fname.values()[i]       
        var_info = create_var_names_from_files(l_f)
        df_var = pd.DataFrame(np.transpose(table))
        df_var.columns = region_id
        df_info = pd.DataFrame(var_info)
        df_var["decade"] = df_info["decade"]
        df_var["month"] = df_info["month"]
        df_var["var"] = df_info["var"]
        tmp_name = ["decade","month","var"] + region_id
        df = df_var[tmp_name] #reorder the columns of the dataframe in python    
        df.to_csv("table_"+var_name+"_df_"+out_suffix+".csv",sep=",") #write out table
        
        #Now reshape to have a 48*16 rows and 4 columns... the last should be ID region
        nb_columns = len(list_rows) # + 1
        tt=table.ravel()
        df_val = pd.DataFrame(tt,columns=["value"])
        df_val["time"] = (list(df["decade"]))*nb_rows
        df_val["valueType"] = (list(df_var["var"]+"_"+df_var["month"]))*nb_rows
        df_val["type"]= (list("c"))*df.value.count() #using the number of count in column value
        del table
        
    #Write function to add back the tables in postgis and join them?   

    return None
    
#Need to add this to run
if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
#
################   END OF SCRIPT  ###############