#!/usr/bin/python
#
######## SUMMARIZE BY ZONES (POLYGON) WITH RASTER FROM POSTGIS  ########
#
#Script to process data for the DSS-SSI website.
#It is aimed at general purpose (i.e. script must strives for generality and automation)
#This script computes average means for a vector polygon layer in postgis.
#using a variable stored in raster format.
#Data inserted in a Postgis database and tables are created. 
##
## Authors: Benoit Parmentier 
# Created on: 04/02/2014
# Updated on: 05/27/2014
# Project: DSS-SSI
#
# TODO:
# - Need to add more functions by breaking out code!!...
# - improve performance...by using in multiprocessing pool
#
######## LOAD LIBRARY/MODULES USED IN THE SCRIPT ###########

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
import pandas as pd          #DataFrame object and other R like features for data munging

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

def create_summary_table(list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix,use_str=False):
    #This works but really needs to be improved!!
    #This function creates a table from the list wiht averages by
    #polygons. We used a numpy array rather than a numpy record array.
    #This only takes the mean right now so could store other information
    #later!!!
    import numpy as np           #Array, matrices and scientific computing

    table = np.ones((nb_rows,nb_columns))
    if use_str==True:
        table = np.array(table,np.dtype("a64")) #convert to string if needed
    
    #add ID first...?
    for i in range(0,nb_columns):
        rows_output = list_rows[i]
        for j in range(0,nb_rows):
            tu = rows_output[j]
            table[j,i] = tu[element_nb]
    fname = "table_"+out_suffix+".txt"       
    outfile1 = os.path.join(out_dir,fname)
    if use_str==False:
        np.savetxt(outfile1, table, fmt='%-7.6f')
    if use_str==True:
        np.savetxt(outfile1, table, fmt="%s")
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
 
## More general function for both NCAR and NLCD data? 
def create_var_names_from_files_NLCD(f_list):
    
    #assume the following string structure
    #tmax_1_clipped_projected_ncar_ccsm3_0_sres_a1b_2030s_tmax_2_5min_no_tile_asc04122014
    list_var_info_dict = []
    for j in range(0,len(f_list)):
        
        filename= os.path.basename(f_list[j])
        list_str= filename.split("_")
        var = list_str[0]
        year = list_str[1]
        category = list_str[2]
        var_info_dict = {}
        var_info_dict["var"] = var
        var_info_dict["year"] = year
        var_info_dict["category"] = category
        
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

def caculate_zonal_statistics(i,out_dir,out_suffix,rast_fname,shp_fname,SRS_EPSG,region_id_col,postgres_path_bin,db_name,user_name,NA_flag_val,tile_size=1):
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
    
    NA_flag_val_str = str(NA_flag_val)
    polygon_input_shp_table = "counties"    #can be any name
                                            # polygons can be any regions
    #polygon_input_shp_table = "towns"
    
    rast_input = rast_fname[i]
    #rast_input = "tmax_7_clipped_projected_ncar_ccsm3_0_sres_a1b_2020s_tmax_2_5min_no_tile_asc05152014.tif"
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
    #conn.set_isolation_level(0) #change isolation level,so that we are not in transaction mode

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
        
    rast_tiling = str(tile_size)+"x"+str(tile_size)
    cmd_str = "".join([postgres_path_bin,
                      "raster2pgsql",                 #We can potentially use ...
                      " ","-s ",SRS_EPSG,             #projection system add as input
                      " ","-I ",                      #build spatial index, overview factor
                      #" ","-N ", NA_flag_val_str,     #Nod data val to use on bands without a NODATA value...
                      " ","-t ",rast_tiling,          #store raster in nxn tile table
                      " ",rast_input,                 #input raster name
                      " ",output_table_rast,
                      " ","| psql",
                      " ","-U ",user_name,
                      " ","-d ",db_name,
                      " >rast.log"])
    os.system(cmd_str)
    
    #First convert raster image into a polygon...
    
    #cur.execute("DROP TABLE IF EXISTS mht_poly_use_0");
    SQL_str = "DROP TABLE IF EXISTS %s;" % (poly_table)
    cur.execute(SQL_str) #Should collect all commands executed in a file (for later)
    
    #Note that rast is column that contains values of tiles imported from the raster
    #rid is the tile id,rast is the column containing the raster object
    #http://postgis.net/docs/RT_ST_SummaryStats.html
    #This is where we need to add the number of pixels per tiles and number of no_data values
    #This will allow weighting of the pixel value (which is the sum...)
    SQL_str = "SELECT ST_ConvexHull(rast) AS pixelgeom, ST_Area(ST_ConvexHull(rast)) AS pixelarea,(ST_SummaryStats(rast)).sum AS pixelval INTO %s FROM %s WHERE (ST_SummaryStats(rast)).sum IS NOT NULL;" % (poly_table,output_table_rast)
    #SQL_str = "SELECT ST_ConvexHull(rast) AS pixelgeom, ST_Area(ST_ConvexHull(rast)) AS pixelarea,(ST_SummaryStats(rast)).sum AS pixelval INTO %s FROM %s WHERE (ST_SummaryStats(rast)).sum IS NOT NULL;" % (poly_table,output_table_rast)
    cur.execute(SQL_str) 
    #NOW create an index
    #cur.execute("CREATE INDEX idx_pixelgeom ON mht_poly_use USING GIST(pixelgeom);")
    
    idx_pixelgeom = "idx_pixelgeom_%s" % (str(i))
    SQL_str = "DROP INDEX IF EXISTS %s;" % (idx_pixelgeom)
    cur.execute(SQL_str) 

    SQL_str = "CREATE INDEX %s ON %s USING GIST(pixelgeom);" % (idx_pixelgeom,poly_table)
    cur.execute(SQL_str)
    #Change cntycode to Name...also save the type column?--> make it a table...
    #region_id_col
    #SQL_str ="WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, cntycode FROM %s, %s WHERE ST_Intersects(pixelgeom,geom)) SELECT cntycode, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by cntycode;" % (poly_table,output_table_shp)
    #SQL_str ="CREATE TABLE zonal_stat_regions AS WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, cntycode FROM %s, %s WHERE ST_Intersects(pixelgeom,geom)) SELECT cntycode, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by cntycode;" % (poly_table,output_table_shp)
    
    zonal_stat_regions = "zonal_stat_regions_%s" % (str(i))
    SQL_str = "DROP TABLE IF EXISTS %s;" % (zonal_stat_regions)
    cur.execute(SQL_str)  
    
    SQL_str = "CREATE TABLE %s AS WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, %s FROM %s, %s WHERE ST_Intersects(pixelgeom,geom)) SELECT %s, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by %s;" % (zonal_stat_regions,region_id_col,poly_table,output_table_shp,region_id_col,region_id_col)
    #cmd_str ="WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, cntycode FROM mht_poly_use, counties WHERE ST_Intersects(pixelgeom,geom)) SELECT cntycode, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by cntycode;"
    cur.execute(SQL_str)       
    
    SQL_str = "SELECT * FROM %s" % (zonal_stat_regions)
    cur.execute(SQL_str)       

    rows_table = cur.fetchall()
    
    #Also export in panda data.frame format if it is available on the system (add default option to function!!!)
    #cur.execute("select instrument, price, date from my_prices")
    #df = DataFrame(cur.fetchall(), columns=['instrument', 'price', 'date'])
    #df.set_index('date', drop=False)
    #import pandas.io.sql as sql
    #sql.read_frame("select * from test",con)
    
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
    #region_name = "CNTYCODE"
    #if region_name is None then create an ID? Add this option later on.
    region_type = "C"    #type for the region entity: C for county, T for town  
    valueType = ["tmin","tmax","tmean"] #temperature the list can be one only...
    zonal_stat = "mean" #This is the statistic extracted from the region

    #polygon_input_shp_file = "metwp24.shp"

    #EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
    #+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
    out_suffix = "05152014"
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
    db_name ="test_ssi_db2"
    
    user_name = "benoit"
    #user_name = "parmentier" 
    postgres_path_bin = "/usr/lib/postgresql/9.1/bin/"  #on SSI Maine
    tile_size = 10  #this set the tile size for the raster postgis table
    NA_flag_val = -9999
    #postgres_path_bin = ""  #on ATLAS NCEAS
    
    #polygon_input_shp_table = "counties"    
    #polygon_input_shp_table = "towns"
    
    ## Get input raster lists containg information to be summarized
    #Temperature info:
    fileglob_pattern = "*projected*ncar*.tif"
    pathglob = os.path.join(out_dir, fileglob_pattern)
    l_f_valueType = glob.glob(pathglob) #this contains the raster variable files that need to be summarized
    l_f_valueType.sort(key=natural_keys) #mixed sorting of files
    #NLCD land cover info:
    fileglob_pattern = "nlcd_*_*_proportion_900_900_*.tif"
    pathglob = os.path.join(out_dir, fileglob_pattern)
    l_f_valueType_NLCD = glob.glob(pathglob) #this contains the raster variable files that need to be summarized
    l_f_valueType_NLCD.sort(key=natural_keys) #mixed sorting of files
       

    ########## START SCRIPT #############
            
    ### FIRST ADD A COLUMN NAME, VALUE TYPE
    
    #get counties names from shapefile,type=c,Time=year,valueType=tmin_1,tmax_1 etc., value
    
    layer_name = os.path.splitext(os.path.basename(shp_fname))[0]
    #os.system("ogrinfo "+shp_fname+" -sql \"SELECT COUNT(*) FROM "+layer_name+"\"")
    nf = get_vct_FieldName(shp_fname)[1] #extract number of feature using user defined function (see above)
    region_type_list = list(region_type)*nf
    dict_region_type = {"Type":region_type_list}
    
    df_region = shp_to_data_frame(shp_fname)
    region_name_list = list(df_region[region_name])
    dict_region_name = {"Name":region_name_list}
    #pdb.runcall(add_field_to_shp,shp_fname,out_fname1)
    out_fname1 = layer_name+"_out_"+out_suffix+".shp"
    add_field_to_shp(shp_fname,out_fname1,dict_region_name)
    df_region1 = shp_to_data_frame(out_fname1)
    #out_fname2 = layer_name+"_out_"+out_suffix+".shp"
    #add_field_to_shp(out_fname1,out_fname2,dict_region_type)
    #df_region2 = shp_to_data_frame(out_fname1)
                            
    ## Temp processing: 2020s,2030s,2040s,2050s for tmin, tmax and tmean (e.g. 4*3 directories with 12 files...)

    #Get raster files that were process previously, sort files by variable...
    dict_rast_fname = {}
    l_f_tmin = filter(lambda x: re.search(r'tmin',x),l_f_valueType)        
    l_f_tmax = filter(lambda x: re.search(r'tmax',x),l_f_valueType)
    l_f_tmean = filter(lambda x: re.search(r'tmean',x),l_f_valueType)        
    dict_rast_fname = {"tmin": l_f_tmin, "tmax": l_f_tmax, "tmean": l_f_tmean }
        
    #This will be parallelized and looped through dict_rast_fname 
    #can make one additional loop to reduce the repitition with tmin, tmax and tmean
    
    list_rows_NLCD = [] # defined lenth right now    
    rast_fname = l_f_valueType_NLCD #copy by refernce!!

    for i in range(1,len(rast_fname)):
    #for i in range(0,2):    
        out_suffix_s = "NLCD_"+out_suffix
        #tile_size = 1
        rows = caculate_zonal_statistics(i,out_dir,out_suffix_s,rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name,NA_flag_val,tile_size)
        list_rows_NLCD.append(rows) #add object rows to list

    list_rows_tmin = [] # defined lenth right now  f  
    rast_fname = l_f_tmin #copy by refernce!!

    for i in range(0,len(rast_fname)):
    #for i in range(0,2):    
        out_suffix_s = "tmin_"+out_suffix
        #tile_size = 1
        rows = caculate_zonal_statistics(i,out_dir,out_suffix_s,rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name,NA_flag_val,tile_size)
        list_rows_tmin.append(rows) #add object rows to list
        #test = pdb.runcall(caculate_zonal_statistics,i,out_dir,out_suffix_s,rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name)

    list_rows_tmax = [] # defined lenth right now
    rast_fname = l_f_tmax #copy by refernce!!
    
    for i in range(0,len(rast_fname)):
    #for i in range(0,2):      
        out_suffix_s = "tmax_"+out_suffix
        rows = caculate_zonal_statistics(i,out_dir,out_suffix_s,rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name)
        #rows = caculate_zonal_statistics(i,out_dir,out_suffix_s,rast_fname,shp_fname,SRS_EPSG,region_id_col,postgres_path_bin,db_name,user_name)
        list_rows_tmax.append(rows) #add object rows to list
        
    list_rows_tmean = [] # defined lenth right now
    rast_fname = l_f_tmean #copy by refernce!!
    
    for i in range(0,len(rast_fname)):
    #for i in range(0,2):     
        out_suffix_s = "tmean_"+out_suffix
        #test = pdb.runcall(caculate_zonal_statistics,i,out_dir,out_suffix_s,rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name)
        rows = caculate_zonal_statistics(i,out_dir,out_suffix_s,rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name)
        #rows = caculate_zonal_statistics(i,out_dir,out_suffix_s,rast_fname,shp_fname,SRS_EPSG,postgres_path_bin,db_name,user_name)
        list_rows_tmean.append(rows) #add object rows to list
        #Note the first columns from row contains: region ID
        #Other columns contain in hte order: region ID, sums, means, counts, maxes, mins, area

    nb_columns = len(list_rows_tmin) # + 1
    nb_rows = len(list_rows_tmin[0])
     
    #Write out results by combining
    #test1 =load_data(fname)
    fname = "list_rows_min_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_rows_tmin,fname)
    ##Write out tmax
    fname = "list_rows_max_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_rows_tmax,fname)
    
    ##Write out tmax
    fname = "list_rows_mean_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_rows_tmean,fname)

    #list_rows_tmin = load_data("list_rows_min_05032014.dat")
    #list_rows_tmax = load_data("list_rows_max_05032014.dat")
    #list_rows_tmean = load_data("list_rows_mean_05032014.dat")
    
    #list_rows_tmean = load_data("list_rows_mean_"+out_suffix)
    
    #li    dict_rast_fname = {"tmin": l_f_tmin, "tmax": l_f_tmax, "tmean": l_f_tmean }
    #list_rows_summary = [list_rows_tmin,list_rows_tmax,list_rows_tmean]
    #dict_rows_summary = {"tmin":list_rows_tmin,"tmax":list_rows_tmax,"tmean":list_rows_tmean}

    dict_rows_summary = {"tmin":list_rows_tmin,"tmax":list_rows_tmax,"tmean":list_rows_tmean}
    dict_table = {}

    ### Create name for variable computed: This part can change significantly...since it is specific to
    ### ncar filename.
    ##quite long...
    
    #There are 16*48 (48 variables and 16 region)...the final data.frame has 16*48 rows...ß
    nb_region = len(df_region1[region_name].unique()) #number of unique regions...16 for counties
    list_val_info =[]
    for i in range(0,len(dict_rows_summary)):
            ##MAKE THIS A SEPARATE PART IN WHICH valueType and Time are added!!
        l_f = dict_rast_fname.values()[i]           
        var_info = create_var_names_from_files(l_f)
        df_info = pd.DataFrame(var_info)
                
        val_info = []
        #name_col = []
        for j in range(0,3): #number of col in df_info
            list_name_col = []
            for k in range(0,df_info.count()[1]): #number of rows in df_info
                l = []  
                l.append(df_info.ix[k,j])
                name_col = l*nb_region
                list_name_col.extend(name_col)
            val_info.append(list_name_col)
        list_val_info.append(val_info)
       
    list_name_var = []
 
    for i in range(0,len(dict_rows_summary)):
        name_var = []
        val_info_tmp = list_val_info[i]
        for k in range(0,len(val_info[0])):
            test = val_info_tmp[0][k]+"_"+val_info_tmp[1][k]+"_"+val_info_tmp[2][k]
            name_var.append(test)
        list_name_var.append(name_var)
        
    ### end of name creation part for ncar temp predictions: might be changed for different input
    
    ### Now add info to create table
    
    #zonal_stat = "mean" #This is set earlier...
    zonal_stat_col = ["region","sum","mean","count","max","min"] #should be set earlier from teh output or rows?
    for i in range(0,len(dict_rows_summary)):
        list_rows = dict_rows_summary.values()[i]
        var_name = dict_rows_summary.keys()[i]
        nb_columns = len(list_rows) # + 1
        nb_rows = len(list_rows[0])
        element_nb = zonal_stat_col.index(zonal_stat)  #find the position of zonal stat in hte list
        #this is the column containing the mean for each entity (third column from postgis table)    
        #table_test = pdb.runcall(create_summary_table,list_rows_tmin,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        table = create_summary_table(list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        dict_table[var_name] = table
        #Now add id back...
        nb_columns = 1 #this regturns only the first
        element_nb = 0 #this is the column containing the region id set in the first part of the script
        id_reg = create_summary_table(list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        region_id = ["%.0d" % number for number in id_reg]
        
        #Now reshape to have a 48*16 rows and 4 columns... the last should be ID region
        nb_columns = len(list_rows) # + 1
        tt=table.ravel()
        df_val = pd.DataFrame(tt,columns=["value"])
        df_val["type"]= (list(region_type))* int(df_val.count()) #using the number of count in column value
        df_val["Name"]= region_id*nb_columns #using the number of count in column value
        df_val["valueType"] = list_name_var[i]
        df_val.ix[1:10,] #print first 10 rows with columns name for quick check
        #Write out results
        df_val.to_csv("table_"+var_name+"_df_"+out_suffix+".csv",sep=",") #write out table        
        
    #Write function to add back the tables in postgis and join them?   
     ############# NLCD SUMMARY TABLE this works for now but need to be improved!!!
    list_rows_NLCD
    l_f_valueType_NLCD
    
    df create_var_names_from_files_NLCD(f_list)   
    
        #There are 16*48 (48 variables and 16 region)...the final data.frame has 16*48 rows...ß
    nb_region = len(df_region1[region_name].unique()) #number of unique regions...16 for counties
    list_val_info_NLCD =[]
    #for i in range(0,len(dict_rows_summary)):
            ##MAKE THIS A SEPARATE PART IN WHICH valueType and Time are added!!
    #l_f = dict_rast_fname.values()[i]           

    var_info_NLCD = create_var_names_from_files_NLCD(l_f_valueType_NLCD)
    df_info_NLCD = pd.DataFrame(var_info_NLCD)
                
    ## Should this be a function similar to NCAR?            
    val_info_NLCD = []
        #name_col = []
    for j in range(0,3): #number of col in df_info:
        list_name_col = []
        for k in range(0,df_info_NLCD.count()[1]): #number of rows in df_info:
            l = []  
            l.append(df_info_NLCD.ix[k,j])
            name_col = l*nb_region
            list_name_col.extend(name_col)
            val_info.append(list_name_col)
        list_val_info_NLCD.append(val_info)
    
    list_name_var = []    
    name_var = []
    val_info_tmp = list_val_info_NLCD
    for k in range(0,len(val_info[0])):
        test = val_info_tmp[0][k]+"_"+val_info_tmp[1][k]+"_"+val_info_tmp[2][k]
        name_var.append(test)
        list_name_var.append(name_var)

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
        




    dict_rows_summary = {"nlcd":list_rows_NLCD} #,"tmax":list_rows_tmax,"tmean":list_rows_tmean}
    dict_table = {}
    #l_f = dict_rast_fname.values()[i]           
    dict_rast_fname = {}
    #l_f_tmin = filter(lambda x: re.search(r'tmin',x),l_f_valueType)        
    #l_f_tmax = filter(lambda x: re.search(r'tmax',x),l_f_valueType)
    #l_f_tmean = filter(lambda x: re.search(r'tmean',x),l_f_valueType)        
    dict_rast_fname = {"nlcd": l_f_valueType_NLCD}  #, "tmax": l_f_tmax, "tmean": l_f_tmean }

    ### Create name for variable computed: This part can change significantly...since it is specific to
    ### ncar filename and nlcd
    ##quite long...
    
    ######### THIS WHOLE SECTION CAN BE MADE AS A FUNCTION.... 
    #There are 16*48 (48 variables and 16 region)...the final data.frame has 16*48 rows...ß
    nb_region = len(df_region1[region_name].unique()) #number of unique regions...16 for counties
    list_val_info =[]
    for i in range(0,len(dict_rows_summary)):
            ##MAKE THIS A SEPARATE PART IN WHICH valueType and Time are added!!
        l_f = dict_rast_fname.values()[i]           
        #var_info = create_var_names_from_files(l_f) #this lines changes
        var_info = create_var_names_from_files_NLCD(l_f)
        
        df_info = pd.DataFrame(var_info)
                
        val_info = []
        #name_col = []
        for j in range(0,3): #number of col in df_info
            list_name_col = []
            for k in range(0,df_info.count()[1]): #number of rows in df_info
                l = []  
                l.append(df_info.ix[k,j])
                name_col = l*nb_region
                list_name_col.extend(name_col)
            val_info.append(list_name_col)
        list_val_info.append(val_info)
       
    list_name_var = []
 
    for i in range(0,len(dict_rows_summary)):
        name_var = []
        val_info_tmp = list_val_info[i]
        for k in range(0,len(val_info[0])):
            test = val_info_tmp[0][k]+"_"+val_info_tmp[1][k]+"_"+val_info_tmp[2][k]
            name_var.append(test)
        list_name_var.append(name_var)
        
    ### end of name creation part for ncar temp predictions: might be changed for different input
    
    ### Now add info to create table
    
    #zonal_stat = "mean" #This is set earlier...
    zonal_stat_col = ["region","sum","mean","count","max","min"] #should be set earlier from teh output or rows?
    for i in range(0,len(dict_rows_summary)):
        list_rows = dict_rows_summary.values()[i]
        var_name = dict_rows_summary.keys()[i]
        nb_columns = len(list_rows) # + 1
        nb_rows = len(list_rows[0])
        element_nb = zonal_stat_col.index(zonal_stat)  #find the position of zonal stat in hte list
        #this is the column containing the mean for each entity (third column from postgis table)    
        #table_test = pdb.runcall(create_summary_table,list_rows_tmin,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        table = create_summary_table(list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        dict_table[var_name] = table
        #Now add id back...
        nb_columns = 1 #this regturns only the first
        element_nb = 0 #this is the column containing the region id set in the first part of the script
        #table_test = pdb.runcall(create_summary_table,list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        id_reg = create_summary_table(list_rows,nb_rows,nb_columns,element_nb,out_dir,out_suffix,use_str=True)
        id_reg = id_reg.ravel()
        region_id = id_reg.tolist() #= region_id
        #region_id = ["%s" % number for number in id_reg]
        
        #Now reshape to have a 48*16 rows and 4 columns... the last should be ID region
        nb_columns = len(list_rows) # + 1
        tt=table.ravel()
        df_val = pd.DataFrame(tt,columns=["value"])
        df_val["type"]= (list(region_type))* int(df_val.count()) #using the number of count in column value
        df_val["Name"]= region_id*nb_columns #using the number of count in column value
        df_val["valueType"] = list_name_var[i]
        df_val.ix[1:10,] #print first 10 rows with columns name for quick check
        #Write out results
        df_val.to_csv("table_"+var_name+"_df_"+out_suffix+".csv",sep=",") #write out table        

    return None
    
#Need to add this to run
if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
#
#        df_var["decade"] = df_info["decade"]
#        df_var["month"] = df_info["month"]
#        df_var["var"] = df_info["var"]

################   END OF SCRIPT  ###############