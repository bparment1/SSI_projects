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
# Updated on: 06/24/2014
# Project: DSS-SSI
#
# TODO:
# - Need to add more functions by breaking out code!!...--> general call to summary
# - improve performance...by using in multiprocessing pool
# - add option to use raster by raster method for town summaries...
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
import pandas.io.sql as sql  #Direct access to database with ouput in DataFrame

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

def caculate_zonal_statistics(i,out_dir,out_suffix,rast_fname,shp_fname,SRS_EPSG,region_id_col,postgres_path_bin,db_name,user_name,NA_flag_val,var_name=None,tile_size=1):
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

    SQL_str = "SELECT * FROM %s" % (zonal_stat_regions)
    df_zonal = sql.read_frame(SQL_str,conn)     
    if var_name==None:
        var_name= "var"+"_"+str(i)
    df_zonal["var"] = var_name #adding the variable name
    df_zonal.to_csv("df_zonal_"+var_name+"_"+out_suffix+".csv")    
    #df_val.to_csv("table_"+var_name+"_df_"+out_suffix+".csv",sep=",") #write out table        

        
    #Also export in panda data.frame format if it is available on the system (add default option to function!!!)
    #cur.execute("select instrument, price, date from my_prices")
    #df = DataFrame(cur.fetchall(), columns=['instrument', 'price', 'date'])
    #df.set_index('date', drop=False)
    #df=sql.read_frame("select * from test",con)
    
    #list_rows.append(rows)
    #writeout list....    
    fname = "rows"+str(i)+"_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(rows_table,fname)
    
    #Close connection with postgres
    cur.close()
    conn.close()
    
    return rows_table,df_zonal
    
    
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
    #shp_fname = county24.shp #run for county
    shp_fname = "metwp24.shp" # run for town...
    shp_fname = os.path.join(in_dir,shp_fname) 
    #region_name = "COUNTY" #name of the column containing the id for each entities...
    region_name = "TOWN" #name of the column containing the id for each entities...
    
    variable_name = "temp" #type of variable can also be nlcd
    #region_name = "CNTYCODE"
    #if region_name is None then create an ID? Add this option later on.
    region_type = "T"    #type for the region entity: C for county, T for town  
    zonal_stat = "means" #This is the statistic extracted from the region

    #polygon_input_shp_file = "metwp24.shp"

    #EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
    #+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
    out_suffix = "06042014"
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

    SRS_EPSG = "26919"
    #Database information       
    db_name ="test_ssi_db3"
    db_name ="test_ssi_db2"
    
    user_name = "benoit"
    #user_name = "parmentier" 
    postgres_path_bin = "/usr/lib/postgresql/9.1/bin/"  #on SSI Maine
    tile_size = 1  #this set the tile size for the raster postgis table
    NA_flag_val = -9999 #this is not the case of ser and se data since it is float64 
    #postgres_path_bin = ""  #on ATLAS NCEAS
    
    #polygon_input_shp_table = "counties"    
    #polygon_input_shp_table = "towns"
    
    ## Get input raster lists containg information to be summarized
    #Temperature info:
    fileglob_pattern = "*projected*ncar*masked*.tif"
    pathglob = os.path.join(out_dir, fileglob_pattern)
    l_f_valueType = glob.glob(pathglob) #this contains the raster variable files that need to be summarized
    l_f_valueType.sort(key=natural_keys) #mixed sorting of files
    #NLCD land cover info:
    fileglob_pattern = "nlcd_*_*_proportion_900_900_*masked*.tif"
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
    lf_rast_fname = []
    l_f_tmin_7=  filter(lambda x: re.search(r'tmin_7_',x),l_f_valueType)       
    l_f_tmin_1 = filter(lambda x: re.search(r'tmin_1_',x),l_f_valueType)
    l_f_tmax_7 = filter(lambda x: re.search(r'tmax_7_',x),l_f_valueType)       
    l_f_tmax_1 = filter(lambda x: re.search(r'tmax_1_',x),l_f_valueType)
    l_f_tmean_7 = filter(lambda x: re.search(r'tmean_7_',x),l_f_valueType)       
    l_f_tmean_1 = filter(lambda x: re.search(r'tmean_1_',x),l_f_valueType)
    lf_rast_fname = l_f_tmin_7 + l_f_tmin_1 + l_f_tmax_1 + l_f_tmax_7 + l_f_tmean_1 + l_f_tmean_7 

    var_info = create_var_names_from_files(lf_rast_fname)
    df_info = pd.DataFrame(var_info)
    type_var_name = df_info['decade'] +"_" +df_info['month']+"_"+df_info["var"] #panda series with name of var
    
    #This will be parallelized and looped through dict_rast_fname 
    #can make one additional loop to reduce the repitition with the variable list of files

    list_rows_var = [] # defined lenth right now
    list_df_var = [] # defined lenth right now
    
    for i in range(0,len(lf_rast_fname)):
    #for i in range(0,2):
        var_name = type_var_name[i]
        out_suffix_s = var_name + "_"+out_suffix
        #test, test_df = pdb.runcall(caculate_zonal_statistics,i,out_dir,out_suffix_s,rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name)
        rows,df = caculate_zonal_statistics(i,out_dir,out_suffix_s,lf_rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name,NA_flag_val,var_name,tile_size)
        list_df_var.append(df) #add object rows to list
        list_rows_var.append(rows) #add object rows to list
        
    ##Write out rows
    fname = "list_rows_"+region_name+"_"+variable_name+"_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_rows_var,fname)
    ##Write out dataframe
    fname = "list_df_"+region_name+"_"+variable_name+"_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_df_var,fname)

    df_var = pd.concat(list_df_var) #problem with the row indices
    nb_rows = df_var['means'].count()
    new_index = np.arange(0,nb_rows).tolist()  
    df_var['ni'] = new_index #change the rows index
    df_var = df_var.set_index('ni')
    #df_var[df_var['means']> 200][['county','means']] #this is an example of subset of data.frame
    #note that there is a problem with name county transformed to lower case
    df_val = df_var[[region_name.lower(),zonal_stat,'var']] #select specific columns in a dataframe
    df_val.columns = ["Name","value","valueType"] #rename columns 
    df_val["type"] = [region_type]*nb_rows #create new column...
    variable_name = "temp" #set up at the beginning
    df_val.to_csv("df_zonal_combined_"+region_name+"_"+variable_name+"_"+out_suffix+".csv")
    df_val.ix[1:10,] #print first 10 rows with columns name for quick check
    df_val["value"].describe() #summary values and range for checking output...
    
    ############ Summarize NLCD data ####################

    #inputs
    var_info_NLCD = create_var_names_from_files_NLCD(l_f_valueType_NLCD)
    df_info = pd.DataFrame(var_info_NLCD) 
    type_var_name = df_info['category'] +"_" +df_info['var']+"_"+df_info["year"] #panda series with name of var
    variable_name =  "NLCD"
    lf_rast_fname = l_f_valueType_NLCD #copy by refernce!!
    #region_name, zonal_stat + all the inputs  for calculate_zonal_statistics
    
    ## Could wrap this into a nice function....
    
    #This will be parallelized and looped through dict_rast_fname 
    #can make one additional loop to reduce the repitition with the variable list of files

    list_rows_var = [] # defined lenth right now, will contains rows from SQL query
    list_df_var = [] # defined lenth right now, will contain dataframe from sql query
    
    for i in range(0,len(lf_rast_fname)):
    #for i in range(0,2):
        var_name = type_var_name[i]
        out_suffix_s = var_name + "_"+out_suffix
        #NA_flag_val_rast = 1.175494351e-38
        NA_flag_val_rast = getNoDataValue(lf_rast_fname[i])
        #test, test_df = pdb.runcall(caculate_zonal_statistics,i,out_dir,out_suffix_s,rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name)
        rows,df = caculate_zonal_statistics(i,out_dir,out_suffix_s,lf_rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name,NA_flag_val_rast,var_name,tile_size)
        list_df_var.append(df) #add object rows to list
        list_rows_var.append(rows) #add object rows to list
        
    ##Write out rows
    fname = "list_rows_"+region_name+"_"+variable_name+"_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_rows_var,fname)
    ##Write out dataframe
    fname = "list_df_"+region_name+"_"+variable_name+"_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_df_var,fname)

    df_var = pd.concat(list_df_var) #problem with the row indices
    nb_rows = df_var['means'].count()
    new_index = np.arange(0,nb_rows).tolist()  
    df_var['ni'] = new_index #change the rows index
    df_var = df_var.set_index('ni')
    #df_var[df_var['means']> 200][['county','means']] #this is an example of subset of data.frame
    #note that there is a problem with name county transformed to lower case
    df_val = df_var[[region_name.lower(),zonal_stat,'var']] #select specific columns in a dataframe
    df_val.columns = ["Name","value","valueType"] #rename columns 
    df_val["type"] = [region_type]*nb_rows #create new column...
    #variable_name = "temp" #set up at the beginning
    df_val.to_csv("df_zonal_combined_"+region_name+"_"+variable_name+"_"+out_suffix+".csv")
    df_val.ix[1:10,] #print first 10 rows with columns name for quick check
    df_val["value"].describe() #summary values and range for checking output...

    ############ Summarize Housing density and population data ####################

    #ser00_housing_clipped_projected_06042014.tif
    lf_var_pop = glob.glob(os.path.join(out_dir,"*pop_clipped_projected*masked*.tif"))
    lf_var_hou = glob.glob(os.path.join(out_dir,"*housing_clipped_projected*masked*.tif"))
    lf_var_hist_temp = glob.glob(os.path.join(out_dir,"*me_hist_*_clipped_projected*masked*.tif"))
    lf_var = lf_var_pop + lf_var_hou + lf_var_hist_temp
      
#['/ssi-dss-data/DSS_SSI_data/output_data_06042014/ser00_pop_clipped_projected_06042014.tif',
# '/ssi-dss-data/DSS_SSI_data/output_data_06042014/ser10_pop_clipped_projected_06042014.tif',
# '/ssi-dss-data/DSS_SSI_data/output_data_06042014/ser90_pop_clipped_projected_06042014.tif',
 #'/ssi-dss-data/DSS_SSI_data/output_data_06042014/ses00_pop_clipped_projected_06042014.tif',
 #'/ssi-dss-data/DSS_SSI_data/output_data_06042014/ser10_housing_clipped_projected_06042014.tif',
 #'/ssi-dss-data/DSS_SSI_data/output_data_06042014/ser90_housing_clipped_projected_06042014.tif',
 #'/ssi-dss-data/DSS_SSI_data/output_data_06042014/ser00_housing_clipped_projected_06042014.tif']
    
    type_var_name = pd.Series(["ser_pop_2000","ser_pop_2010","ser_pop_1990","ses_pop_2000",
                              "ser_housing_2000","ser_housing_2010","ser_housing_1990",
                              "hist_tmin","hist_tmean","hist_tmax"])
            
    #inputs
    #var_info_NLCD = create_var_names_from_files_NLCD(l_f_valueType_NLCD)
    #df_info = pd.DataFrame(var_info_NLCD) 
    #Ser_var_name = df_info['category'] +"_" +df_info['var']+"_"+df_info["year"] #panda series with name of var
    variable_name =  "ser"
    #make a series from the list:
    lf_rast_fname = lf_var #copy by refernce!!
    #region_name, zonal_stat + all the inputs  for calculate_zonal_statistics
    
    ## Will be wrap in a function at later stage of the code !!
    
    #This will be parallelized and looped through dict_rast_fname 
    #can make one additional loop to reduce the repitition with the variable list of files

    list_rows_var = [] # defined lenth right now, will contains rows from SQL query
    list_df_var = [] # defined lenth right now, will contain dataframe from sql query

    for i in range(0,len(lf_rast_fname)):
    #for i in range(0,2):
        var_name = type_var_name[i]
        out_suffix_s = var_name + "_"+out_suffix
        #test, test_df = pdb.runcall(caculate_zonal_statistics,i,out_dir,out_suffix_s,rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name)
        rows,df = caculate_zonal_statistics(i,out_dir,out_suffix_s,lf_rast_fname,out_fname1,SRS_EPSG,region_name,postgres_path_bin,db_name,user_name,NA_flag_val,var_name,tile_size)
        list_df_var.append(df) #add object rows to list
        list_rows_var.append(rows) #add object rows to list
        
    ##Write out rows
    fname = "list_rows_"+region_name+"_"+variable_name+"_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_rows_var,fname)
    ##Write out dataframe
    fname = "list_df_"+region_name+"_"+variable_name+"_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_df_var,fname)

    df_var = pd.concat(list_df_var) #problem with the row indices
    nb_rows = df_var['means'].count()
    new_index = np.arange(0,nb_rows).tolist()  
    df_var['ni'] = new_index #change the rows index
    df_var = df_var.set_index('ni')
    #df_var[df_var['means']> 200][['county','means']] #this is an example of subset of data.frame
    #note that there is a problem with name county transformed to lower case
    df_val = df_var[[region_name.lower(),zonal_stat,'var']] #select specific columns in a dataframe
    df_val.columns = ["Name","value","valueType"] #rename columns 
    df_val["type"] = [region_type]*nb_rows #create new column...
    #variable_name = "temp" #set up at the beginning
    df_val.to_csv("df_zonal_combined_"+region_name+"_"+variable_name+"_"+out_suffix+".csv")
    df_val.ix[1:10,] #print first 10 rows with columns name for quick check
    df_val["value"].describe() #summary values and range for checking output...
    

    return None
    
#Need to add this to run
if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-

################   END OF SCRIPT  ###############