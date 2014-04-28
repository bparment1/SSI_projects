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
# Updated on: 04/28/2014

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
import numpy                 #Array, matrices and scientific computing
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

def create_summary_table(list_rows,nb_rows,nb_columns, element_nb,out_dir,out_suffix):
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
   
def caculate_zonal_statistics(i,out_dir,out_suffix,rast_fname,shp_fname,
    SRS_EPSG,postgres_path_bin,db_name,user_name):
    
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

    SQL_str ="WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, cntycode FROM %s, %s WHERE ST_Intersects(pixelgeom,geom)) SELECT cntycode, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by cntycode;" % (poly_table,output_table_shp)

    #cmd_str ="WITH temp_table AS (SELECT ST_Area((ST_Intersection(pixelgeom,geom)).geometry) AS intersectarea, pixelval, pixelarea AS origarea, cntycode FROM mht_poly_use, counties WHERE ST_Intersects(pixelgeom,geom)) SELECT cntycode, SUM(pixelval*intersectarea/origarea) AS sums, SUM(pixelval*intersectarea)/SUM(intersectarea) AS means, COUNT(*) AS counts, MAX(pixelval) AS maxes, MIN(pixelval) AS mins, SUM(intersectarea)/1000000 AS area from temp_table group by cntycode;"
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
    
    return(rows_table)
    
    
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
    #polygon_input_shp_file = "metwp24.shp"

    #EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
    #+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
    out_suffix = "04122014"
    os.chdir(in_dir) #set current dir
    #os.getcwd() #check current dir
    out_dir = in_dir

    #out_path<-"/data/project/layers/commons/data_workflow/output_data"
    out_dir = "output_data_"+out_suffix
    out_dir = os.path.join(in_dir,out_dir)
    create_dir_and_check_existence(out_dir)
    
    #EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
    #+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
    #file_format = ".tif"

    SRS_EPSG = "2037"             
    #Database information       
    db_name ="test_ssi_db2"
    #user_name = "benoit"
    user_name = "parmentier" 
    #postgres_path_bin = "/usr/lib/postgresql/9.1/bin/"
    postgres_path_bin = ""
    
    #polygon_input_shp_table = "counties"    
    #polygon_input_shp_table = "towns"
    
    ########## START SCRIPT #############
         
    ## Temp processing: 2020s,2030s,2040s,2050s for tmin, tmax and tmean (e.g. 4*3 directories with 12 files...)

    #loop through files...
    fileglob = "ncar*projected.tif"
    pathglob = os.path.join(in_dir, fileglob)
    l_f_temp = glob.glob(pathglob)
    l_f_temp.sort(key=natural_keys) #mixed sorting of files
    #Get tmin files
    l_f_tmin = filter(lambda x: re.search(r'tmin',x),l_f_temp)
    #Get tmax files
    
    list_rows_tmin = [] # defined lenth right now
    #i=1
    #This will be parallielized
    #debug   
    rast_fname = l_f_tmin #copy by refernce!!
    #list_rows_tmin.append(i) = pdb.runcall(caculate_zonal_statistics,i,out_dir,out_suffix,rast_fname,shp_fname,

    #for i in range(0,len(rast_fname)):
    for i in range(0,2):      
        rows = caculate_zonal_statistics(i,out_dir,out_suffix,rast_fname,shp_fname,SRS_EPSG,postgres_path_bin,db_name,user_name)
        list_rows_tmin.append(rows) #add object rows to list

    #Write out results by combining
    #test1 =load_data(fname)
    fname = "list_rows_min_"+out_suffix+".dat"
    fname = os.path.join(out_dir,fname)
    save_data(list_rows,fname)
    
    out_dir = "/ssi-dss-data/DSS_SSI_data/output_data_04122014"
    list_rows_mean = load_data(os.path.join(out_dir,"list_rows_mean_04122014.dat"))
    
    nb_columns=len(list_rows_mean)
    nb_columns=48

    nb_rows=16
    #test = pdb.runcall(create_summary_table,list_rows_mean,nb_rows,nb_columns,out_dir,out_suffix)
    element_nb=2 #this is the column containing the mean for each entity (third column from postgis table)    
    test = create_summary_table(list_rows_mean,nb_rows,nb_columns,element_nb,out_dir,out_suffix)
        
    element_nb=0 #this is the column containing the mean for each entity (third column from postgis table)    
    id_reg = create_summary_table(list_rows_mean,nb_rows,nb_columns=1,element_nb,out_dir,out_suffix)
        
    table_mean = np.hstack(id_reg,)    
    return None
    
#Need to add this to run
if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
#
################   END OF SCRIPT  ###############