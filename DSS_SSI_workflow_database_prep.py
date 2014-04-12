#!/usr/bin/python
#
#Script to process data for the DSS-SSI website.
#Raster Data inserted in a Postgis database and tables are created. 
#clipped to match the Maine study area.
#
# TODO:
#  - functionalize to encapsulate high level procedural steps
#
# Authors: Benoit Parmentier 
# Created on: 04/02/2014
# Updated on: 04/10/2014

import os, glob, sys
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
import psycopg2
import numpy

#------------------
# Functions used in the script 
#------------------


def create_dir_and_check_existence(path):
    try:
        os.makedirs(path)
        
    except:
        print "directory already exists"
            
def save_data(data,fname):
    import pickle
    with open(fname, "wb") as f:
        pickle.dump(data, f)
        
def load_data(fname):
    import pickle
    try:
        with open(fname) as f:
            data_obj = pickle.load(f)
    except:
        print "No data"
        data_obj=[]
    return data_obj
      
      
#import re
#http://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def create_summary_table(list_rows,nb_rows,nb_columns, out_dir,out_suffix):
    #
    table = numpy.ones((16,48))
    #add ID first...?
    for i in range(0,48):
        rows_output = list_rows[i]
        for j in range(0,16):
            tu = rows_output[j]
            table[j,i] = tu[2]
    fname = "table_"+out_suffix+".txt"       
    outfile1 = os.join.path(in_dir,fname)
    numpy.savetxt(outfile1, table, fmt='%-7.6f')

        
    
## Need to add more functions by breaking out code!!...

## add function for use in multiprocessing pool

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
        
    in_dir = "/home/parmentier/Data/IPLANT_project/Maine_interpolation/DSS_SSI_data/"
    #in_dir = "/ssi-dss-data/DSS_SSI_data/"
    polygon_input_shp_file = "county24.shp"
    #polygon_input_shp_file = "metwp24.shp"

    #EPSG: http://spatialreference.org/ref/epsg/26919/proj4/ -->  
    #+proj=utm +zone=19 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
    #file_format = ".tif"
    out_suffix = "03242014"
    os.chdir(in_dir) #set current dir
    #os.getcwd() #check current dir
    out_dir = in_dir

    #out_path<-"/data/project/layers/commons/data_workflow/output_data"
    out_dir = "output_data_"+out_suffix
    out_dir = os.path.join(in_dir,out_dir)
    create_dir_and_check_existence(out_dir)
            
    #Database information       
    db_name ="test_ssi_db2"
    #user_name = "benoit"
    user_name = "parmentier" 
    #postgres_path_bin = "/usr/lib/postgresql/9.1/bin/"
    postgres_path_bin = ""
    
    polygon_input_shp_table = "counties"    
    #polygon_input_shp_table = "towns"
    
    ########## START SCRIPT #############
         
    ## Temp processing: 2020s,2030s,2040s,2050s for tmin, tmax and tmean (e.g. 4*3 directories with 12 files...)

    #loop through files...
    fileglob = "ncar*projected.tif"
    pathglob = os.path.join(out_dir, fileglob)
    l_f_temp = glob.glob(pathglob)
    l_f_temp.sort(key=natural_keys) #mixed sorting of files
    l_f_tmin = filter(lambda x: re.search(r'tmin',x),l_f_temp)

    l_f = l_f_tmin
    #def grep_fun(l,s):
    #    return [i for i in l if re.search(r'aet',i)]
        
    #def grep_fun(pattern,word_list):
    #    expr = re.compile(pattern)
    #    return [elem for elem in word_list if expr.match(elem)]

    #l_f_tmean = grep_fun(l_f,"tmean")

    #Need to sort this by inputs...tmin, tmax and tmean

    #Write a function to excute the next lines in a loop (for paralelization)
    #Multiprocessing module requires individual db connections for each process 
    
    #INPUT PARAMETERS

    #postgres_path_bin = "/usr/lib/postgresql/9.1/bin/"
    #Database information       
    #db_name ="test_ssi_db2"
    #user_name = "benoit"
    #in_dir
    
    list_rows = []
    #i=1
    SRS_EPSG = "2037" 
        
    for i in range(0,len(l_f)):

        rast_input = l_f[i]
        #output_table_rast = "mtemp"
        output_table_rast ="".join(["mtemp_",str(i)])
        poly_table = "".join(["mht_poly_use_",str(i)]) #This is the name of the table that will store the values of 
        #output_table_rast ="mtemp"
        #poly_table = "mht_poly_use" #This is the name of the table that will store the values of 
        
        #polygon_input_shp_file
        shp_input = os.path.join(in_dir,polygon_input_shp_file)
        output_table_shp = "".join([polygon_input_shp_table,str(i)])
        #output_table_shp = polygon_input_shp_table

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
        rows = cur.fetchall()
        
        list_rows.append(rows)
        #writeout list....
        
        fname = "rows"+str(i)+"_"+out_suffix+".dat"
        fname = os.path.join(out_dir,fname)
        save_data(rows,fname)
        
        #test1 =load_data(fname)
        fname = "list_rows_"+out_suffix+".dat"
        fname = os.path.join(out_dir,fname)
        save_data(list_rows,fname)
        
        #return rows
                
        #Parse in Python or have a table in Postgres...        
        #return
        #end of function
    
    #Write out results by combining
    
    create_summary_table(list_rows,in_dir,out_suffix)
        

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

# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 07:08:47 2014

@author: parmentier
"""

