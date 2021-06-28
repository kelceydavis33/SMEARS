#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
import glob
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from math import e
import os
import sys
from scipy import io
from astropy.time import Time
from astropy import wcs
import pandas as pd
import skimage
import scipy as sp
from skimage.feature import register_translation 
from matplotlib import cm, colors
from matplotlib.collections import PatchCollection
import scipy.misc
import pyradiosky
from pyradiosky import SkyModel
from astropy.coordinates import Longitude, Latitude
from IPython.display import clear_output
def match_to_gleam(directory, GLEAM_path, index_lim = None, extra_flux_table = False, missed_source_table = False):
    """
    Matches all the .sav files in a given directory with GLEAM files and creates a Pandas DataFrame from 
    the data
    
    Parameters
    ---
    
    directory: string
                an absolute path to the data files that have been proccessed through FHD
    
    GLEAM_path: string
                an absolute path to a saved GLEAM .sav file
    
    index_lim:  integer, optional
                index cut-off of the last observation you want processed. This is helpful in limiting
                time code takes to run and exploring small chunks of your data
    
    Returns
    ---
    df
        A pandas DataFrame containing matching information for sources. Each row is an object in the input
        GLEAM file. The columns represent the input fies matched to GLEAM. The names of the rows indicate 
        which data file is being matched, indicated by an integer index based on it's numerical order in the
        input directory. The columns indicate the following information:
        
            RA:
                original right asention coordinate value in decimal degrees specified in GLEAM file
            
            DEC:
                original declination coordinate value in decimal degrees specified in GLEAM file
            
            Mag GLEAM:
                original magnitude in Jansky of the object as specified in the GLEAM file
            
            RA EO Gleam:
                right asention in degrees for extended components of the object as specified in the GLEAM
                file. If no extended components are specified, this values is zero.
            
            DEC EO GLEAM
                declination in degrees for extended components of the object as specified in the GLEAM
                file. If no extended components are specified, this values is zero.
            
            Mag EO GLEAM
                magnitude in Jansky for extended components of the object as specified in the GLEAM
                file. If no extended components are specified, this values is zero.
                
            The remaining columns include the integer identifier of the data file being matched to gleam.
            These columns have the following form:
            
            Mag {integer identifier}:
                The point-like high-level value for the magnitude of the source in the data file that best
                matched the GLEAM source specified by the row. If no match was made, a 0 is added in this
                position.
            
            Distance {integer identifier}:
                This is the matching distance in decimal degrees for the match indicated in the data. For
                cases where no match was made, the closest matchin distance is added
            
            RA {integer identifier}:
                The point-like high-level value for the right asention in degrees of the source in the data 
                file that best matched the GLEAM source specified by the row. If no match was made, a 0 
                is added in this position.
            
            DEC {integer identifier}:
                The point-like high-level value for the declination in degrees of the source in the data 
                file that best matched the GLEAM source specified by the row. If no match was made, a 0 
                is added in this position.
            
            EO RA {integer identifier}:
                The right ascention in decimal degrees of the low-level extended components of the object 
                in the data files being processed that best matched the GLEAM object specefied by the row. 
                If no match was made or no extended components are present, a 0 is added in this position.
            
            EO DEC {integer identifier}:
                The declination in decimal degrees of the low-level extended components of the object in 
                the data files being processed that best matched the GLEAM object specefied by the row. If 
                no match was made or no extended components are present, a 0 is added in this position.
            
            EO MAG {integer identifier}:
                The magnitude in Janksy of the low-level extended components of the object in 
                the data files being processed that best matched the GLEAM object specefied by the row. If 
                no match was made or no extended components are present, a 0 is added in this position.
            
            STON {index identifier}:
                The signal to noise ratio of the high-level point-like source indicated by the data file 
                for the indicated best match source.
            
            STON EO {index identifier}:
                The signal to noise ratio of the low-level extended points of the source indicated by the 
                data file for the indicated best match source.
            
    """
    #Create paths to all the .sav files in the specified directory
    paths = glob.glob(directory + '*source_array.sav')
    #Specify the matching distance
    match_dist = 0.1
    #Specify the fraction of flux of the GLEAM magnitude that a source must be to be considered a match
    match_frac = 3/4
    #Define empty lists to store the GLEAM data
    ra_gleam = []
    dec_gleam = []
    imag_gleam = []
    #Print an indicarion that the GLEAM file is being opened 
    print('Opening GLEAM File', end = "\r")
    #Create the path to the GLEAM file
    gpath = glob.glob(GLEAM_path)
    #Create an empty dictionary to store the GLEAM data
    GLEAM_data = {}
    #Open the data from the GLEAM file
    GLEAM_data['data'] =  [scipy.io.readsav(gpath[i], python_dict=True)
            for i in range(len(gpath))]
    #Find the length of the GLEAM file
    size = GLEAM_data['data'][0]['catalog'].shape[0]    #Create empty lists to store the extended objects in the GLEAM file
    ra_gleam_ext = []
    dec_gleam_ext = []
    imag_gleam_ext = []
    #Create a list to store the names associated with the GLEAM data
    gnams = []
    #Check if the data type is correct for interpretation
    if type(GLEAM_data['data'][0]['catalog'][0][0]) == bytes:
        #Loop over the sources in GLEAM
        for z in range(size):
            #Pull out the identifier for the source
            s = GLEAM_data['data'][0]['catalog'][z]['ID']
            #Get the important part of the identifier and make it a string
            nam = str(s)[2:-1]
            #Append the name to the GLEAM identifiers
            gnams.append(nam)
    #Load the GLEAM data
    for i in range(size):
        #Check if there are no extended components
        if GLEAM_data['data'][0]['catalog'][i]['EXTEND'] is None:
            #append RA to appropriate gleam list
            ra_gleam.append(GLEAM_data['data'][0]['catalog'][i]['RA'])
            #append DEC to appropriate gleam list
            dec_gleam.append(GLEAM_data['data'][0]['catalog'][i]['DEC'])
            #append Magnitude to appropriate gleam list
            imag_gleam.append(GLEAM_data['data'][0]['catalog'][i]['FLUX']['I'])
            #append a 0 to indicate there are no extened components to apropriate lists
            ra_gleam_ext.append(0)
            dec_gleam_ext.append(0)
            imag_gleam_ext.append(0)
        else:
            #If there are extended sources
            #append RA ,DEC, and magnitude to appropriate GLEAM list
            ra_gleam.append(GLEAM_data['data'][0]['catalog'][i]['RA'])
            dec_gleam.append(GLEAM_data['data'][0]['catalog'][i]['DEC'])
            imag_gleam.append(GLEAM_data['data'][0]['catalog'][i]['FLUX']['I'])
            #create an empty lsit to store the extended magnitude components
            glm_mag_eo = []
            #Loop through and pick out the extended points in I flux then append these to the list
            for j in range(0, GLEAM_data['data'][0]['catalog'][i]['FLUX'].shape[0]):
                glm_mag_eo.append(GLEAM_data['data'][0]['catalog'][i]['EXTEND']['FLUX'][j]['I'])
            imag_gleam_ext.append(glm_mag_eo)
            #append the values for extended RA and DEC components
            ra_gleam_ext.append(GLEAM_data['data'][0]['catalog'][i]['EXTEND']['RA'])
            dec_gleam_ext.append(GLEAM_data['data'][0]['catalog'][i]['EXTEND']['DEC'])
    #Create a Pandas Data Frame with the RA, DEC, and GLEAM Magnitudes
    n=0
    df = pd.DataFrame({'RA': ra_gleam,'Mag GLEAM': imag_gleam,  'DEC' : dec_gleam, 'RA EO GLEAM':ra_gleam_ext, 
                      'DEC EO GLEAM': dec_gleam_ext, 'Mag EO GLEAM': imag_gleam_ext})
    #Look at each path in the directory
    fl = 0
    #Check if the optional keyword to grab only a few files has been entered, if so apropriately 
    #index the data
    iter_paths = paths[0:index_lim] if index_lim else paths
    for path in iter_paths:
        fl +=1
        os.system('clear')
        print(f"Generating table for file {fl} of {index_lim}", end = "\r")
        #Collect the data for each path
        n = n + 1
        #Identify path to this data file
        datpath = glob.glob(path)
        #Create a dictionary to store the data
        data = {}
        #Open the data
        data['data'] = [scipy.io.readsav(datpath[i], python_dict=True)
            for i in range(len(datpath))]
        #Create some empty dictionaries to store the data
        eo = []
        eo_ra = []
        eo_dec = []
        ps_RA = []
        ps_DEC = []
        i_mag = []
        EO_imag = []
        ps_ston = []
        eo_ston = []
        #Try using the 'catalog' keyword to open the data
        try:
            d_s = data['data'][0]['catalog']
        #If the 'catalog' keyword does not exist, assume that the keyword is 'source_array'
        except:
            d_s = data['data'][0]['source_array']
        #Loop over the sources in the data file
        for d in d_s:
            #Check if the source has extended, or multiple, components
            #If the source is not extended, store only the high-level
            #source information and disregard extended keywords, adding 
            #these high-level values to the lists containing extended information
            if d['EXTEND'] is None:
                ps_RA.append(d['RA'])
                ps_DEC.append(d['DEC'])
                #Flux in I is stored and other polarized fluxes are ignored
                EO_imag.append(d['FLUX']['I'])
                eo_ra.append(d['RA'])
                eo_dec.append(d['DEC'])
                #Flux in I is stored and other polarized fluxes are ignored
                i_mag.append(d['FLUX']['I'])
                ps_ston.append(d['STON'])
                eo_ston.append(d['STON'])
            #If the source does have extended components, add the high level point-
            #like values to high level lists and add extended components to
            #the extended lists
            else:
                ps_RA.append(d['RA'])
                ps_DEC.append(d['DEC'])
                #Create empty lists to store information about extended flux and signal to noise
                EOmags = []
                eoston = []
                #Loop over the extended flux to pull out only the relavent flux in I
                for i in range(0, d['EXTEND']['FLUX'].shape[0]):
                    #Flux in I is stored and other polarized fluxes are ignored
                    EOmags.append(d['EXTEND']['FLUX'][i]['I'])
                    eoston.append(d['EXTEND']['STON'][i])
                #Append these lists of 
                eo_ston.append(np.array(eoston))
                EO_imag.append(np.array(EOmags))
                eo_ra.append(d['EXTEND']['RA'])
                eo_dec.append(d['EXTEND']['DEC'])
                #Flux in I is stored and other polarized fluxes are ignored
                i_mag.append(d['FLUX']['I'])
                ps_ston.append(d['STON'])   
        #Match this path with the GLEAM catalog
        #idx: an array of indices corresponding to matches
        #d2d: the two dimensional distances between these matches
        #d3d: three dimensional distances between matches. This array is blank becasue we do 
        #not have 3 dimensional data, but the match_to_catalog_sky function requires it anyway
        catalog = SkyCoord(ra=ps_RA*u.deg, dec=ps_DEC*u.deg)  
        c = SkyCoord(ra=ra_gleam*u.deg, dec=dec_gleam*u.deg)  
        idx, d2d, d3d = c.match_to_catalog_sky(catalog)    
        #Only return matches within one degree 
        #Create an empty list to store the matched data
        mags = []
        #Sort the flux array with the the idx index array so that the magnitudes are ordered by match
        imags = np.array(i_mag)[idx]
        # create number indicies to loop over
        nums = np.arange(0, len(idx))
        #Loop through each index
        for num in nums:
            #Check if the matching distance is whithin the specefied distance
            #Check if the flux is withing a specified fraction of the GLEAM magnitude 
            #to specify a match
            if (d2d[num] < (match_dist*u.deg)) and ((imags[num] > imag_gleam[num] + (match_frac)*imag_gleam[num]) or (imags[num] > (match_frac)*imag_gleam[num])):
                mags.append(imags[num])
                #If there is no match, add a 0 in place of the flux vlaue
            else: 
                mags.append(0)
        #Add a new column to the data frame with the information from these observations 
        #Create a pandas data series of the flux array and the 2-dimensional distance arrays
        #which have already been ordered correctly
        s_mag = pd.Series(mags)
        s_dist = pd.Series(d2d)
        #Add a new column to the data frame for each list and fill it with the relavent data
        df['Mag {}'.format(n)] = s_mag
        df['Distance {}'.format(n)] = s_dist
        #Order the remaining information corectly by indexing them with the idx arrays
        #create a pandas dara series from these arrays to ensure that they are corectly formatted
        df['RA {}'.format(n)] = np.array(ps_RA)[idx]
        df['DEC {}'.format(n)] = np.array(ps_DEC)[idx]
        df['EO RA {}'.format(n)] = np.array(pd.Series(eo_ra))[idx]
        df['EO DEC {}'.format(n)] = np.array(pd.Series(eo_dec))[idx]
        df['EO Mag {}'.format(n)] = np.array(pd.Series(EO_imag))[idx]
        df['STON {}'.format(n)] = np.array(ps_ston)[idx]
        df['EO STON {}'.format(n)] = np.array(pd.Series(eo_ston))[idx]
    #change the index of the data frames to reflect the names of the objects specified by GLEAM
    df.index = gnams
    return df
def save_df(table, sav_path):
    """
    Parameters
    ---
    table: DataFrame
            a generated pandas DataFrame from the match_to_GLEAM function
    sav_path: string 
            a string absolute path to a folder where the file should be saved
    
    Returns
    ---
    A .pickle file saved to the specefied directory
    """
    table.to_pickle(sav_path + 'SMEARS_table.pickle')
def load_df(df_path):
    """
    Parameters
    ---
    df_path: string
            absolute path to where a pandas DataFrame generated by the match_to_GLEAM function
             has been generated
    Returns
    ---
    table: DataFrame
            the loaded pandas DataFrame
    """
    table = pd.read_pickle(df_path)
    return table
def gaussian(sigma, x_0, y_0, x, y, power):
    """
    Returns the value at a given x and y position for a gaussian surface centered at x_0, y_0
    
    Parameters
    ---
    
    sigma: float
           standard deviation of the gaussain beam
    
    x_0: float
         x value of peak 
    
    y_0: float 
        y value of peak
    
    x: float 
      the x value of the coordinate to be calculated
    
    y: flaot
       the y vlaue of the coordinate to be calculated
    
    power: float
           amplitude of gaussian at peak
    
    Returns
    ---
    gauss_val: float
                the value of the gaussian beam at the specified x, y position
    """
    #xs is the squared difference between x and x_0
    xs = (x-x_0)**2
    #ys is the squared difference between y and y_0
    ys = (y-y_0)**2
    #gauss_val is the vlaue of a gaussin beam at the specified point
    gauss_val = power * (e** ((-(xs + ys))/(2*(sigma**2))))
    return gauss_val
def shift_arrays(table, index, save_plots = True):
    """
    Shifts all arrays for a given observation towards an unweighted mean center. This is a good way to
    quickly check what a source looks like without running the smoothing which is more time consuming
    
    Parameters
    ---
    
    table: DataFrame
           a table generated by the match_to_gleam function
    
    index: string or integer
            the index of the target source and is an integer index or a string 
            
    save_plots: bool, optional
                Saves plots to current working directory
    
    Returns
    ---
    Displays and saves plot of source to the current working directory
    
    """
    #Settilng matplotlib to display a clean plot
    plt.rcParams['figure.figsize'] = (10, 10)
    plt.rc('axes', labelsize=14)
    plt.rc('axes', labelweight='bold')
    plt.rc('axes', titlesize=16)
    plt.rc('axes', titleweight='bold')
    plt.rc('font', family='sans-serif')
    #Pull relavent data from the pandas DataFrame
    mags = table.loc[index][6::9]
    mags_eo = table.loc[index][12::9]
    ra_center = table.loc[index][8::9]
    dec_center = table.loc[index][9::9]
    ext_ras = table.loc[index][10::9]
    ext_decs = table.loc[index][11::9]
    #Create empty lists to store shifted data
    ras = []
    decs = []
    mag_change = []
    adj_obs = []
    imnum = []
    true_eo_mags = []
    #Create a number to track iterations
    num = 0
    #Loop through the individual observations of the source
    for i in range(0, len(mags)):
        #Add to the iteration tracker
        num = num+1
        #Check if the mags array is 0. If it is, this 0, then it was identified
        #as not a match and we pass over the data
        if mags[i] !=0:
            #Append the values to the proper lists if recognized as a match
            mag_change.append(mags[i])
            ras.append(ra_center[i])
            decs.append(dec_center[i])
            imnum.append(num)
    #Find the center for the x and y values by finding the mean, unweighted
    x_center = np.mean(ras)
    y_center = np.mean(decs)
    #Loop through the observations again
    for i in range(0, len(mags)):
        #Check agian if the observation is a match
        if mags[i] !=0:
            #Grab the relavent RA and DEC values from the lists
            rs = ext_ras[i]
            ds = ext_decs[i]
            #Calculate the distance the RA and DEC values must be shifted
            ra_roll = x_center - ra_center[i]
            dec_roll = y_center - dec_center[i]
            #Create a new array containing the shifted RA and DEC vlaues
            new_array = rs+ra_roll, ds+dec_roll
            #Add the adjusted array to the list of new observations
            adj_obs.append(new_array)
            #Add the flux aray to the list of flux arrays
            true_eo_mags.append(mags_eo[i])
    #Loop over the identified and shifted observations
    for r in range(len(adj_obs)):
        #plot a scatter plot of the newly shifted observations
        plt.scatter(adj_obs[r][0], adj_obs[r][1],# s = true_eo_mags[r]*5, 
            label = 'Observation {}, {} Jy'.format(imnum[r],mag_change[r]))
    #Create X and Y labels, create title
    plt.xlabel('RA (Degrees)')
    plt.ylabel('DEC (Degrees)')
    plt.title('Source {}, Center Calculation'.format(index))
    #Identify the RA and DEC values originally identified in GLEAM
    g_ra = table.loc[index]['RA']
    g_dec = table.loc[index]['DEC']
    #Use these to adjust the X and Y limits of the plot
    plt.xlim(g_ra-.25, g_ra+.25)
    plt.ylim(g_dec-.25, g_dec+.25)
    #Add a grid to the plot
    plt.grid()
    #Invert the x axis
    plt.gca().invert_xaxis()
    #Add a legend
    plt.legend()
    #Save the figure to the current directory
    plt.savefig(f'source{index}scatterplots.png');
def create_points(table, index):
   """
   Returns point arrays from available observations for the identified source that give relavent
   information about the source. This is intended to be given to the clust_points function.
   
   Parameters
   ---
   
   table: DataFrame
          a pandas data frame generated by the match_to_gleam function
   
   index: integer or string 
          an integer index or string identifier for the target source
   
   Returns
   ---
   clust_ras: list
              list of lists. Each list is a series of Right Assention values in decimal degrees
              of the extended points for this source. Each list represents a seperate observation
   
   clust_decs: list
              list of lists. Each list is a series of Declination values in decimal degrees
              of the extended points for this source. Each list represents a seperate observation
   
   clust_mags: list
              list of lists. Each list is a series of flux values in  Jansky
              of the extended points for this source. Each list represents a seperate observation
              
   pointmags: list
              list of point-like flux of the source in Jansky 
   
   x_axis: numpy array
           Array representing the x-axis for the modeling area
   
   y_axis  numpy array
           Array representing the y-axis for the modeling area
   """
   #Set a vavriable to count the number of observations to zero
   obscount = 0
   #Create empty lists to store data
   immags = []
   imras = []
   imdecs = []
   upras = []
   updecs = []
   upmags = []
   ra_obs = []
   dec_obs = []
   mag_obs = []
   obs_data = []
   ra_obs_smth = []
   dec_obs_smth = []
   mag_obs_smth = []
   mag_change = []
   image_ra = []
   image_dec = []
   image_mag = []
   image_ston = []
   observation_number = []
   ps_mag = []
   upper_ra = []
   upper_dec = []
   h_angs = []
   otimes = []
   true_o_num = []
   adj_obs = []
   raz = []
   decz = []
   true_eo_mags = []
   #Identify the RA and DEC values of the original GLEAM objects
   gleam_ra = table.loc[index]['RA']
   gleam_dec = table.loc[index]['DEC']
   #mags are point source magnitudeds
   mags = table.loc[index][6::9]
   #mags_eo are the extended magnitude components
   mags_eo = table.loc[index][12::9]
   #ra for the extended components
   ras = np.array(table.loc[index][10::9])
   #dec for extended components
   decs = np.array(table.loc[index][11::9])
   #Singal to noise as a  point source value
   STON_ps = table.loc[index][13::9]
   #Grab the extended RA and DEC for the source
   high_ra = table.loc[index][8::9]
   high_dec = table.loc[index][9::9]
   #Create a variable to flag large sources
   large_source_flag = 0
   #Loop over the observations in the table
   for i in range(0, len(mags)):
       #Check if the observation identified a match
       if mags[i] !=0:
           #Add the point-like magnitude to a list
           mag_change.append(mags[i])
           #Increase the conter for the number of observations
           obscount+=1
           #Add the extended RA and DEC values to the list
           raz.append(high_ra[i])
           decz.append(high_dec[i])
           #Check if the source spills out of the 0.5 by 0.5 degree modeling area, flag if so
           if abs(np.max(ras[i]) - np.min(ras[i]))> 0.5:
               large_source_flag +=1
           if abs(np.max(decs[i]) - np.min(decs[i]))> 0.5:
               large_source_flag +=1
   #Create the x and y axes, depending on if the source requires a larger modeling area
   if large_source_flag == 0:
       x_stretch = .25
       y_stretch = .25
   else:
       x_stretch = .35
       y_stretch = .35
   #Create min and max x and y axis values
   xmin = gleam_ra - x_stretch
   xmax = gleam_ra + x_stretch
   ymin = gleam_dec - y_stretch
   ymax = gleam_dec + y_stretch
   #Create x and y axes 
   x_axis = np.linspace(xmin - 0.001, xmax + 0.001, num=100, endpoint=True)
   y_axis = np.linspace(ymin - 0.001, ymax + 0.001, num=100, endpoint=True)
   #Find the source center in degrees from the identified observations
   x_center = np.mean(raz)
   y_center = np.mean(decz)
   #Loop over the seperate obersvations
   for i in range(0, len(mags)):
       #Check if the observation is a match
       if mags[i] !=0:
           #Append the RA and DEC values to the apropriate list
           rs = ras[i]
           ds = decs[i]
           #Move the arrays towards the calculated center
           ra_roll = x_center - high_ra[i]
           dec_roll = y_center - high_dec[i]
           #Create an array with the adjusted observations
           new_array = rs+ra_roll, ds+dec_roll
           #Add array with RA and DEc adjustments to master list
           adj_obs.append(new_array)
           #Add the flux arrays to a list
           true_eo_mags.append(mags_eo[i])
   #Loop over the observations
   for n in range(0, len(mags)):
       #Check if the observation is identified as a match
       if mags[n] !=0:
           #Store the data for each match. 
           #the image arrays are arrays made of lists with one list for each match 
           image_ra.append(ras[n])
           image_dec.append(decs[n])
           image_mag.append(mags_eo[n])
           observation_number.append(n)
           image_ston.append(STON_ps[n])
           ps_mag.append(mags[n])
           upper_ra.append(high_ra[n])
           upper_dec.append(high_dec[n])
   #Calculate the mean image center
   center_ra = np.mean(upper_ra)
   center_dec = np.mean(upper_dec)
   #Create empty lists to store the data
   clust_ras = []
   clust_decs = []
   clust_mags = []
   pointmags  = []
   #Loop through the identified data and add relavent data to the apropriate list
   for m in range(0, len(adj_obs)):
       clust_ras.append(adj_obs[m][0])
       clust_decs.append(adj_obs[m][1])
       clust_mags.append(image_mag[m])
       pointmags.append(ps_mag[m])
   return clust_ras, clust_decs, clust_mags, pointmags, x_axis, y_axis
def source_search(table, n_points = None, min_mg = None, min_o = None, ra_min= None, ra_max = None, 
                  dec_min = None, dec_max = None, exclude_already_unique = True):
    """
    Searches for sources to be targeted for modeling
    
    Parameters
    ---
    table: DataFrame
          DataFrame generated by match_to_GLEAM function
    
    n_points: int, optional
                the minimum number of points in any observation of the source to flag it
    
    min_mg: float, optional
            the minimum brightness in Jy a source must meet to be identified
    
    min_o: int, optional
            minimum number of observations a source must meet to be identified. By default, the source
            must also be diffuse in a minimum of half of observations to be considered for modeling
    
    ra_min: float, optional
            minimum RA positional value in decimal degrees
            
    dec_min: float, optional
            minimum DEC positional value in decimal degrees
            
    ra_max: float, optional
            maximum RA positional value in decimal degrees
            
    dec_max: float, optional
            maximum DEC positional value in decimal degrees
            
    exclude_already_unique: bool, optional
                            When set to True, all identifiers that are not formated with the generic
                            'J' at the begining of their string identified will be excluded in searches.
                            By default, this is set to True
            
    Returns
    ---
    
    source_idents: lsit of strings or integers
                   list contiaining identifiers for those sources meeting the specefied criteria
    
    observation_numbers: list of integers
                         list contining the integer number of observatios that contributes to the sources
                         identified in source_idents. Has length matching source_idents
    
    is_big: list of integer 0s and 1s
            list containing either a 0 for sources that do not require a larger modeling area or a 1
            for sources that do require a larger modeling area. Has length matching source_idents
    
    
    in_area: integer
             the ammount of sources in the given area that passed all checks except for that it 
             was not considered diffuse. This can be helpful in determining what fraction of sources
             in a given area or under certain parameters are identified as diffuse.
    
    """
    #Grab the indices of the input DataFrame
    inds = table.index
    if exclude_already_unique == True:
        for i in inds:
            if list(i)[0][0] != 'J':
                inds.remove(i)
    #Create empty lists for returns
    is_big = []
    observation_numbers  = []
    source_idents = []
    #Create variables to be increased during iterations
    sourcenums = 0
    in_area = 0
    #Print an indication that the search is begining
    print('Begining Search')
    #Check if optional keywords have been filled, if not replace with default values
    num_points = n_points if n_points else 15
    min_mag = min_mg if min_mg else 0.01
    min_obs = min_o if min_o else 4
    has_position_limits = 1 if ra_min else 0
    #For cases where the search is not concentrated on one part of the sky
    if has_position_limits == 1:
        #Loop through the string index values
        for index in inds:
            sourcenums +=1
            extended_obs = 0
            observation = 0
            #Create a print statement that indicates the loop's position in the 
            #source list
            print(f'Considering source {sourcenums} of {len(inds)}', end = "\r")
            #Get the information about the source in this iteration from the DataFrame
            mags = table.loc[index][6::9]
            mags_eo = table.loc[index][12::9]
            ra_center = table.loc[index][8::9]
            dec_center = table.loc[index][9::9]
            ext_ras = table.loc[index][10::9]
            ext_decs = table.loc[index][11::9]
            ext_ras = tab.loc[index][10::9]
            ext_decs = tab.loc[index][11::9]
            #Loop through the observations of the source
            for i in range(0, len(mags)):
                #Confirm that this observatoin was identified as a match
                if mags[i] !=0:
                    #Increase the variable that tracks observation number
                    observation +=1
                    #Check if the source is above the minimum indicated flux
                    if (mags[i] > min_mag):
                        #Check if the ammount of points in the observation is greater than the cutoff
                        if (len(mags_eo[i]) > num_points):
                            #If this case is satisfied, increase the iteration tracker for 
                            #extended observations
                            extended_obs +=1
                            #Check if the source is already in the list of identified sources
                            if index not in source_idents:
                                already_in = 0
                            #Get the values for RA and DEC positional coordinates from the DataFrame
                            #for this source from the original input GLEAM catalog
                            ind_ra = table.loc[index]['RA']
                            ind_dec = table.loc[index]['DEC']
                            #Check if the RA value is greater than 180, subtract 360 if so to 
                            #center the data at 0
                            if ind_ra > 180:
                                ind_ra = ind_ra - 360
                            #Check that the RA and DEC values are within the specefied values
                            if ra_min <ind_ra < ra_max:
                                if dec_min<ind_dec< dec_max:
                                    #Loop through the sources already in the list of identified sources
                                    #to check if the source currently being iterated over overlaps with
                                    #any of these
                                    for source in source_idents:
                                        #Get the RA and DEC values of the source from the original GLEAM
                                        #input catalog
                                        source_ra = table.loc[source]['RA']
                                        source_dec = table.loc[source]['DEC']
                                        #Adjust the RA coordinate so that it is also centered on 0
                                        if source_ra > 180:
                                            source_ra = source_ra - 360
                                        #Check if the RA and DEC values fall within the range of the modeling
                                        #area for a source already in the list of sources
                                        if ((source_ra -0.25)< ind_ra <(source_ra+0.25)):
                                            if ((source_dec -0.25)< ind_dec <(source_dec+0.25)):
                                                #If so, increase the variable that identifies if a 
                                                #source is already modeled by 1
                                                already_in +=1
                                        #Check if the source spills out of the 0.5 degree by 0.5 degree
                                        #modeling area and check if the larger modeling areas clash with
                                        #already identified sources
                                        if abs(np.max(ext_ras[i]) - np.min(ext_ras[i]))>.5:
                                            if ((source_ra -0.35)< ind_ra <(source_ra+0.35)):
                                                if ((source_dec -0.5)< ind_dec <(source_dec+0.35)):
                                                    already_in +=1
                                        elif abs(np.max(ext_decs[i]) - np.min(ext_decs[i]))>.5:
                                            if ((source_ra -0.35)< ind_ra <(source_ra+0.35)):
                                                if ((source_dec -0.35)< ind_dec <(source_dec+0.35)):
                                                    already_in +=1
                                    #If the source has passed these checks and is not already covered by
                                    #another source, proceed to other checks
                                    if already_in == 0:
                                        #Check if the source appears in the minimum number of observations
                                        if observation >= min_obs:
                                            #If so, increase the tracker for the number of sources in 
                                            #the area that satisfy this requirement
                                            in_area+=1
                                            #Check if more than half of the observations are extended
                                            if extended_obs > (0.5*observation):
                                                #If all these and above are satisfied, add the index 
                                                #to the master list of identified sources and the integer
                                                #number of observations to the tracking list
                                                source_idents.append(index) 
                                                observation_numbers.append(observation)
                                            #Check if the source spills out of the modeling area and add a 1
                                            #to the is_big list if this is true
                                            if abs(np.max(ext_decs[i]) - np.min(ext_decs[i]))>.5:
                                                if abs(np.max(ext_ras[i]) - np.min(ext_ras[i]))>.5:
                                                    is_big.append(1)
                                            #If the source does not spill out of the area, add a 0
                                            #to the is_big list
                                            else:
                                                is_big.append(0)
    #For cases where there is no specification for RA and DEC position limitations                                        
    else:
        #Loop through the string index values
        for index in inds:
            sourcenums +=1
            extended_obs = 0
            observation = 0
            #Create a print statement that indicates the loop's position in the 
            #source list
            print(f'Considering source {sourcenums} of {len(inds)}', end = "\r")
            #Get the information about the source in this iteration from the DataFrame
            mags = table.loc[index][6::9]
            mags_eo = table.loc[index][12::9]
            ra_center = table.loc[index][8::9]
            dec_center = table.loc[index][9::9]
            ext_ras = table.loc[index][10::9]
            ext_decs = table.loc[index][11::9]
            ext_ras = tab.loc[index][10::9]
            ext_decs = tab.loc[index][11::9]
            #Loop through the observations of the source
            for i in range(0, len(mags)):
                #Confirm that this observatoin was identified as a match
                if mags[i] !=0:
                    #Increase the variable that tracks observation number
                    observation +=1
                    #Check if the source is above the minimum indicated flux
                    if (mags[i] > min_mag):
                        #Check if the ammount of points in the observation is greater than the cutoff
                        if (len(mags_eo[i]) > num_points):
                            #If this case is satisfied, increase the iteration tracker for 
                            #extended observations
                            extended_obs +=1
                            #Check if the source is already in the list of identified sources
                            if index not in source_idents:
                                already_in = 0
                            #Get the values for RA and DEC positional coordinates from the DataFrame
                            #for this source from the original input GLEAM catalog
                            ind_ra = table.loc[index]['RA']
                            ind_dec = table.loc[index]['DEC']
                            #Loop through the sources already in the list of identified sources
                            #to check if the source currently being iterated over overlaps with
                            #any of these
                            for source in source_idents:
                                #Get the RA and DEC values of the source from the original GLEAM
                                #input catalog
                                source_ra = table.loc[source]['RA']
                                source_dec = table.loc[source]['DEC']
                                #Check if the RA and DEC values fall within the range of the modeling
                                #area for a source already in the list of sources
                                if ((source_ra -0.25)< ind_ra <(source_ra+0.25)):
                                    if ((source_dec -0.25)< ind_dec <(source_dec+0.25)):
                                        #If so, increase the variable that identifies if a 
                                        #source is already modeled by 1
                                        already_in +=1
                                #Check if the source spills out of the 0.5 degree by 0.5 degree
                                #modeling area and check if the larger modeling areas clash with
                                #already identified sources
                                if abs(np.max(ext_ras[i]) - np.min(ext_ras[i]))>.5:
                                    if ((source_ra -0.35)< ind_ra <(source_ra+0.35)):
                                        if ((source_dec -0.5)< ind_dec <(source_dec+0.35)):
                                            already_in +=1
                                elif abs(np.max(ext_decs[i]) - np.min(ext_decs[i]))>.5:
                                    if ((source_ra -0.35)< ind_ra <(source_ra+0.35)):
                                        if ((source_dec -0.35)< ind_dec <(source_dec+0.35)):
                                            already_in +=1
                            #If the source has passed these checks and is not already covered by
                            #another source, proceed to other checks
                            if already_in == 0:
                                #Check if the source appears in the minimum number of observations
                                if observation >= min_obs:
                                    in_area +=1
                                    #Check if more than half of the observations are extended
                                    if extended_obs > (0.5*observation):
                                        #If all these and above are satisfied, add the index 
                                        #to the master list of identified sources and the integer
                                        #number of observations to the tracking list
                                        source_idents.append(index) 
                                        observation_numbers.append(observation)
                                    #Check if the source spills out of the modeling area and add a 1
                                    #to the is_big list if this is true
                                    if abs(np.max(ext_decs[i]) - np.min(ext_decs[i]))>.5:
                                        if abs(np.max(ext_ras[i]) - np.min(ext_ras[i]))>.5:
                                            is_big.append(1)
                                    #If the source does not spill out of the area, add a 0
                                    #to the is_big list
                                    else:
                                        is_big.append(0)

    print('Finished')
    return source_idents, observation_numbers, is_big, in_area
def clust_points(clust_ras, clust_decs, clust_mags, pointmags, radius, percent_total):
    """
    Creates clustered source models given specified input data. This function is meant to be 
    called by the create_plots function.
    
    This function clusters points to identify portions of an object that have potential to be point-like. 
    Each observation of a source is considered individualy. The function will loop through all points within 
    an obesrvation and place a circle around these points of the speciified radius. Other points in the 
    observation that fit within this radius are considered together. If this cluster of points has a 
    total flux in Jy greater than some percent of the total source flux, specified by percent_total, all 
    points within this radius are clustered together into a signle point and removed from the loop.
    
    Parameters
    ---
    
    clust_ras: list of lists. Each list is a series of Right Assention values in decimal degrees
               of the extended points for this source. Each list represents a seperate observation
    
    clust_decs:list of lists. Each list is a series of Declination values in decimal degrees
               of the extended points for this source. Each list represents a seperate observation
    
    clust_mags:list of lists. Each list is a series of flux values in  Jansky
               of the extended points for this source. Each list represents a seperate observation
               
    pointmags: list of point-like flux of the source in Jansky 
    
    x_axis: Array representing the x-axis for the modeling area
    
    y_axis  Array representing the y-axis for the modeling area
    
    
    radius: the radius of the circle that will determine if objects are point like in decimal degrees
    
    percent_total: the percent cutoff, in decimal percentage, of the total object brightness 
    that will determine if an object is point-like
    
    Returns
    ---
    pointlike_dict: a dictionary containing information about the clustered points. The keys have the 
                    following format
                    
                    obsx:The highest-level keys have format obsx where x is an integer indicating the 
                         observation numer to which the points within that key belong
                    
                    component_number_x: These keys are contained within the observation keywords. x indicates 
                                        the component number associated with this point within the observation
                                        and is an integer
                                        
                    The remaining keywords are a dictionary stored wihtin the component_number_x keys
                    
                    sum_match: the total flux in Jy that went into the particular component
                    
                    component_ras: the RA positions in decimal degrees of the points that atributed to
                                   the component 
                                   
                    component_decs: the DEC positions in decimal degrees of the points that atributed to
                                   the component 
                                   
                    component_mags: the Mag positions in Jy of the points that atributed to
                                   the component 
                    
                    component_ind: the integer indices of the points that atributed to
                                   the component 
                                   
                    weight_ra: the weighted position in decimal degrees of the RA for the component
                    
                    weight_dec: the weighted position in decimal degrees of the DEC for the component
                    
    untouched_ra: the original input RA values in decimal degrees
    
    untouched_dec: the original input DEC values in decimal degrees
    
    untouched_mag: the original input flux values in Jy
    """
    #Creating several empty lists
    onums = []
    untouched_ra = []
    untouched_dec = []
    untouched_mag = []
    colorlims = []
    obsdatasmth = []
    untouched_ra = []
    untouched_dec = []
    untouched_mag = []
    #Try
    try:
        #Create an empty dictionary
        pointlike_dict = {}
        #Loop over the number of total observations of the source identified
        #by the previously called functions
        for m in range(0, len(clust_ras)):
            #Add an integer to indicate the number of observations
            onums.append(m)
            #Identify the RA, DEC, and flux values from the input data
            ra = clust_ras[m]
            dec = clust_decs[m]
            mag = clust_mags[m]
            #Identify the high-level point-like flux for the source
            point_mag = pointmags[m]
            #Save these input arrays so that the function can output the original
            #input point clusters 
            untouched_ra.append(ra)
            untouched_dec.append(dec)
            untouched_mag.append(mag)
            #Identify the point-like flux value
            upper_mag = pointmags[m]
            #Create empty lists to store the data
            match_ra = []
            match_dec = []
            match_mag = []
            match_ind = []
            smth_ra = []
            smth_dec = []
            smth_mag = []
            m_ras = []
            m_decs = []
            m_mags = []
            s_ras = []
            s_decs = []
            s_mags = []
            idx = []
            #Create two variables that can be adjusted to track iterations
            n = 0
            iterations = 0
            #Add the first key to the dictionary, indicated by the observation number
            pointlike_dict[f'obs{m}'] = {}
            #Try here will catch all cases except those that are point-like observations
            try:
                #Loop through individual points in this observation
                for r in range(0, len(ra)):
                    #Add the flux value to a list for tracking the flux limits
                    colorlims.append(mag[r][0])
                    #Check if this point is in the list match_ra
                    if ra[r] in match_ra:
                        #If it is, do nothing.
                        s = 0
                    #If the point does not already exist in the match_ra list, we must attepmt to
                    #cluster it
                    else:
                        #Increase n to track the number of iterations for clustering
                        n +=1
                        #Define the first coordinate to be considered
                        #x and y are the positions in decimal degrees
                        x = ra[r]
                        y = dec[r]
                        #power is the flux of the point
                        power = mag[r][0]
                        #xmatch, ymatch, and magmatch are empty lists to track matches
                        xmatch = []
                        ymatch = []
                        magmatch = []
                        #idxq tracks the indices of these points
                        idxq = []
                        #iterations is increased by one to track
                        iterations = iterations + 1
                        #We next loop through all other points in the observation
                        for q in range (0, len(ra)):
                            #Check if this point has already been pulled out by a previous iteration of the loop
                            if q in match_ind:
                                #If it has been identified, do nothing.
                                s = 0
                            #If it has not been identified, the point is considered
                            else:
                                #Calculate the distance between the x,y point defined above and this particular 
                                #point in the rest of the data
                                dist = np.sqrt((((ra[q]) - x)**2) + (((dec[q]) - y)**2))
                                #Check if the distance is less than the clustering radius
                                if dist < radius:
                                    #Add information about points that pass this check to the relavent lists
                                    #idxq tracks the index of these matches
                                    idxq.append(q)
                                    #xmatch and ymatch track the RA and DEC values in decimal degrees
                                    xmatch.append(ra[q])
                                    ymatch.append(dec[q])
                                    #magmatch tracks the flux value in Jy 
                                    magmatch.append(mag[q][0])
                            #Check if the sum of the flux of all points identified in the previous loop 
                            #are greater than the percentage cutoff indicated at function call to indicate
                            #that the cluster should be pulled out
                            if (np.sum(magmatch) > (percent_total*(np.sum(upper_mag)))):
                                #If these points will be pulled from the list, their indices as tracked by
                                #the idxq array are added to the master index array idx
                                for num in idxq:
                                    idx.append(num)
                                #The position coordinates, flux arrays, and indices are added to
                                #lists specifically tracked by this loop iteration
                                match_ra.extend(xmatch)
                                match_dec.extend(ymatch)
                                match_mag.extend(magmatch)
                                match_ind.extend(idxq)
                                #The mean position in RA is calculated, weighted by the point flux values
                                mean_pos_ra = np.average(xmatch, weights = magmatch)
                                #The mean position in DEC is calculated, weighted by the point flux values
                                mean_pos_dec = np.average(ymatch, weights = magmatch)
                                #A keyword within the key for this observation is created for this point.
                                #The point dictionary has information about the total flux clustered in
                                #this point as well as the position and flux of the components that went
                                #into it and the indices of these points in the larger array. The mean
                                #weighted position is also indicated.
                                pointlike_dict[f'obs{m}'][f'component_number_{r}'] = {'sum_match': np.sum(magmatch), 
                                                                            'component_ras' : xmatch, 
                                                                            'component_decs' : ymatch,
                                                                            'component_mags' : magmatch,
                                                                            'component_ind':idxq, 
                                                                            'weight_ra':mean_pos_ra,
                                                                            'weight_dec': mean_pos_dec}
            #Except case here is for when the observation has a single point, and so nothing will be 
            #clustered but the single point is still added to the master dictionary for this source
            except:
                #This will be the only component in this observation as it is point-like.
                #The information about the point is directly stored with 0s in the 
                #component_x keywords indicating that no points went into this cluster
                r = 0
                pr = ra
                pd = dec
                power = mag[0]
                pointlike_dict[f'obs{m}'][f'component_number_{r}'] = {'sum_match': power, 
                                                            'component_ras' : [0],
                                                            'component_decs' : [0],
                                                            'component_mags' : [0],
                                                            'component_ind':[0], 'weight_ra':pr,
                                                            'weight_dec': pd}
            #Try here catched all cases where the observation is not already point-like
            try:
                #create an empty list for indices to be stored
                match_idx = []
                #Loop through the components in this observation's dictionary
                for key in pointlike_dict[f'obs{m}'].keys():
                    #Extend the list with the indices of the components that went in to the point
                    match_idx.extend(pointlike_dict[f'obs{m}'][key]['component_ind'])
                #Create a master variable with all indices of all points that went in to point-like
                #structure for the observation
                n = list(set(match_idx))
                #Print an indication of how many loop iterations ran, how many points from the original
                #observation were extracted, and the total original points
                print('Total Points', len(ra), 'Number Extracted', len(ra[n]), end = "\r")
                print( 'Total Iterations', np.sum(iterations), end = "\r")
            #Except catched the case where the observation was already a single point
            except:
                #Print an indication that the observation was not iterated over and is point-like
                print('Point Like Observation', end = "\r")
    #Except catches            
    except:
        #Create a variable m that can be used to track total observations
        m = 0
        #Define an empty dictionary to store information about clustering
        pointlike_dict = {}
        #Add the current integer observation number to the number of observations
        onums.append(m)
        #Identify RA, DEC, and flux of the clusters of points within this observation
        ra = clust_ras
        dec = clust_decs
        mag = clust_mags
        #Identify the high-level point-like flux value in Jy
        point_mag = pointmags
        #Story the untouched RA, DEC, and Flux values
        untouched_ra.append(ra)
        untouched_dec.append(dec)
        untouched_mag.append(mag)
        #Store high-level point-lke flux
        upper_mag = pointmags
        #Define empty lists for storing the data
        match_ra = []
        match_dec = []
        match_mag = []
        match_ind = []
        smth_ra = []
        smth_dec = []
        smth_mag = []
        m_ras = []
        m_decs = []
        m_mags = []
        s_ras = []
        s_decs = []
        s_mags = []
        idx = []
        #Create variables to track iterations
        n = 0
        iterations = 0
        #Create a keywords in the dictionary for this observation
        pointlike_dict[f'obs{m}'] = {}
        #Try here will catch all cases except those that are point-like observations
        try:
            #Loop through individual points in this observation
            for r in range(0, len(ra)):
                #Check if this point is in the list match_ra
                if ra[r] in match_ra:
                    #If it is, do nothing
                    s = 0
                #If it has not been identified, the point is considered
                else:
                    #Increase n to track iterations for clusterig
                    n +=1
                    #Define the first coordinate to be considered
                    #x and y are the positions in decimal degrees
                    x = ra[r]
                    y = dec[r]
                    #power is the flux of the point
                    power = mag[r]
                    #xmatch, ymatch, and magmatch are empty lists to track matches
                    xmatch = []
                    ymatch = []
                    magmatch = []
                    #idxq tracks the indices of these points
                    idxq = []
                    #iterations is increased by one to track
                    iterations = iterations + 1
                    #We next loop through all other points in the observation
                    for q in range (0, len(ra)):
                        #Check if this point has already been pulled out by a previous iteration of the loop
                        if q in match_ind:
                            #If it has been identified, do nothing.
                            s = 0
                        else:
                            #Calculate the distance between the x,y point defined above and this particular 
                            #point in the rest of the data
                            dist = np.sqrt((((ra[q]) - x)**2) + (((dec[q]) - y)**2))
                            #Check if the distance is less than the clustering radius
                            if dist < radius:
                                #Add information about points that pass this check to the relavent lists
                                #idxq tracks the index of these matches
                                idxq.append(q)
                                #xmatch and ymatch track the RA and DEC values in decimal degrees
                                xmatch.append(ra[q])
                                ymatch.append(dec[q])
                                magmatch.append(mag[q])
                        #Check if the sum of the flux of all points identified in the previous loop 
                        #are greater than the percentage cutoff indicated at function call to indicate
                        #that the cluster should be pulled out
                        if (np.sum(magmatch) > (percent_total*(np.sum(upper_mag)))):
                            #If these points will be pulled from the list, their indices as tracked by
                            #the idxq array are added to the master index array idx
                            for num in idxq:
                                idx.append(num)
                            #The position coordinates, flux arrays, and indices are added to
                            #lists specifically tracked by this loop iteration
                            match_ra.extend(xmatch)
                            match_dec.extend(ymatch)
                            match_mag.extend(magmatch)
                            match_ind.extend(idxq)
                            #The mean position in RA is calculated, weighted by the point flux values
                            mean_pos_ra = np.average(xmatch, weights = magmatch)
                            #The mean position in DEC is calculated, weighted by the point flux values
                            mean_pos_dec = np.average(ymatch, weights = magmatch)
                            #A keyword within the key for this observation is created for this point.
                            #The point dictionary has information about the total flux clustered in
                            #this point as well as the position and flux of the components that went
                            #into it and the indices of these points in the larger array. The mean
                            #weighted position is also indicated.
                            pointlike_dict[f'obs{m}'][f'component_number_{r}'] = {'sum_match': np.sum(magmatch), 
                                                                        'component_ras' : xmatch, 
                                                                        'component_decs' : ymatch,
                                                                        'component_mags' : magmatch,
                                                                        'component_ind':idxq, 
                                                                        'weight_ra':mean_pos_ra,
                                                                        'weight_dec': mean_pos_dec}
        #Except case here is for when the observation has a single point, and so nothing will be 
        #clustered but the single point is still added to the master dictionary for this source
        except:
            #This will be the only component in this observation as it is point-like.
            #The information about the point is directly stored with 0s in the 
            #component_x keywords indicating that no points went into this cluster
            r = 0    
            pr = ra
            pd = dec
            power = mag[0]
            pointlike_dict[f'obs{m}'][f'component_number_{r}'] = {'sum_match': power, 
                                                        'component_ras' : [0],
                                                        'component_decs' : [0],
                                                        'component_mags' : [0],
                                                        'component_ind':[0], 'weight_ra':pr,
                                                        'weight_dec': pd}
        #Try here catched all cases where the observation is not already point-like
        try:    
            #create an empty list for indices to be stored
            match_idx = []
            #Loop through the components in this observation's dictionary
            for key in pointlike_dict[f'obs{m}'].keys():
                #Extend the list with the indices of the components that went in to the point
                match_idx.extend(pointlike_dict[f'obs{m}'][key]['component_ind'])
            #Create a master variable with all indices of all points that went in to point-like
            #structure for the observation
            n = list(set(match_idx))
            #Print an indication of how many loop iterations ran, how many points from the original
            #observation were extracted, and the total original points
            print('Total Points', len(ra), 'Number Extracted', len(ra[n]), end = "\r")
            print( 'Total Iterations', np.sum(iterations), end = "\r")
            #Except catched the case where the observation was already a single point
        except:
            #Print an indication that the observation was not iterated over and is point-like
            print('Point Like Observation', end = "\r")
    return pointlike_dict, untouched_ra, untouched_dec, untouched_mag
def create_plots(table, index):
    """
    Generates model with relavent plot for indicated source
    
    Parameters
    ---
    
    table: pandas DataFrame generated by the match_to_GLEAM function
    
    index: integer or string index identifying the source a model is to be generated for 
    
    Returns
    ---
    
    """
    #Formatting Matplotlib
    plt.rcParams['figure.figsize'] = (10, 10)
    plt.rc('axes', labelsize=14)
    plt.rc('axes', labelweight='bold')
    plt.rc('axes', titlesize=16)
    plt.rc('axes', titleweight='bold')
    plt.rc('font', family='sans-serif')
    #Seting the radius for the circle to determine what is pointlike
    radius = 0.000277778*20
    #Determining the percentage of total source brightness used to determine what is pointlike
    percent_total = .15
    #Creating the necessary point arrays and x and y axis dimmensions for this source
    cx, cy, cz, cp, x_axis, y_axis = create_points(table, index)
    #Getting the pointlike sources and untouched point arrays
    pointlike_dict, untouched_ra, untouched_dec, untouched_mag = clust_points(cx, cy, cz, cp, radius, percent_total)
    #Creating a check to determine what was pulled out
    point_check = 0
    #Check for source type and indicate if source falls outside of print statements implemented by
    #previously called functions
    for c in cx:
        if type(c) != np.ndarray:
            point_check +=1
    if pointlike_dict == {}:
        return print('Source Not Observed', end = "\r")
    elif point_check == len(cx):
        return print('Source is Observed to be Pointlike', end = "\r")
    else:
        #Defining some empty lists to store observation data
        obs_data = []
        upper_obs_data = []
        #Create two variables that can be increased with the loops
        n = -1
        numobs = 0
        # Going through each observation of the source
        for dat in range(len(pointlike_dict)):
            #Defining empty lists to store information
            upper_ra_points = []
            upper_dec_points = []
            upper_mag_points = []
            contribra = []
            contribdec = []
            contribmag = []
            imgras = []
            imgdecs = []
            imgmags = []
            n = n+1
            match_idx = []
            #Open each of the individual components of the observation that was pulled into the 
            #point-like dictionaries generated by the clust_points function
            for key in pointlike_dict[f'obs{n}'].keys():
                upper_ra_points.append(pointlike_dict[f'obs{n}'][key]['weight_ra'])
                upper_dec_points.append(pointlike_dict[f'obs{n}'][key]['weight_dec'])
                upper_mag_points.append(pointlike_dict[f'obs{n}'][key]['sum_match'])
                contribra.extend(pointlike_dict[f'obs{n}'][key]['component_ras'])
                contribdec.extend(pointlike_dict[f'obs{n}'][key]['component_decs'])
                contribmag.extend(pointlike_dict[f'obs{n}'][key]['component_mags'])
            #Try here catches cases where there are points recorded for the observation 
            try:
                #loop through the untouched data
                for coord in range(len(untouched_ra[dat])):
                    #check if the untouched component contributed to the point sources and add to
                    # the data arrays if they did not
                    if untouched_ra[dat][coord] not in np.array(contribra):
                        imgras.append(untouched_ra[dat][coord])
                        imgdecs.append(untouched_dec[dat][coord])
                        imgmags.append(untouched_mag[dat][coord])
                #Create a dictionary containing the weighted clustered points
                pointobs = {'ra': upper_ra_points, 'dec': upper_dec_points, 'mag': upper_mag_points}
                #append the dictionary to a list of dictionaries
                upper_obs_data.append(pointobs)
                #increase the number of observations
                numobs += 1
                #Create an empty 100x100 matrix to store the data
                data = np.matrix(np.zeros((100,100)))
                #loop through the points that were ot pulled into clusters
                for val in range (0, len(imgdecs)):
                        #Identify the point positional values and flux
                        smooth_x = imgras[val]
                        smooth_y = imgdecs[val]
                        smooth_m = imgmags[val]
                        #Enumerate the x and y axes, allowing us to loop through all the points in the data 
                        #matrix
                        for xind, xval in enumerate(x_axis):
                            for yind, yval in enumerate(y_axis):
                                #Set the value for the data matrix at this point to include the value of the 
                                #gaussian beam created by the point we are currently looping over 
                                data[xind, yind] += gaussian(1.5*(x_axis[1]-x_axis[0]), smooth_x, smooth_y, xval, yval, smooth_m)
            #Except here catches the case when there are no components to smooth out and returns an empty 
            #data matrix
            except:
                data = np.matrix(np.zeros((100,100)))
            #Append the data matrix, whether it be blank or containing smoothed gaussian surfaces,
            #to a list containing the data arrays from the observations
            obs_data.append(data)
            #Create a matplotlib figure
            fig, ax = plt.subplots(figsize = (10,10))
            #Plot the data matrix generated by the most recent observation
            implot = plt.imshow(np.flip(data.T, axis = 0), cmap = 'pink', interpolation = None, origin = 'upper', 
                               extent = [np.min(x_axis), np.max(x_axis), np.min(y_axis), np.max(y_axis)]) #, origin = 'lower', interpolation = None)#'gist_ncar')
            #Add a colorbar
            plt.colorbar()
            #Create a scatter plot over this image containing all points that have been clustered out of
            #the observation
            plt.scatter(upper_ra_points, upper_dec_points, color = 'gold', s = upper_mag_points)
            #Invert the RA values
            plt.gca().invert_xaxis()
            #Create a title for the observation plot indicating the observation number
            plt.title(f'Observation {dat}, First Iteration')
            #Add x and y labels
            plt.xlabel('RA (Deg)')
            plt.ylabel('DEC (Deg)')
            #Create a unique file name containing the source index, observation number, clustering radius,
            #and flux percentage cut-off. The image is displayed and saved to the current working directory
            plt.savefig(f'source{index}observation{dat}radius{radius}bright{percent_total}scatter.png')
            #Display figure
            plt.show()
            #Create a fresh figure
            plt.figure()
        #Grag the values from the data frame
        mags = table.loc[index][6::9]#upper level magnitude values
        ra_center = table.loc[index][8::9]#upper level ra values
        dec_center = table.loc[index][9::9]#upper level dec values
        ext_ras = table.loc[index][10::9] #lower level RA values
        ext_decs = table.loc[index][11::9]#lower level DEC values
        #Find the dimmensions of the square pixel
        pixel_dim = x_axis[1]-x_axis[0]
        #create some empty lists to store data
        shifts = []
        newdata = []
        #Create a stacking variable and innitiate it as None
        imstack = None 
        #Loop through the data matricies
        for img in range(len(obs_data)):
            #Grab the observation
            shifted_obs = np.expand_dims(obs_data[img], axis = 0)
            #Check if this is the first iteration
            if imstack is None:
                #If it is the first iteration, simply add this to the data stack
                imstack = shifted_obs
            #If this is not the first iteration, stack it with the others using vstack
            else:
                imstack = np.vstack((imstack, shifted_obs))
        #Create a combination of the images by taking the median pixel value of the observations
        #This adds some accounting for outliers without requireing a complex analysis of the variables
        #that can contribute to differences between observations
        raw_smoothed_image = np.median(imstack, axis = 0)
        #Create a new figure
        fig, ax = plt.subplots(figsize = (10,10))
        #Plot the image with the median pixel generated by the stacked observations
        #CHECK IN: see if division by observation number is necessary, prove that it returns propper fux vals
        implot = plt.imshow(np.flip((raw_smoothed_image/numobs).T, axis = 0), cmap = 'pink', interpolation = None, origin = 'upper', 
                            extent = [np.min(x_axis), np.max(x_axis), np.min(y_axis), np.max(y_axis)]) 
        #Add a colorbar
        plt.colorbar()
        #Loop through the point-like structure pulled out in each individual observation, plot each
        #observation in a different color and indicate this in the legend
        for scat in range (len(upper_obs_data)):
            plt.scatter(upper_obs_data[scat]['ra'], upper_obs_data[scat]['dec'], s = upper_obs_data[scat]['mag'], 
                        label = f'Observation {scat}')
        #Display legend
        plt.legend()
        #Invert RA values
        plt.gca().invert_xaxis()
        plt.title(f'Median Observation and Point-like structure from individual observations')
        #Add x and y labels
        plt.xlabel('RA (Deg)')
        plt.ylabel('DEC (Deg)')
        #Save the figure to the current working directory and create a unique filename based on the
        #index identifyer, clustering radius, and flux percentage cut-off
        plt.savefig(f'source{index}medianobsradius{radius}bright{percent_total}scatternoclust.png')
        #Display the figure
        plt.show()
        #MASTER CLUSTERING PULLED POINTS
        #Create empty lists to store the data
        mastr_ras = []
        mastr_decs = []
        mastr_mags = []
        #Loop through the pulled point-like structure for the observations
        for k in range(0, len(upper_obs_data)):
            #Identify the RA, DEC, and flux values for these points, put them all together into
            #a master list
            mastr_ras.extend(upper_obs_data[k]['ra'])
            mastr_decs.extend(upper_obs_data[k]['dec'])
            mastr_mags.extend(upper_obs_data[k]['mag'])
        #Define a radius and flux percentage cutoff for this clustering
        radius = 0.000277778*30
        percent_total = .15
        #Call the clust_points function to cluster the point-like structure from all observations
        clustdict = clust_points(mastr_ras, mastr_decs, mastr_mags, np.mean(cp), radius, percent_total)
        #Create empty lists to store information about what points need to be removed
        true_keys = []
        false_keys = []
        #Define an empty dictionary to store information about the clustered points
        obscontriblen = {}
        #Try catches all cases except those where no points were clustered
        try:
            #Loop through the components within the dictionary. Only one observation key is generated
            #as the points were fed into the function as if they were a single observation
            for key in clustdict[0]['obs0'].keys():
                #Defing an empty list to store the infomation for this component
                contribobsclust = []
                #Loop through the obseration data
                for o in range(0, len(upper_obs_data)):
                    #Pull out the RA and DEC for this observation
                    observation_ra = upper_obs_data[o]['ra']
                    observation_dec = upper_obs_data[o]['dec']
                    #Loop through the ammount of individual components identified in this obesrvation
                    for indv in range(len(clustdict[0]['obs0'][key]['component_ras'])):
                        #Check if the component is in the untouched RA
                        if clustdict[0]['obs0'][key]['component_ras'][indv] in observation_ra:
                            #Check if the component is also in the untouched DEC
                            if clustdict[0]['obs0'][key]['component_decs'][indv] in observation_dec:
                                #if so, append the integer loop number
                                contribobsclust.append(o)
                #If the point has contributions from more than half of the observations of the source
                #it is kept by adding the key to the "true" key list. If not, the key is added to the 
                #false list
                if len(set(contribobsclust)) > (len(upper_obs_data) * 0.5):
                    true_keys.append(key)
                    obscontriblen[key] = len(set(contribobsclust))
                else:
                    false_keys.append(key)
                    obscontriblen[key] = len(set(contribobsclust))
        #Except catches cases where the observation has no contributing clustered points
        except:
            contribobsclust = []
        #This next step in the code identifies points that were pulled out of their original observations
        #but were clustered into a point that did not pass the checks to be put into the final model.
        #These points must be thrown back into the models for their original observations so that they
        #are still considered in the final model, just not as part of a point source.
        #Create empty lists to store the data for the points that need to be thrown back
        backras = []
        backdecs = []
        backmags = []
        #Loop through the observations 
        for obz in range(len(upper_obs_data)):
            #Create empty lists to store the information for this observation
            obsbackra = []
            obsbackdec = []
            obsbackmag = []
            #Create an integer that can be increased to track iterations
            numbertest = 0
            #Loop through the clustered points in this observaiton
            for key in pointlike_dict[f'obs{obz}'].keys():
                #Try here catches cases where multiple clustered points were identified for the observation
                try:
                    #Loop through the points
                    for upperpoint in range(0, len(pointlike_dict[f'obs{obz}'][key]['weight_ra'])):
                            #Grab the weighted RA and DEC positions from the dictionary
                            tstra = pointlike_dict[f'obs{obz}'][key]['weight_ra'][upperpoint]
                            tstdec = pointlike_dict[f'obs{obz}'][key]['weight_dec'][upperpoint]
                            #Loop through the keys identified as belonging to points that do not pass the 
                            #checks and need to be thrown back
                            for k in false_keys:
                                #Get the RA and DECs of the points that contributed to the clustered points
                                #for the key we are iterating over
                                cras = clustdict[0]['obs0'][k]['component_ras']
                                cdecs = clustdict[0]['obs0'][k]['component_decs'] 
                                #Check if the RA and DEC appear in the cluster
                                if tstra in cras:
                                    if tstdec in cdecs:
                                        #If they do, extend the original dictionary keys to include the 
                                        #previously tossed out points
                                        numbertest+=1
                                        obsbackra.extend(pointlike_dict[f'obs{obz}'][key]['component_ras'])
                                        obsbackdec.extend(pointlike_dict[f'obs{obz}'][key]['component_decs'])
                                        obsbackmag.extend(pointlike_dict[f'obs{obz}'][key]['component_mags'])
                #Except catches all cases where singular or no points were captured in the clustering
                except:
                    #Set the RA and DEC coordinates for the single point if it exists
                    #upperpoint = 0
                    tstra = pointlike_dict[f'obs{obz}'][key]['weight_ra']
                    tstdec = pointlike_dict[f'obs{obz}'][key]['weight_dec']
                    #Loop through the keys for points that need to be thrown back and check if
                    #any match this coordinate
                    for k in false_keys:
                        cras = clustdict[0]['obs0'][k]['component_ras']
                        cdecs = clustdict[0]['obs0'][k]['component_decs'] 
                        if tstra in cras:
                            if tstdec in cdecs:
                                #Increase the iteration tracker and extend the original observation arrays
                                #if they pass the checks
                                numbertest+=1
                                obsbackra.extend(pointlike_dict[f'obs{obz}'][key]['component_ras'])
                                obsbackdec.extend(pointlike_dict[f'obs{obz}'][key]['component_decs'])
                                obsbackmag.extend(pointlike_dict[f'obs{obz}'][key]['component_mags'])
            #Check if there are points that need to be thrown back into the original observations
            if np.sum(numbertest) != 0:
                #Throw these data points into the apropriate arrays
                backras.append(obsbackra)
                backdecs.append(obsbackdec)
                backmags.append(obsbackmag)
            else:
                #Add 0s if the observation has nothing to be thrown back
                backras.append(0)
                backdecs.append(0)
                backmags.append(0)
        #Create 2 iteration trackers and create an empty list to store the new observation data        
        obs_data = []
        n = -1
        numobs = 0
        #Loop through the observations in the dictionary
        for dat in range(len(pointlike_dict)):
            #Create empty lists to store the modified observation data
            upper_ra_points = []
            upper_dec_points = []
            upper_mag_points = []
            contribra = []
            contribdec = []
            contribmag = []
            imgras = []
            imgdecs = []
            imgmags = []
            #Increase the iteration tracker
            n = n+1
            match_idx = []
            #Loop through the clustered points in the observations
            for key in pointlike_dict[f'obs{n}'].keys():
                #Add the data from the components to the master lists
                upper_ra_points.append(pointlike_dict[f'obs{n}'][key]['weight_ra'])
                upper_dec_points.append(pointlike_dict[f'obs{n}'][key]['weight_dec'])
                upper_mag_points.append(pointlike_dict[f'obs{n}'][key]['sum_match'])
                contribra.extend(pointlike_dict[f'obs{n}'][key]['component_ras'])
                contribdec.extend(pointlike_dict[f'obs{n}'][key]['component_decs'])
                contribmag.extend(pointlike_dict[f'obs{n}'][key]['component_mags'])
            #Try catches cases where there are points that need to be added to original
            #observations and extends the observation data apropriately
            try:
                if np.sum(backras[dat]) != 0:
                    imgras.extend(backras[dat])
                    imgdecs.extend(backdecs[dat])
                    imgmags.extend(backmags[dat])
            #Except catches cases where there is nothing to extend and adds a 0 to indicate nothing
            #was thrown back from these observations
            except:
                imgras.append(0)
                imgdecs.append(0)
                imgmags.append(0)
            #Try catches cases where there exists data that needs to be "smoothed" and the 
            #observation is not point-like
            try:
                #Loop through the coordinates in the RA and DEC coordinates
                for coord in range(len(untouched_ra[dat])):
                    #Check if the coordinate is already in the matrix
                    #CHECK delete and: if error
                    if (untouched_ra[dat][coord] not in np.array(contribra)) and (untouched_dec[dat][coord] not in np.array(contribdec)):
                        #If so, add the coordinates
                        imgras.append(untouched_ra[dat][coord])
                        imgdecs.append(untouched_dec[dat][coord])
                        imgmags.append(untouched_mag[dat][coord])
                #Increase the iteration tracker and create an empty data matrix for the observation
                numobs += 1
                data = np.matrix(np.zeros((100,100)))
                #Loop through the list to get the RA and DEC positional values and the flux value
                for val in range (0, len(imgdecs)):
                        smooth_x = imgras[val]
                        smooth_y = imgdecs[val]
                        smooth_m = imgmags[val]
                        #Enumerate the x and y axes, loop through these to calculate the values for the 
                        #gaussian spreading across the data matrix
                        for xind, xval in enumerate(x_axis):
                            for yind, yval in enumerate(y_axis):
                                data[xind, yind] += gaussian(1.5*(x_axis[1]-x_axis[0]), smooth_x, smooth_y, xval, yval, smooth_m)
            #Except catches cases where the observation in point-like and the data matrix should have no vlaues
            except:
                data = np.matrix(np.zeros((100,100)))
            #Add the data matrix to the list of diffuse observation components
            obs_data.append(data)
            #Initialize a new figure
            fig, ax = plt.subplots(figsize = (10,10))
            #Plot the data matrix for this observation with the necessary components thrown back
            implot = plt.imshow(np.flip(data.T, axis = 0), cmap = 'pink', interpolation = None, origin = 'upper', 
                               extent = [np.min(x_axis), np.max(x_axis), np.min(y_axis), np.max(y_axis)])
            #Add x and y labels
            plt.xlabel('RA (Deg)')
            plt.ylabel('DEC (Deg)')
            #Display colorbar scale
            plt.colorbar()
            #Invert the RA axis
            plt.gca().invert_xaxis()
            #Add a title to the plot
            plt.title(f'Observation {dat}, Second Iteration With Points Returned')
            #Save the figure with a unique identifier
            plt.savefig(f'source{index}observation{dat}radius{radius}bright{percent_total}diffusereturned.png')
            plt.show()
            plt.figure()
        #MAKE FINAL PLOT
        #Grab the original data from the DataFrame
        mags = table.loc[index][6::9]#upper level magnitude values
        ra_center = table.loc[index][8::9]#upper level ra values
        dec_center = table.loc[index][9::9]#upper level dec values
        ext_ras = table.loc[index][10::9]
        ext_decs = table.loc[index][11::9]
        #Find dimension of one square pixel
        pixel_dim = x_axis[1]-x_axis[0]
        #SHIFTING THE OBSERVAITONS
        #Create empty lists for shifting the observations
        shifts = []
        newdata = []
        #Create empty lists to store the pulled clustered points
        truepulled_ra = []
        truepulled_dec = []
        truepulled_mag = []
        #Create a None variable to stack the observations
        imstack = None    
        #Loop through the data matricies
        for img in range(len(obs_data)):
            #Expand the dimensions of the matrix
            shifted_obs = np.expand_dims(obs_data[img], axis = 0)
            #Check if this is the first observation
            if imstack is None:
                #If it is, use it to initialize the data stack
                imstack = shifted_obs
            else:
                #If it is not, add it to the image stack
                imstack = np.vstack((imstack, shifted_obs))
        #Create a smoothed image by taking the median value for all positions in the
        #imstack variable
        raw_smoothed_image = np.median(imstack, axis = 0)
        #Create a new figure
        fig, ax = plt.subplots(figsize = (10,10))
        #Divide the pixel values by the number of observations
        #CHECK that this produces apropriate flux values
        displaydata = np.flip((raw_smoothed_image/numobs).T, axis = 0)
        #Plot the final data matrix generated by the combined observations
        implot = plt.imshow(displaydata, cmap = 'pink', interpolation = None, origin = 'upper', 
                            extent = [np.min(x_axis), np.max(x_axis), np.min(y_axis), np.max(y_axis)]) 
        plt.colorbar()
        #Loop through the clustered points identified for this observation and grab the weighter RA and 
        #DEC positions. Divide the combined flux observed by the total number of observations
        #CHECK should this be total observations, or total for observations that contained this point
        for k in true_keys:
            truepulled_ra.append(clustdict[0]['obs0'][k]['weight_ra'])
            truepulled_dec.append(clustdict[0]['obs0'][k]['weight_dec'])
            numberobs = obscontriblen[k]
            truepulled_mag.append(np.sum(clustdict[0]['obs0'][k]['component_mags'])/numberobs)
        #Create a scatter plot on top of the plotted data matrix to display the final point-like structure
        for point in range(0, len(truepulled_ra)):
            plt.scatter(truepulled_ra[point], truepulled_dec[point], 
                        label = f'Point Like Structure, {truepulled_mag[point]:.4f} Jy')
        #Display the legend
        plt.legend()
        #Invert the RA axis
        plt.gca().invert_xaxis()
        #Create x and y axis labels
        plt.xlabel('RA (Deg)')
        plt.ylabel('DEC (Deg)')
        #Create a title for the final plot
        plt.title(f'Final Model for Source {index} With Point Like Structure')
        #Save the figure with a unique name and display it
        plt.savefig(f'source{index}medianobsradius{radius}bright{percent_total}falsereturned.png')
        plt.show()
        return truepulled_ra, truepulled_dec, truepulled_mag, displaydata, x_axis, y_axis
def gen_SkyModel(table, source_ids, is_big, sav_path, source_freq):
    """
    table: pandas DataFrame generated by the match_to_GLEAM function
    
    source_ident: List of source identifyers or list single source identifyer to be replaced in
                  the generated catalog. This list can be generated by the source_search function.
    
    is_big: list of length matching source_ident containing a 0 for every source that is not 
            too large to fit in the modeling area and a 1 for every source that is too large
            for the modeling area. This list can be generated by the source_search function.
            
    sav_path: absolute path to the folder where the files should be saved
    
    source_freq: an array of frequency values with units like Hz of length matching that of source_ident.
                For example, if the frequency of the observations is 1.842e8 Hz, the input should look
                like the below:
                source_freq = np.repeat(1.824e8*u.Hz, len(source_ident))
    Returns
    ---
    For a more complete descripion of the plots and print statements generated as this function runs, see
    the doccumentation for the create_plots function
    
    Upon completion, the function will print: Modeling has completed. A SkyModel has been saved to 
    savpath + Combined_models.txt
    
    Combined_models.txt is a SkyModel object. This SkyModel contains the models for the source_ident sources
    and the original models from table. Sources in source_ident and in the surrounding area modeled by those
    sources have been removed from the sources in table and are not modeled
    
    The SkyModel has several columns.
    
    SOURCE_ID contains string names for each of the generated sources. Diffuse sources are broken into 
            individual sources. Each name has an original GLEAM identifier followed by an integer number.
            If there only exists one string name with a zero following, this represents a point source.
            If there are several names with identical strings but different integer numbers, these 
            components belong to the same source
            
    RA_J2000 [deg] contains decimal degree coordinates for the RA positions of the points
    
    DEC_J2000 [deg] contains decimal degree coordinates for the DEC positions of the points
    
    Flux [Jy] contains the flux values for these points 
    
    Frequency [Hz] contains the frequency values for the points
    
    """
    #Cut is the minimum flux value in Jy a pixel must satisfy to be comsidered part of the source and
    #added to the final catalog.
    cut = 1e-4
    #Generate empty lists to store the data in
    idents = []
    sourceras = []
    sourcedecs = []
    sourcestokes = []
    sourcefs = []
    extended_model_groups = []
    num = -1
    remove_inds = []
    #Loop through the sources in the input list
    for i in range(len(source_ids)):
        #Generate a list of the indices in the table
        inds = table.index.to_list()
        #Grab the RA and DEC values for all the sources in the table
        r= table['RA']
        d = table['DEC']
        #Grab the identifyer, RA, and DEC for the current source in the loop, check the value from the 
        #is_big list
        identif = source_ids[i]
        big = is_big[i]
        source_r = tab.loc[identif]['RA']
        source_d = tab.loc[identif]['DEC']
        #If the source is not too big for the modeling area
        if big == 0:
            #Grab the limits on the modeling area and identify sources within the larger catalog that 
            #fit within this area
            tab_ind = (r<(source_r +.25)) & (r>(source_r -.25))& (d<(source_d +.25)) & (d>(source_d -.25))
            inds_arround = tab[tab_ind].index.to_list()
            remove_inds.extend(inds_arround)
        else:
            #Do the same for a sightly larger area should the source be flagged as big
            tab_ind = (r<(source_r +.35)) & (r>(source_r -.35))& (d<(source_d +.35)) & (d>(source_d -.35))
            inds_arround = tab[tab_ind].index.to_list()
            remove_inds.extend(inds_arround)
        #Increase the variable that tracks iterations of this loop
        num +=1
        #Generate the RA, DEC, and flux arrays for the data with the x and y axes from the create_plots function
        truepulled_ra, truepulled_dec, truepulled_mag, displaydata, x_axis, y_axis = create_plots(table, source_ids[i])
        #Generate some empty lists to store the data for the pixel values
        pix_mag = []
        pix_ra = []
        pix_dec = []
        #Loop thorugh the enumerated x and y axes to locate the pixel values
        for xind, xval in enumerate(x_axis):
            for yind, yval in enumerate(y_axis):
                #Check if the data matrix pixel value at a given location has a flux value greater than
                #the specefied cutoff value
                if displaydata[xind, yind] > cut:
                    #If the pixel passes the check, append the RA and DEC coordinate values along with the 
                    #flux values to the appropriate lists
                    pix_mag.append(displaydata[xind, yind])
                    pix_ra.append(xval)
                    pix_dec.append(yval)
        #Add the RA, DEC, and flux values for the clustered points that were pulled from the data matrix to 
        #the larger pixel lists
        pix_ra.extend(truepulled_ra)
        pix_dec.extend(truepulled_dec)
        pix_mag.extend(truepulled_mag)
        #Create a list of SkyCoord objects from these RA and DEC vlaues, specifying that they are in 
        #units of decimal degrees
        coords = SkyCoord(pix_ra*u.degree, pix_dec*u.degree)
        ras = coords.ra
        decs = coords.dec
        #Identify the lenght of these coordinate lists
        compind = len(ras)
        #Grab the ientifyer for the source being iterated over
        sourceid = source_ids[i]
        #Create an empty list to store component names
        compnames = []
        #Generate a unique name for every point in the data arrays from the original source 
        #identifyer and the integer index of the point in the list
        for ind in range(compind):
            indname = sourceid + '_' + str(ind)
            compnames.append(indname)
        #Create an array of the same length as the pixel lists containing the value for the observation
        #frequency specified at function call
        source_freqs = np.repeat(source_freq[i],len(pix_ra))
        #Create an array for the stokes parameter from the flux list
        #IMPROVE COMMENTS FOR THIS LINE
        stokes = np.zeros((4, 1, len(pix_ra)), dtype = float)
        stokes[0, :, :] = pix_mag
        #Create an extended model group by repeating the iteration tracker for the lengths of the pixel values
        extended_model_group = np.repeat(num, len(ras))
        #Concatenate the generated component names and component names from previous iterations
        #into one array
        idents = np.concatenate((idents, compnames), axis = 0)
        #Create Longitude and Latitude arrays from the RA and DEC coordinates of the source and concatenate
        #this with the same array from previous loops
        sourceras = np.concatenate((Longitude(sourceras, unit = u.deg), ras), axis = 0)
        sourcedecs = np.concatenate((Latitude(sourcedecs, unit = u.deg), decs), axis = 0)
        #Concatenate the frequency array with the previous iterations
        sourcefs = np.concatenate((sourcefs, source_freqs), axis = 0)
        #Concatenate the extended model group array with the array from the previous loops
        extended_model_groups = np.concatenate((extended_model_groups, extended_model_group), axis = 0)
        #Check if this is the first iteration
        if len(sourcestokes)<1:
            #If sourcestokes has not yet experienced a concatenation, it must be initialized by simply
            #extending soutcestoked with the stokes array
            sourcestokes.extend(stokes) 
        else:
            #If this is not the first iteration, concatenate stokes with the previous iteration
            sourcestokes= np.concatenate((sourcestokes, stokes), axis = 2)
    #Create an untouched version of the index array for the input DataFrame
    inds = table.index.to_list()
    #Loop through the indices and remove any names that appear in either the list of input sources
    #or in the sources identified as existing within the modeling area
    for name in inds:
        if name in source_ids:
            inds.remove(name)
        elif name in remove_inds:
            inds.remove(name)
    #Check the intersection of the resulting list and original input to ensure no repeats are present
    st = set(inds) & set(source_ids)
    if len(st) !=0:
        for nm in st:
            #If the set is not empty, remove these sources from the original list
            inds.remove(nm)
    #Create a section of the original DataFrame tha only includes the remaining portion of the original 
    #catalog that was not extracted 
    remaining_GLEAM = table[table.index.isin(inds)]
    #Create some empty lists to store the data and define an integer iteration tracker to increase
    #as the function loops
    idents_gleam = []
    sourceras_gleam = []
    sourcedecs_gleam = []
    sourcestokes_gleam = []
    sourcefs_gleam = []
    extended_model_groups_gleam = []
    num_gleam = -1
    #Loop through these remaining sources
    for source_name in remaining_GLEAM.index:
        #Check if the original input catalog identified this source as diffuse, or containing multiple
        #components. If it did, the source needs to be treated like the other modeled diffuse sources
        if np.sum(tab.loc[source_name]['RA EO GLEAM']) >0:
            #Increase the integer iteration tracker and pull the information for the diffuse
            #components from the original DataFrame
            num_gleam +=1
            g_ra = tab.loc[source_name]['RA EO GLEAM']
            g_dec = tab.loc[source_name]['DEC EO GLEAM']
            g_mag = tab.loc[source_name]['Mag EO GLEAM']
            #Create SkyCoord objects for these points
            coords_g = SkyCoord(g_ra*u.degree, g_dec*u.degree)
            ras_g = coords_g.ra
            decs_g = coords_g.dec
            #Check the length of the coordinates
            compind = len(g_ra)
            #Create empty list to store component names
            compnames = []
            #Generate a unique name for every point in the data arrays from the original source 
            #identifyer and the integer index of the point in the list
            for ind in range(compind):
                indname = source_name + '_' + str(ind)
                compnames.append(indname)
            #Create an array of the same length as the pixel lists containing the value for the observation
            #frequency specified at function call
            source_freqs_gleam = np.repeat(source_freq[0],len(g_ra))
            #Create an array for the stokes parameter from the flux list
            stokes_gleam = np.zeros((4, 1, len(g_ra)), dtype = float)
            stokes_gleam[0, :, :] = g_mag
            #Create an extended model group by repeating the iteration tracker for the lengths of the pixel values
            extended_model_group_gleam = np.repeat(num_gleam, len(g_ra))
            #Concatenate the generated component names and component names from previous iterations
            #into one array
            idents_gleam = np.concatenate((idents_gleam, compnames), axis = 0)
            #Create Longitude and Latitude arrays from the RA and DEC coordinates of the source and concatenate
            #this with the same array from previous loops
            sourceras_gleam = np.concatenate((Longitude(sourceras_gleam, unit = u.deg), ras_g), axis = 0)
            sourcedecs_gleam = np.concatenate((Latitude(sourcedecs_gleam, unit = u.deg), decs_g), axis = 0)
            #Concatenate the frequency array with the previous iterations
            sourcefs_gleam = np.concatenate((sourcefs_gleam, source_freqs_gleam), axis = 0)
            #Concatenate the extended model group array with the array from the previous loops
            extended_model_groups_gleam = np.concatenate((extended_model_groups_gleam, extended_model_group_gleam), axis = 0)
            #Check if this is the first iteration
            if len(sourcestokes_gleam)<1:
                #If sourcestokes has not yet experienced a concatenation, it must be initialized by simply
                #extending soutcestoked with the stokes array
                sourcestokes_gleam.extend(stokes_gleam) 
            else:
                #If this is not the first iteration, concatenate stokes with the previous iteration
                sourcestokes_gleam= np.concatenate((sourcestokes_gleam, stokes_gleam), axis = 2)
        #For casees where the source in the original catalog does not have diffuse components
        else:
            #Increase the integer that tracks the loop iterations
            num_gleam +=1
            #Get the RA and DEC positions and flux values from the original input DataFrame 
            g_ra = tab.loc[source_name]['RA']
            g_dec = tab.loc[source_name]['DEC']
            g_mag = tab.loc[source_name]['Mag GLEAM']
            #Create SkyCoord objects for these positions
            coords_g = SkyCoord(g_ra*u.degree, g_dec*u.degree)
            ras_g = coords_g.ra
            decs_g = coords_g.dec
            #Create a singular component name for the object
            compnames = [source_name + '_' + str(0)]
            #Create an array of the same length as the pixel lists containing the value for the observation
            #frequency specified at function call
            source_freqs_gleam = np.repeat(source_freq[0],1)
            #Create an array for the stokes parameter from the flux list
            stokes_gleam = np.zeros((4, 1, 1), dtype = float)
            stokes_gleam[0, :, :] = g_mag
            #Create an extended model group by repeating the iteration tracker for the lengths of the pixel values
            extended_model_group_gleam = np.repeat(num_gleam, 1)
            #Concatenate the generated component names and component names from previous iterations
            #into one array
            idents_gleam = np.concatenate((idents_gleam, compnames), axis = 0)
            #Create Longitude and Latitude arrays from the RA and DEC coordinates of the source and concatenate
            #this with the same array from previous loops
            sourceras_gleam = np.concatenate((Longitude(sourceras_gleam, unit = u.deg), [ras_g]), axis = 0)
            sourcedecs_gleam = np.concatenate((Latitude(sourcedecs_gleam, unit = u.deg), [decs_g]), axis = 0)
            #Concatenate the frequency array with the previous iterations
            sourcefs_gleam = np.concatenate((sourcefs_gleam, source_freqs_gleam), axis = 0)
            #Concatenate the extended model group array with the array from the previous loops
            extended_model_groups_gleam = np.concatenate((extended_model_groups_gleam, extended_model_group_gleam), axis = 0)
            #Check if this is the first iteration
            if len(sourcestokes_gleam)<1:
                #If sourcestokes has not yet experienced a concatenation, it must be initialized by simply
                #extending soutcestoked with the stokes array
                sourcestokes_gleam.extend(stokes_gleam) 
            else:
                #If this is not the first iteration, concatenate stokes with the previous iteration
                sourcestokes_gleam= np.concatenate((sourcestokes_gleam, stokes_gleam), axis = 2)
    #GENERATING THE SKYMODELS
    #SkyModel for input objects
    skyobj = SkyModel(
            name=idents,  # string names in a list or array
            ra=sourceras,      # astropy Longitude array
            dec=sourcedecs,    # astropy Latitude array
            stokes=sourcestokes,  # astopy Quantity units like Jy, shape (4, 1, Nsrcs)
            spectral_type="flat", #string spectral type
            reference_frequency=sourcefs, # astropy Quantity, units like Hz
            extended_model_group=np.array(list(map(str, extended_model_groups))))  # integer array where components of an extended source all have the same value, sources with only one component should have a -1
    #skyobj.write_text_catalog(savpath + 'skyobj_EoR0_test_new.txt')
    #SkyModel for remaining catalog objects
    skyobj_g = SkyModel(
            name=idents_gleam,  # string names in a list or array
            ra=sourceras_gleam,      # astropy Longitude array
            dec=sourcedecs_gleam,    # astropy Latitude array
            stokes=sourcestokes_gleam,  # astopy Quantity units like Jy, shape (4, 1, Nsrcs)
            spectral_type="flat", #string spectral type
            reference_frequency=sourcefs_gleam, # astropy Quantity, units like Hz
            extended_model_group=np.array(list(map(str, extended_model_groups_gleam))))  # integer array where components of an extended source all have the same value, sources with only one component should have a -1
    #skyobj_g.write_text_catalog(savpath + 'new_gleam.txt')
    skyobj.concat(skyobj_g)
    skyobj.write_text_catalog(savpath + 'Combined_models.txt')
    #print(f'Modeling has completed. A SkyModel has been saved to {savpath + 'Combined_models.txt'}')
    return print(f'Modeling has completed. A SkyModel has been saved to {savpath}Combined_models.txt')

