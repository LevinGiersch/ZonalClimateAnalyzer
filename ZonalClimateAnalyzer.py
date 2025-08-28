#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Import Packages

from pathlib import Path
import os
import json
import shutil
import zipfile
import gzip
import filetype

import numpy as np

from tqdm import tqdm
from pprint import pprint
from itertools import chain

from urllib.parse import urljoin, urlparse
import requests
from bs4 import BeautifulSoup

import folium
import geopandas as gpd
from pyproj import CRS
import rasterio
from rasterio.crs import CRS
from rasterstats import zonal_stats

import matplotlib
import matplotlib.pyplot as plt


# # Get Shapefile

# In[2]:


def check_crs(shapefile:str):
    '''
    Takes a path to a Shapefile (as String) as Input.
    Reads shapefile.
    Returns True if it has a valid CRS,
    Returns False if it die NOT have a valid CRS.
    '''

    # Load the shapefile
    gdf = gpd.read_file(shapefile)
    
    # Check if CRS is defined
    valid_crs = gdf.crs
    if valid_crs:
        print("\nCRS is defined:", gdf.crs)
        return True
    else:
        print("\nCRS is NOT defined.")
        return False


# In[3]:


def get_shp():
    '''
    Lets the user input a path to the shapefile in the terminal.
    Checks if the path leads to a shapefile.
    Checks if the shapefile has a valid CRS.
    Gives feedback depending on the user input.
    Returns path to the shapefile if Valid shp and CRS are found.
    '''

    print('\n'+'#'*64)
    print('This Program lets you analyze the climate history of any area within Germany.')
    print('You only need a shapefile defining the area you want to analyze.')

    while True:  # This function runs until input is valid
        shp_input = input('\nEnter the path to the shapefile here: ').strip()
        shp_path = Path(shp_input)

        # If the input is a file path
        if shp_path.is_file() and shp_path.suffix.lower() == '.shp':
            try:
                gdf = gpd.read_file(shp_path)
                if gdf.crs:
                    print('\nValid shapefile with valid CRS found.\n')
                    return shp_path
                else:
                    print('\nShapefile found, but CRS is not defined.\n')
                    continue
            except Exception as e:
                print(f'Error reading shapefile: {e}\n')
                continue

        # If the input is a folder path, search for any .shp file inside
        elif shp_path.is_dir():
            shp_files = list(shp_path.glob("*.shp"))
            if shp_files:
                print(f'Found shapefiles: {[f.name for f in shp_files]}. \nPlease append the filename to the path and try again.\n')
                continue
            else:
                print('\nNo shapefile found in the folder.')
                continue

        else:
            print('\nInvalid path or not a shapefile.')
            continue


# # Download Data

# In[4]:


def list_of_dwd_data(file_types=['.asc.gz', '.pdf', '.zip']):
    '''
    Creates a list containing all dwd files to download.

    Parameters:
        file_types: string, ending that correspondes to the filetype that we want to download

    Returns:
        links: list containing links to all files to download
    '''

    # Get list of all download locations containing the data to download:
    base_download_location = 'https://opendata.dwd.de/climate_environment/CDC/grids_germany/annual/'
    folder_download_locations = [
        'air_temperature_max/',
        'air_temperature_mean/',
        'air_temperature_min/',
        'drought_index/',
        #'erosivity/',
        'frost_days/',
        'hot_days/', 
        'ice_days/',
        'phenology/',
        'precipGE10mm_days/',
        'precipGE20mm_days/',
        'precipGE30mm_days/',
        'precipitation/',
        #'radiation_diffuse/',
        #'radiation_direct/',
        #'radiation_global/',
        'snowcover_days/',
        'summer_days/',
        'sunshine_duration/',
        'vegetation_begin/',
        'vegetation_end/'
    ]
    download_locations = [base_download_location+f for f in folder_download_locations]

    # Create a list containing all links
    links = []
    for location in download_locations:
        response = requests.get(location)
        if response.status_code != 200:
            raise Exception(f'\nFailed to retrieve the webpage: {location}')

        soup = BeautifulSoup(response.text, 'html.parser')

        # Build absolute URLs
        found = [
            location + a['href']
            for a in soup.find_all('a', href=True)
            if a['href'].endswith(tuple(file_types))
        ]
        links.append(found)

    return list(chain.from_iterable(links))
        
    return links


# In[5]:


def download_dwd_data(links:list[str], dest_dir:str, timeout:int=30):
    """
    Download files from a list of full URLs into a target directory.

    Parameters:
        links (list[str]): List of full download URLs.
        dest_dir (str or Path): Local directory where the files will be saved.
        timeout (int, optional): Maximum number of seconds to wait for a server response. Defaults to 30.

    Returns:
        None
    """

    # TODO: Implement a test if data is already in folder and then skip
    
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)

    for file_url in tqdm(links,
                     desc='',
                     bar_format='{l_bar}{bar:40}| ({n_fmt}/{total_fmt}) Downloading files. This may take a few minutes.',
                     ncols=120):
        filename = Path(urlparse(file_url).path).name
        file_path = dest / filename

        with requests.get(file_url, stream=True, timeout=timeout) as r:
            if r.status_code != 200:
                print(f'\nFehler {r.status_code}: {file_url}')
                continue
            with file_path.open('wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)


# # Process the Data

# In[6]:


def list_of_files(folder:str, file_type='.gz'):
    '''
    Returns list containing all filesnames in folder with the ending file_type.

    Args:
        folder: string, path in filesystem including target folder
        file_type: string, ending that correspondes to the filetype that we want to download
    Return:
        List of all files of file_type within folder
    '''

    files = sorted([os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(file_type)])
    return files


# In[7]:


def rename_dwd_file(file:str):
        """
        Takes filename as string as input.
        Removes everything except for the core name and the year from the dwd filename.
        Returns new filename as string.
        """
        file = file.replace('grids_germany_annual_', '')
        file = file.removesuffix('.asc.gz')
        if file.endswith('_1917') or file.endswith('_2017'):
            pass
        else:
            file = file.removesuffix('17')
        file = file.removesuffix('_')
        file = file+'.asc'

        return file


# In[8]:


def decompress_file(file:str):
    """
    Decompress files and saves a copy in the same folder.

    Args:
        file (str): path to file (input file)
    Return:
        decompressed_file (str): Path to decompressed file (output file)
    """

    # Check Filetype:
    ft = filetype.guess(file)
    ft_ext = ft.extension

    if ft_ext == 'asc' or ft_ext == 'tif':
        print(f'Filetype is already {ft_ext}, no decompression needed')

    elif ft_ext == 'gz':
        decompressed_file = Path(rename_dwd_file(file))         # Name for output file
        with gzip.open(file, mode='rb') as f_in:                # Open .gz file and decompress it
            with open(decompressed_file, mode='wb') as f_out:   # Create decompressed file
                shutil.copyfileobj(f_in, f_out)                 # Copy content of compressed file to decompressed file

    # Unpack the mis‑labelled “…asc.gz” archive (which is really a ZIP)
    elif ft_ext == 'zip':
        decompressed_file = Path(rename_dwd_file(file))         # Name for output file
        zip_path = Path(file).expanduser().resolve()
    
        with zipfile.ZipFile(zip_path) as zf:
            # find the first .asc inside the archive
            asc_members = [m for m in zf.namelist() if m.lower().endswith(".asc")]
            if not asc_members:
                raise ValueError("No .asc file found in archive.")
            # write it directly to the desired location
            with zf.open(asc_members[0]) as src, decompressed_file.open("wb") as dst:
                dst.write(src.read())
                
    else:
        print(f'Error while decompressing File: {file}. \nFiletype {ft_ext} is not supported. \nSupported filetypes are: .gz, .zip, .asc or .tif\n')
        raise TypeError(f'Error while decompressing File. Filetype is not supported. Supported filetypes are: .gz, .zip, .asc or .tif')

    #print(f'decompress_file: Successfully decompressed \n"{file}" to \n"{decompressed_file}".\n')

    return decompressed_file


# In[9]:


def asc_to_tif_add_crs(asc_input:str, prj_txt:str):
    """
    Takes asc_input, adds crs (from prj.txt), saves it as .tif in the same folder.

    Args:
        asc_input (str): path to asc file (input file)
        prj_txt (str): path to prj file (text file that contains projection information)
    
    Returns:
        tif_output (tif): path/to/output.tif file
    """

    # CRS from .prj_file
    with open(prj_file, 'r') as f:
        prj_txt = f.read()
    crs = CRS.from_wkt(prj_txt)

    # Read the .asc file
    with rasterio.open(asc_input) as src:
        data = src.read(1)
        profile = src.profile

    # Update Profile with CRS
    profile.update({
        'driver': 'GTiff',
        'crs': crs,
        'dtype': rasterio.float32,
        'compress': 'lzw'
    })

    tif_output = asc_input.replace('.asc', '')+'.tif'

    # Write to a new GeoTIFF file with CRS assigned
    with rasterio.open(tif_output, 'w', **profile) as dst:
        dst.write(data.astype(rasterio.float32), 1)
    
    #print(f'asc_to_tif_add_crs: Successfully transformed \n"{asc_input}" to \n"{tif_output}" \nand added {crs}.\n')

    return tif_output


# In[10]:


def delete_raster_files(folder_path:str):
    """
    Deletes all .asc, .asc.gz and .zip files within the folder..

    Args:
        folder_path (str): Path to folder containing the files.
    """
    
    folder = Path(folder_path)
    patterns = ['*.asc', '*.asc.gz', '*.zip']
    deleted_files = 0

    for pattern in patterns:
        for file in folder.glob(pattern):
            file.unlink()
            deleted_files += 1


# In[11]:


def change_shp_crs(shp_input:str, prj_txt:str):
    """
    Takes shp_input, transforms to crs (from prj.txt), outputs as shp_output

    Args:
        shp_input (str): path to shp file (input file)
        prj_txt (str): path to prj file (text file that contains projection information)
    
    Returns:
        shp_output (str): path/to/output.tif file
    """

    # CRS from .prj_file
    with open(prj_file, 'r') as f:
        prj_txt = f.read()

    target_crs = CRS.from_wkt(prj_txt)

    # Read Shapefile
    gdf = gpd.read_file(shp_input)

    # Check if CRS is defined
    if gdf.crs is None:
        raise ValueError('Input shapefile crs is undefined. Set correct crs. (eg. gdf.set_crs())')

    # Transform shp_input to target_crs
    gdf_transformed = gdf.to_crs(target_crs)

    # Create Output Folder
    shp_folder_path = Path.cwd() / 'shp'
    shp_folder_path.mkdir(parents=True, exist_ok=True)

    # Define output name
    shp_output = str(shp_folder_path / Path(str(shp_input).replace('.shp', '')+'_'+str(target_crs).replace(':','')+'.shp').name)

    # Save transformed shp to shp_output
    gdf_transformed.to_file(shp_output, encoding='utf-8')

    #print(f'change_shp_crs: Successfully copied \n"{shp_input}" to \n"{shp_output}" \nand added {target_crs}.\n')

    return shp_output


# In[12]:


def dissolve_shp(shp_input:str):
    """
    Takes shp_input, dissolve all polygon features into one, outputs as dissolved_shp

    Args:
        shp_input (str): path to shp file (input file)
    
    Returns:
        shp_output (str): path/to/output.tif file
    """

    # Read Shapefile
    gdf = gpd.read_file(shp_input)

    # Dissolve features in gdf
    gdf_dissolved = gdf.dissolve()

    # Create Output Folder
    shp_folder_path = Path.cwd() / 'shp'
    shp_folder_path.mkdir(parents=True, exist_ok=True)

    # Define output name
    shp_output = str(shp_folder_path / Path(str(shp_input).replace('.shp', '')+'_dissolved'+'.shp').name)

    # Save transformed shp to shp_output
    gdf_dissolved.to_file(shp_output, encoding='utf-8')
    
    return shp_output


# In[13]:


def calculate_zonal_stats(shp:str, tif:str):
    """
    Calculates zonal stats of the tif for each feature in the shp.

    Args:
        shp (str): path to shp file
        tif (str): path to tif file
    
    Returns:
        stats (list[dict]): list of dictionarys which contain min, max, mean and count of the raster data for each poly.
    """

    # Set all_touched to False if you want to include only raster-cells that are completely within the shapefile.
    stats = zonal_stats(shp, tif, all_touched=True)

    return stats


# In[14]:


def zonal_climate_analysis(shp_input:str, raster_folder:str, prj_file:str):
    """
    Perform rasterstats calculation on shp_input with every raster file in the raster_folder.

    Args:
        shp_input (str): path to shp file (input file) to perform calculations on
        raster_folder (str): path to folder containing all raster files to perform the rasterstats calculations with. has to be in .asc.gz file format
        prj_txt (str): path to prj file (text file that contains projection information)
    
    Creates:
        rasterstats_dict (dict{str:[{}]}): dict containing the name of the raster file as key and the corresponding rasterstats as a list of dicts as values.
    Returns:
        json_output_path_name (str): path to the created json file conatining rasterstats calculations.
        shp_crs_dissolved (str): path to the dissolved shapefile with transformed crs the rasterstats where calculated on.
    """

    # Prepare shapefile:
    shp_crs = change_shp_crs(shp_input, prj_file)
    shp_crs_dissolved = dissolve_shp(shp_crs)
    
    # Create list of compessed .asc.gz rasterfiles:
    files_asc_gz = list_of_files(raster_folder, file_type='.asc.gz')

    # Decompress rasterfiles:
    for f in tqdm(files_asc_gz,
                  desc='',
                  bar_format='{l_bar}{bar:40}| ({n_fmt}/{total_fmt}) Decompressing files.',
                  ncols=120):
        decompress_file(f)

    # Create list of decompressed .asc rasterfiles:
    files_asc = list_of_files(raster_folder, file_type='.asc')

    # Transform decompressed files to tif and add crs:
    for f in tqdm(files_asc,
                  desc='',
                  bar_format='{l_bar}{bar:40}| ({n_fmt}/{total_fmt}) Transforming files to the right format.',
                  ncols=120):
        asc_to_tif_add_crs(f, prj_file)

    # Create list .tif rasterfiles
    files_tif = list_of_files(raster_folder, file_type='.tif')

    # Create list containing rasterstats:
    rasterstats_list = []

    # Iterate over files_tif and perform rasterstats calculations on each rasterfile and the shapefile:
    for f in tqdm(files_tif,
                  desc='',
                  bar_format='{l_bar}{bar:40}| ({n_fmt}/{total_fmt}) Calculating rasterstats.',
                  ncols=120):
        rasterstats_list.append(calculate_zonal_stats(shp_crs_dissolved, f)) # Append rasterstats to rasterstats_list

    # Combine rasterstats and the name of the raster the stats are calculated with
    raster_path = str(Path.cwd() / 'climate_environment_CDC_grids_germany_annual')
    filenames = filenames = [fn.replace(raster_path + '/', '').replace('.tif', '') for fn in files_tif] # Ceate list with filenames without path and type  

    # Delete deprecated files
    delete_raster_files(raster_path)

    # Create dict
    rasterstats_dict = dict(zip(filenames, rasterstats_list))

    # Convert rasterstats_dict to better json format:
    rasterstats_json = {}

    for key, value in rasterstats_dict.items():
        name = key[:-5]
        year = key[-4:]
        if name not in rasterstats_json:
            rasterstats_json[name] = {}
        rasterstats_json[name][year] = value

    #pprint(rasterstats_json)

    # Export dict as json:
    path_to_shp = Path(shp_crs_dissolved)
    shp_name = path_to_shp.name
    json_output_path_name = shp_name.replace('.shp','')+'_rasterstats.json'
    
    with open(json_output_path_name, 'w') as rs_json:
        json.dump(rasterstats_json, rs_json)

    return json_output_path_name, shp_crs_dissolved


# # Visualize

# In[15]:


def years_values(parameter_name:str):
    """
    Arguments:
        parameter_name (str): key in rasterstats.json dictionary
    Returns:
        title (str): parameter_name
        years (list): list of years
        values_max (list): list of max values
        values_mean (list): list of mean values
        values_min (list): list of min values
    """

    # Title:
    title = parameter_name

    # Years:
    years = []
    for year in rs[title]:
        years.append(year)

    # Values max:
    values_max = []
    for year in rs[parameter_name].values():
        for entry in year:
            values_max.append(entry['max'])

    # Values mean:
    values_mean = []
    for year in rs[parameter_name].values():
        for entry in year:
            values_mean.append(entry['mean'])

    # Values mean:
    values_min = []
    for year in rs[parameter_name].values():
        for entry in year:
            values_min.append(entry['min'])
    
    return title, years, values_max, values_mean, values_min


# In[16]:


def create_map(shapefile:str):
    '''
    Takes path to shapefile as string as input.
    Creates interactive map as html.
    Adds area and perimeter as tooltips on hover in the html map.
    '''
    
    shp_path = Path(shapefile)
    gdf = gpd.read_file(shapefile)
    if gdf.crs is None:
        raise ValueError("CRS is missing. Set a CRS before running.")

    # Compute in a local metric CRS
    gdf_m = gdf.to_crs(gdf.estimate_utm_crs())
    is_poly = gdf_m.geom_type.str.contains("polygon", case=False, na=False)
    gdf["area_km2"] = (gdf_m.area.where(is_poly)) / 1_000_000
    gdf["perim_km"] = (gdf_m.length.where(is_poly)) / 1_000
    gdf["shapefile"] = shp_path.stem

    # interactive map
    m = gdf.to_crs(4326).explore(
        tooltip=["shapefile", "area_km2", "perim_km"],
        popup=False
    )
    
    # Save Map
    mapname = 'map.html'
    map_folder_path = Path.cwd() / 'output'
    map_folder_path.mkdir(parents=True, exist_ok=True)
    map_path = str(map_folder_path / mapname)
    m.save(map_path)
    print(f'Successfully created and saved map: {mapname}')


# In[17]:


def plot_air_temp_min_mean_max():
    # Air Temp min mean max
    plt.close()

    # Create a figure containing a single Axes.
    fig, ax = plt.subplots()

    # List oft startyears
    startyears = []
    
    # Max Temp
    title, years, values_max, t_max, values_min = years_values('air_temp_max')
    t_max = [t_max[i]/10 for i in range(len(t_max))]                            # 1/10 so it is in degrees noch in degrees/10
    startyears.append(years[0])

    # Mean Temp
    title, years, values_max, t_mean, values_min = years_values('air_temp_mean')
    t_mean = [t_mean[i]/10 for i in range(len(t_mean))]
    startyears.append(years[0])

    # Min Temp
    title, years, values_max, t_min, values_min = years_values('air_temp_min')
    t_min = [t_min[i]/10 for i in range(len(t_min))]
    startyears.append(years[0])

    # Crop to the same start-year
    common_startyear = int(max(startyears))
    t_max  = [v for y, v in zip(map(int, years), t_max)  if y >= common_startyear]
    t_mean = [v for y, v in zip(map(int, years), t_mean) if y >= common_startyear]
    t_min  = [v for y, v in zip(map(int, years), t_min)  if y >= common_startyear]

    # Plot
    ax.plot(years, t_max, color='red', label='Maximale Lufttemperatur')
    ax.plot(years, t_mean, color='Black', label='Mittlere Lufttemperatur')
    ax.plot(years, t_min, color='blue', label='Minimale Lufttemperatur')

    # Fill between lines
    ax.fill_between(years, t_max, t_mean, color='red', alpha=0.25)
    ax.fill_between(years, t_mean, t_min, color='blue', alpha=0.25)

    # Gridlines:
    ax.grid(color='lightgrey', linewidth=0.5)

    # Plot Customization
    fig.set_size_inches(6.3*2, 3.15*2)
    ax.set_ylim([0, max(t_max)*1.2])
    ax.set_xlim(['1954', max(years)])
    #ax.set_title('Eis- und Frosttage')
    ax.set_xlabel('Jahre')
    ax.set_ylabel('Temperatur in °C')
    ax.legend(loc='upper right')
    plt.xticks(rotation=45)
    for label in ax.xaxis.get_ticklabels():  # Iterate over all ticklabels
        if int(label.get_text()) % 5 == 0:   # Check if ticklabel is dividable by 5
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Save Plot
    plotname = 'min_mean_max_temp'+'_plot.png'
    plot_folder_path = Path.cwd() / 'output'
    plot_folder_path.mkdir(parents=True, exist_ok=True)
    plot_path = str(plot_folder_path / plotname)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    print(f'Successfully created and saved plot: {plotname}')


# In[18]:


def plot_frost_ice_days():
    # Frost and Ice Days
    plt.close()

    # Create a figure containing a single Axes.
    fig, ax = plt.subplots()

    # Ice Days
    title, years, values_max, values_mean_id, values_min = years_values('ice_days')

    # Frost Days
    title, years, values_max, values_mean_fd, values_min = years_values('frost_days')

    # List containing 365 (days per year) as many times as there are years (upper limit)
    days_in_year_max = [365]*len(years)

    # List containing 365 (days per year) as many times as there are years (lower limit)
    days_in_year_min = [0]*len(years)

    # Plot
    ax.plot(years, values_mean_fd, color='lightblue', label='Frosttage (min 0°C)')
    ax.plot(years, values_mean_id, color='darkblue', label='Eistage (max 0°C)')

    # Fill between lines
    ax.fill_between(years, values_mean_fd, values_mean_id, color='lightblue', alpha=0.25)
    ax.fill_between(years, values_mean_id, days_in_year_min, color='darkblue', alpha=0.25)

    # Gridlines:
    ax.grid(color='lightgrey', linewidth=0.5)

    # Plot Customization
    fig.set_size_inches(6.3*2, 3.15*2)
    ax.set_ylim([0, max(values_mean_fd)*1.2])
    ax.set_xlim(['1954', max(years)])
    #ax.set_title('Eis- und Frosttage')
    ax.set_xlabel('Jahre')
    ax.set_ylabel('Tage')
    ax.legend(loc='upper right')
    plt.xticks(rotation=45)
    for label in ax.xaxis.get_ticklabels():  # Iterate over all ticklabels
        if int(label.get_text()) % 5 == 0:   # Check if ticklabel is dividable by 5
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Save Plot
    plotname = 'ice_frost_days'+'_plot.png'
    plot_folder_path = Path.cwd() / 'output'
    plot_folder_path.mkdir(parents=True, exist_ok=True)
    plot_path = str(plot_folder_path / plotname)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    print(f'Successfully created and saved plot: {plotname}')


# In[19]:


def plot_snowcover_days():
    # Snowcover Days
    plt.close()

    # Create a figure containing a single Axes.
    fig, ax = plt.subplots()

    # Snowcover Days
    title, years, values_max, values_mean_snd, values_min = years_values('snowcover_days')

    # List containing 365 (days per year) as many times as there are years (upper limit)
    days_in_year_max = [365]*len(years)

    # List containing 365 (days per year) as many times as there are years (lower limit)
    days_in_year_min = [0]*len(years)

    # Plot
    ax.plot(years, values_mean_snd, color='lightblue', label='Tage mit > 1cm Schneehöhe')

    # Fill between lines
    ax.fill_between(years, values_mean_snd, days_in_year_min, color='lightblue', alpha=0.25)
    
    # Gridlines:
    ax.grid(color='lightgrey', linewidth=0.5)

    # Plot Customization
    fig.set_size_inches(6.3*2, 3.15*2)
    ax.set_ylim([0, max(values_mean_snd)*1.2])
    ax.set_xlim(['1954', max(years)])
    #ax.set_title('Eis- und Frosttage')
    ax.set_xlabel('Jahre')
    ax.set_ylabel('Tage')
    ax.legend(loc='upper right')
    plt.xticks(rotation=45)
    for label in ax.xaxis.get_ticklabels():  # Iterate over all ticklabels
        if int(label.get_text()) % 5 == 0:   # Check if ticklabel is dividable by 5
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Save Plot
    plotname = 'snowcover_days'+'_plot.png'
    plot_folder_path = Path.cwd() / 'output'
    plot_folder_path.mkdir(parents=True, exist_ok=True)
    plot_path = str(plot_folder_path / plotname)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    print(f'Successfully created and saved plot: {plotname}')


# In[20]:


def plot_summer_hot_days():
    # Summer and Hot Days
    plt.close()

    # Create a figure containing a single Axes.
    fig, ax = plt.subplots()

    # Ice Days
    title, years, values_max, values_mean_sd, values_min = years_values('summer_days')

    # Frost Days
    title, years, values_max, values_mean_hd, values_min = years_values('hot_days')

    # List containing 365 (days per year) as many times as there are years (upper limit)
    days_in_year_max = [365]*len(years)

    # List containing 365 (days per year) as many times as there are years (lower limit)
    days_in_year_min = [0]*len(years)

    # Plot
    ax.plot(years, values_mean_sd, color='orange', label='Sommertage (max 25°C)')
    ax.plot(years, values_mean_hd, color='red', label='Heiße Tage (max 30°C)')

    # Fill between lines
    ax.fill_between(years, values_mean_hd, values_mean_sd, color='orange', alpha=0.25)
    ax.fill_between(years, values_mean_sd, days_in_year_min, color='red', alpha=0.25)

    # Gridlines:
    ax.grid(color='lightgrey', linewidth=0.5)

    # Plot Customization
    fig.set_size_inches(6.3*2, 3.15*2)
    ax.set_ylim([0, max(values_mean_sd)*1.2])
    ax.set_xlim(['1954', max(years)])
    #ax.set_title('Heiße- und Sommertage')
    ax.set_xlabel('Jahre')
    ax.set_ylabel('Tage')
    ax.legend(loc='upper right')
    plt.xticks(rotation=45)
    for label in ax.xaxis.get_ticklabels():  # Iterate over all ticklabels
        if int(label.get_text()) % 5 == 0:   # Check if ticklabel is dividable by 5
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Save Plot
    plotname = 'summer_hot_days'+'_plot.png'
    plot_folder_path = Path.cwd() / 'output'
    plot_folder_path.mkdir(parents=True, exist_ok=True)
    plot_path = str(plot_folder_path / plotname)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    print(f'Successfully created and saved plot: {plotname}')


# In[21]:


def plot_precipitaion():
    # Precipitation
    plt.close()

    # Create a figure containing a single Axes.
    fig, ax = plt.subplots()

    # Precipitation
    title, years, values_max, values_mean_pp, values_min = years_values('precipitation')

    # Precipitation
    title, yearsdi, values_max, values_mean_di, values_min = years_values('drought_index')
    values_mean_di = [di*10 for di in values_mean_di]

    # Plot
    ax.plot(years, values_mean_pp, color='blue', label='Niederschlag in mm')
    ax.fill_between(years, values_mean_pp, color='blue', alpha=0.25)
    ax.plot(yearsdi, values_mean_di, color='orange', label='Trockenheitsindex (mm/°C)')

    # Gridlines:
    ax.grid(color='lightgrey', linewidth=0.5)

    # Plot Customization
    fig.set_size_inches(6.3*2, 3.15*2)
    ax.set_ylim([0, max(values_mean_pp)*1.2])
    ax.set_xlim(['1954', max(years)])
    #ax.set_title('Niederschlag und Trockenheitsindex')
    ax.set_xlabel('Jahre')
    ax.set_ylabel('Niederschlag (mm)')
    ax.legend(loc='upper right')
    plt.xticks(rotation=45)
    for label in ax.xaxis.get_ticklabels():  # Iterate over all ticklabels
        if int(label.get_text()) % 5 == 0:   # Check if ticklabel is dividable by 5
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Save Plot
    plotname = 'precipitation_drought'+'_plot.png'
    plot_folder_path = Path.cwd() / 'output'
    plot_folder_path.mkdir(parents=True, exist_ok=True)
    plot_path = str(plot_folder_path / plotname)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    print(f'Successfully created and saved plot: {plotname}')


# In[22]:


def plot_precipitaion_days():
    # Precipitation Days
    plt.close()

    # Create a figure containing a single Axes.
    fig, ax = plt.subplots()

    # 10mm
    title, years, values_max, p10, values_min = years_values('precipGE10mm_days')

    # 20mm
    title, years, values_max, p20, values_min = years_values('precipGE20mm_days')
    #p1020 = [p10[i]+p20[i] for i in range(len(p10))]

    # 30mm
    title, years, values_max, p30, values_min = years_values('precipGE30mm_days')
    #p102030 = [p10[i]+p20[i]+p30[i] for i in range(len(p10))]

    # List containing 365 (days per year) as many times as there are years (lower limit)
    days_in_year_min = [0]*len(years)

    # Plot
    ax.plot(years, p10, color='lightblue', label='Anzahl der Tage mit Niederschlagshöhe >= 10 mm')
    ax.plot(years, p20, color='blue', label='Anzahl der Tage mit Niederschlagshöhe >= 20 mm')
    ax.plot(years, p30, color='darkblue', label='Anzahl der Tage mit Niederschlagshöhe >= 30 mm')

    # Fill between lines
    ax.fill_between(years, p10, p20, color='lightblue', alpha=0.25)
    ax.fill_between(years, p20, p30, color='blue', alpha=0.25)
    ax.fill_between(years, p30, days_in_year_min, color='darkblue', alpha=0.25)

    # Gridlines:
    ax.grid(color='lightgrey', linewidth=0.5)

    # Plot customization
    fig.set_size_inches(6.3*2, 3.15*2)
    ax.set_ylim([0, max(p10)*1.2])
    ax.set_xlim(['1954', max(years)])
    #ax.set_title('Anzahl der Niederschlagstage')
    ax.set_xlabel('Jahre')
    ax.set_ylabel('Tage')
    ax.legend(loc='upper right')
    plt.xticks(rotation=45)
    for label in ax.xaxis.get_ticklabels():  # Iterate over all ticklabels
        if int(label.get_text()) % 5 == 0:   # Check if ticklabel is dividable by 5
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Save Plot
    plotname = 'precip_days'+'_plot.png'
    plot_folder_path = Path.cwd() / 'output'
    plot_folder_path.mkdir(parents=True, exist_ok=True)
    plot_path = str(plot_folder_path / plotname)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    print(f'Successfully created and saved plot: {plotname}')


# In[23]:


def plot_sunshine_duration():
    # Sunshine Duration
    plt.close()

    # Create a figure containing a single Axes.
    fig, ax = plt.subplots()

    # Sunshine Duration
    title, years, values_max, values_mean_sd, values_min = years_values('sunshine_duration')
    values_mean_sd = [sd/365 for sd in values_mean_sd]

    # Plot
    ax.plot(years, values_mean_sd, color='orange', label='Durchschnittliche Sonnenstunden pro Tag')
    ax.fill_between(years, values_mean_sd, color='orange', alpha=0.25)

    # Gridlines:
    ax.grid(color='lightgrey', linewidth=0.5)

    # Plot Customization
    fig.set_size_inches(6.3*2, 3.15*2)
    ax.set_ylim([0, max(values_mean_sd)*1.2])
    ax.set_xlim([min(years), max(years)])
    #ax.set_title('Niederschlag und Trockenheitsindex')
    ax.set_xlabel('Jahre')
    ax.set_ylabel('Sonnenstunden pro Tag')
    ax.legend(loc='upper right')
    plt.xticks(rotation=45)
    for label in ax.xaxis.get_ticklabels():  # Iterate over all ticklabels
        if int(label.get_text()) % 5 == 0:   # Check if ticklabel is dividable by 5
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Save Plot
    plotname = 'sunshine_duration'+'_plot.png'
    plot_folder_path = Path.cwd() / 'output'
    plot_folder_path.mkdir(parents=True, exist_ok=True)
    plot_path = str(plot_folder_path / plotname)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    print(f'Successfully created and saved plot: {plotname}')


# In[24]:


def plot_vegetation_begin_end():
    # Vegetation begin and vegetation end
    plt.close()

    # Create a figure containing a single Axes.
    fig, ax = plt.subplots()

    # Vegetation begin line
    title, years, values_max, values_mean_b, values_min = years_values('vegetation_begin')

    # Vegetation end line
    title, years, values_max, values_mean_e, values_min = years_values('vegetation_end')

    # List containing 365 (days per year) as many times as there are years (upper limit)
    days_in_year_max = [365]*len(years)

    # List containing 365 (days per year) as many times as there are years (lower limit)
    days_in_year_min = [0]*len(years)

    # Plot
    ax.plot(years, values_mean_e, color='red', label='Ende der vegetativen Phase')
    ax.plot(years, values_mean_b, color='green', label='Begin der vegetativen Phase')

    # Marking the beginning of the seasons
    #ax.axhline(y=335, color='lightblue', linestyle='--', label='Winterbeginn')
    #ax.axhline(y=244, color='orange', linestyle='--', label='Herbstbeginn')
    #ax.axhline(y=152, color='darkgreen', linestyle='--', label='Sommerbeginn')
    #ax.axhline(y=60, color='lightgreen', linestyle='--', label='Frühlingsbeginn')

    # Fill between lines
    ax.fill_between(years, values_mean_e, values_mean_b, color='green', alpha=0.25, label='Vegetative Phase')

    # Gridlines:
    ax.grid(color='lightgrey', linewidth=0.5)

    # Plot Customization
    fig.set_size_inches(6.3*2, 3.15*2)
    ax.set_ylim([0, 365])
    ax.set_xlim([min(years), max(years)])
    #ax.set_title('Vegetative Phase')
    ax.set_xlabel('Jahre')
    ax.set_ylabel('Tage')
    ax.legend(loc='upper right')
    plt.xticks(rotation=45)
    for label in ax.xaxis.get_ticklabels():  # Iterate over all ticklabels
        if int(label.get_text()) % 5 == 0:   # Check if ticklabel is dividable by 5
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Add month ticks on right side
    ax2 = ax.twinx()
    month_days = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    month_labels = ["Jan", "Feb", "Mär", "Apr", "Mai", "Jun", 
                    "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(month_days)
    ax2.set_yticklabels(month_labels)
    ax2.set_ylabel("Monatsbeginn")

    # Save Plot
    plotname = 'vegetativ_phase'+'_plot.png'
    plot_folder_path = Path.cwd() / 'output'
    plot_folder_path.mkdir(parents=True, exist_ok=True)
    plot_path = str(plot_folder_path / plotname)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    print(f'Successfully created and saved plot: {plotname}')


# In[25]:


def plot_vegetation_phase_length():
    # Vegetation phase length
    plt.close()

    # Create a figure containing a single Axes.
    fig, ax = plt.subplots()

    # Vegetation begin line
    title, years, values_max, values_mean_b, values_min = years_values('vegetation_begin')

    # Vegetation end line
    title, years, values_max, values_mean_e, values_min = years_values('vegetation_end')

    # Vegetation phase length
    veg_len = [values_mean_e[i]-values_mean_b[i] for i in range(len(years))]

    # List containing 365 (days per year) as many times as there are years (upper limit)
    days_in_year_max = [365]*len(years)

    # List containing 365 (days per year) as many times as there are years (lower limit)
    days_in_year_min = [0]*len(years)

    # Plot
    ax.plot(years, veg_len, color='green', label='Vegetative Phase')

    # Fill between lines
    ax.fill_between(years, veg_len, days_in_year_min, color='green', alpha=0.25)

    # Gridlines:
    ax.grid(color='lightgrey', linewidth=0.5)

    # Plot Customization
    fig.set_size_inches(6.3*2, 3.15*2)
    ax.set_ylim([0, 365])
    ax.set_xlim([min(years), max(years)])
    #ax.set_title('Länge der vegetativen Phase')
    ax.set_xlabel('Jahre')
    ax.set_ylabel('Tage')
    ax.legend(loc='upper right')
    plt.xticks(rotation=45)
    for label in ax.xaxis.get_ticklabels():  # Iterate over all ticklabels
        if int(label.get_text()) % 5 == 0:   # Check if ticklabel is dividable by 5
            label.set_visible(True)
        else:
            label.set_visible(False)

    # Save Plot
    plotname = 'vegetativ_phase_len'+'_plot.png'
    plot_folder_path = Path.cwd() / 'output'
    plot_folder_path.mkdir(parents=True, exist_ok=True)
    plot_path = str(plot_folder_path / plotname)
    plt.savefig(plot_path, bbox_inches='tight', dpi=300)
    print(f'Successfully created and saved plot: {plotname}')


# # Run the Program

# In[26]:


# Get the shapefile to analyze

shp = get_shp()


# In[27]:


# Download the data

pdf_foldername = 'data_info'
print(f'\nDownload the PDF files containing informations about the used data:')
pdf_links = list_of_dwd_data(file_types=['.pdf'])
download_dwd_data(pdf_links, pdf_foldername)

raster_foldername = 'climate_environment_CDC_grids_germany_annual'
print(f'\nDownload the Rasterfiles:')
raster_links = list_of_dwd_data(file_types=['.asc.gz', '.zip'])
download_dwd_data(raster_links, raster_foldername)


# In[28]:


# Create JSON

print('\nProcess the Data:')

raster_path = str(Path.cwd() / raster_foldername)
prj_file = 'gk3.prj'

rasterstats_json, shp_crs_dissolved = zonal_climate_analysis(shp, raster_path, prj_file) # shp comes from get_shp() in the beginning of this program


# In[29]:


# Load Rasterstats JSON

input_file = rasterstats_json # Created with zonal_climate_analysis(shp_input, raster_folder, prj_file)

# Open rasterstats_dict.json file
with open(input_file) as json_file:
    rs = json.load(json_file)


# In[30]:


# Create Maps and Plots

print('\nCreating Map and Plots:')
matplotlib.use('Agg')  # Use a non-GUI backend. Prevents "QSocketNotifier: Can only be used with threads started with QThread" message in cmd.

# Create Map
create_map(shp_crs_dissolved)

# Create Plots
plot_air_temp_min_mean_max()
plot_frost_ice_days()
plot_snowcover_days()
plot_summer_hot_days()
plot_precipitaion()
plot_precipitaion_days()
plot_sunshine_duration()
plot_vegetation_begin_end()
plot_vegetation_phase_length()

print('\nFinished!')
print(f'Map and plots saved here: \n{Path.cwd() / 'output'}\n')

