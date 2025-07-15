# ZonalClimateAnalyzer
This Program will make it easy to look at the climate history of any Place in Germany. All you need is this GitHub repository and a shapefile (filename.shp) with a corresponding projection file (filename.prj). This Program will download climate information for every point in germany from the DWD (https://opendata.dwd.de/climate_environment/CDC/grids_germany/annual/), calculte the climate parameters for your shapefile and return some nice plots 

### Prerequesites
- Stable internet Connection
- Approx.  7GB of available diskspace
- All python packages in the pip_installer file


### How this Program works:
1. Asks you for a shapefile
2. Creates a new python venv
3. Installs all needed packages in the new python venv
4. Downloads Raster data containing Climate information of Germany
5. Calculates climate properties for the provided shapefile
6. Returns a folder full of Plots and Maps showing the climate history for your shapefile from the years 1901 to 2024
7. Removes the downloaded rasterfiles
8. Optional: you can use the provided Latex-File to automatically create a PDF containing all created Plots and Maps