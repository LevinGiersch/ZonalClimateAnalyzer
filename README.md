# Lokal Climate Analysis  
*Analyze the climate history of any area inside Germany with nothing more than a shapefile.*

---

## 1 | Project Purpose
This script automates an end‑to‑end workflow to  

1. **Download** annual gridded climate rasters (1951 → latest) from the **DWD Climate Data Center (CDC)**.  
2. **Clip & analyse** those rasters for a user‑supplied area (polygon shapefile).  
3. **Summarise** the results as zonal statistics (min / mean / max) in a tidy JSON file.  
4. **Visualise** long‑term trends with ready‑made PNG plots.

The goal is to give municipalities, researchers and students a quick way to quantify and visualise local climate change indicators without manual GIS work.

---

## 2 | Key Features
| Stage | What happens | Where in code |
|-------|--------------|---------------|
| **Input** | Interactive prompt asks for your *.shp* path | `get_shp()` |
| **Data download** | All relevant **.asc.gz / .zip** rasters fetched from the DWD open‑data mirror | `download_dwd_data()` |
| **Pre‑processing** | Decompress → rename → add CRS → re‑encode to GeoTIFF | `decompress_file()`, `asc_to_tif_add_crs()` |
| **Clean‑up** | Source archives removed to keep disk footprint small | `delete_raster_files()` |
| **Analysis** | Per‑polygon zonal stats with [`rasterstats`](https://pythonhosted.org/rasterstats/) | `calculate_zonal_stats()` |
| **Output** | `area_rasterstats.json` plus 9 publication‑ready plots in */plots* | visualiser section |


---

## 3 | Quick‑start
```bash
# 1. clone & enter
git clone https://github.com/LevinGiersch/ZonalClimateAnalyzer
cd ZonalClimateAnalyzer

# 2. create venv
python3 -m venv venv
source venv/bin/activate

# 3. install Python ≥3.9 deps (recommended: venv/conda)
pip install -r requirements.txt

# 4. run the script
python lokal_climate_analysis.py
# follow the prompt:
# Enter the path to the shapefile: /path/to/my_area.shp

# 5. wait until the process is finished
```

Runtime hint: the first execution downloads ≈ 1 GB of rasters and can take 10–20 min (depending on your connection). Subsequent runs use the cached GeoTIFFs and finish in seconds.
## 4 | Required Python Packages

    geopandas ‒ vector I/O & reprojection

    rasterio ‒ raster I/O

    rasterstats ‒ zonal statistics

    requests + beautifulsoup4 ‒ HTML scraping

    tqdm ‒ progress bars

    matplotlib ‒ plotting

    pyproj, filetype, numpy, pandas (implicit)

Install them via pip install -r requirements.txt or adapt to your environment (e.g. conda‑forge).


## 5 | Data Source & Coordinate System
| Aspect | Details |
|--------|---------|
| **Provider** | Deutscher Wetterdienst (DWD) – Climate Data Center |
| **URL root** | <https://opendata.dwd.de/climate_environment/CDC/grids_germany/annual/> |
| **Parameters pulled** | air_temperature\_\*, frost_days, hot_days, ice_days, drought_index, precipitation, snowcover_days, precipGE{10,20,30}mm\_days, sunshine_duration, vegetation\_{begin,end} |
| **Spatial grid** | 1 × 1 km **GK3 / DHDN Zone 3** (EPSG 31467) |
| **Temporal coverage** | 1951 – present (updated yearly) |

A small `gk3.prj` file ships with the repo; shapefiles are re‑projected into this CRS so raster overlays line up exactly.

---

## 6 | Outputs Explained

| File | What it shows |
|------|---------------|
| `min_mean_max_temp_plot.png` | Annual **maximum, mean and minimum air temperature**. Filled bands visualise the spread between the three series. |
| `ice_frost_days_plot.png` | Counts of **frost days** (T<sub>min</sub> < 0 °C) and **ice days** (T<sub>max</sub> < 0 °C). |
| `snowcover_days_plot.png` | **Days per year with snow depth > 1 cm**. |
| `summer_hot_days_plot.png` | Counts of **summer days** (T<sub>max</sub> ≥ 25 °C) and **hot days** (T<sub>max</sub> ≥ 30 °C). |
| `precipitation_drought_plot.png` | **Annual precipitation totals** (bars) with **drought index** overlay (line, scaled). |
| `precip_days_plot.png` | **Number of days with heavy precipitation** ≥ 10 mm, ≥ 20 mm, ≥ 30 mm. |
| `sunshine_duration_plot.png` | **Average daily sunshine hours** for each year. |
| `vegetativ_phase_plot.png` | **Start and end dates of the vegetative season** plus reference lines for astronomical seasons. |
| `vegetativ_phase_len_plot.png` | **Length of the vegetative season** (days between start and end) for each year. |


## 7 | Extending the Script
- **Add more parameters** – append new sub‑folders to `folder_download_locations`.  
- **Different country** – swap the DWD URLs for your own open‑grid source and adjust the CRS.  
- **Batch mode** – loop `lokal_climate_analysis()` over multiple shapefiles.  
- **Package it** – move the functions into `src/` and expose a CLI entry‑point.  
Contributions are welcome – open an issue or PR!

---

## 8 | Known Limitations
- Designed for **Germany‑wide 1 km grids** only.  
- Memory‑intensive if the shapefile contains *many* polygons.  
- Re‑downloads rasters each year; archive copies yourself for full reproducibility.  

---

## 9 | License
Include an OSI‑approved license here (e.g. **MIT**).  
Always credit the **DWD Climate Data Center** when publishing derived work.

---

## 10 | Citation
If you use this tool in an academic context, please cite:

> Deutscher Wetterdienst (2025): *Grids Germany – Annual*.  
https://opendata.dwd.de/climate_environment/CDC/grids_germany/annual/

---
