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
git clone https://github.com/<you>/lokal-climate-analysis.git
cd lokal-climate-analysis

# 2. install Python ≥3.9 deps (recommended: venv/conda)
pip install -r requirements.txt

# 3. run the script
python lokal_climate_analysis.py
# follow the prompt:
# Enter the path to the shapefile: /path/to/my_area.shp
