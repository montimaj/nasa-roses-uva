# NASA ROSES UVA Proposal Codes 
Codes for remote sensing data acquisition and preprocessing

Author: Sayantan Majumdar

Email: smxnv@mst.edu

PI: Molly Lipscomb, Venkataraman Lakshmi

Other member(s): Brendan Novak

## Running nasa-roses-uva

### <i>Setup the Conda environment</i>
```
conda env create -f environment.yml
```

### <i>Execute uva_roses.py (change the paths and flags accordingly) on Linux or Mac </i>
```
python uva_roses.py \
--input-shp '../Data/gis shapefiles/nga_admbnda_adm2_osgof_20190417.shp' \
--load-files True \
--output-dir ../Outputs/ \
--start-year 2012 \
--end-year 2021 \
--data-list MODIS_ET SMOS_SMAP GPM LANDSAT_NDWI \
--target-res 1000 \
--gdal-path /usr/bin/gdal/ \
--skip-download True \
--remove-na False \
--use-hpc False \
--num-chunks 100
```

### <i>Execute uva_roses.py (change the paths and flags accordingly) on Windows powershell</i>
```
python uva_roses.py `
--input-shp '../Data/gis shapefiles/nga_admbnda_adm2_osgof_20190417.shp' `
--load-files False `
--output-dir ../Outputs/ `
--start-year 2012 `
--end-year 2021 `
--data-list All `
--target-res 1000 `
--gdal-path C:/OSGeo4W64/ `
--skip-download True `
--remove-na False `
--use-hpc False `
--num-chunks 100
```

## Usage Notes

Once the data sets are downloaded, the program reprojects them based on the administrative boundary shapefile 
(converted to a raster in the pipeline). All the rasters are in 1 km x 1 km grid (0.09 deg x 0.09 deg, EPSG:4326).
The CSV files are created individually for each raster in a Dask parallel processing pipeline. In order to get the 
county index associated with each pixel, the raster CSV files have to be inner-joined with the admintration boundary CSV
(created during raster CSV generation) based on the 'idx' fields.