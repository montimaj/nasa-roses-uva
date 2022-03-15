# NASA ROSES UVA Proposal Codes 
Codes for remote sensing data acquisition and preprocessing

Author: Sayantan Majumdar

Email: smxnv@mst.edu

PI: Molly Lipscomb, Venkataraman Lakshmi

Other member(s): Brendan Novak

## Running nasa-roses-uva

Execute uva_roses.py (change the paths and flags accordingly) on Linux or Mac
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
--use-hpc True \
--num-chunks 100
```

Execute uva_roses.py (change the paths and flags accordingly) on Windows powershell
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
--use-hpc True `
--num-chunks 100
```