# Author: Sayantan Majumdar
# Email: smxnv@mst.edu
# NASA ROSES Proposal Codes
# PI: Molly Lipscomb, Venkataraman Lakshmi

# --------------------------------------- Running nasa-roses-uva ------------------------------------------------------

# Execute map_ml.py (change the paths and flags accordingly) on Linux or Mac
# python uva_roses.py \
# --input-shp '../Data/gis shapefiles/nga_admbnda_adm2_osgof_20190417.shp' \
# --load-files True \
# --output-dir ../Outputs/ \
# --start-year 2012 \
# --end-year 2021 \
# --data-list MOD16 SM_IDAHO SMOS_SMAP RO DEF GPM \
# --target-res 1000 \
# --gdal-path /usr/bin/gdal/ \
# --skip-download True \
# --remove-na False \
# --use-hpc True

# Execute map_ml.py (change the paths and flags accordingly) on Windows powershell
# python uva_roses.py `
# --input-shp '../Data/gis shapefiles/nga_admbnda_adm2_osgof_20190417.shp' `
# --load-files False `
# --output-dir ../Outputs/ `
# --start-year 2012 `
# --end-year 2021 `
# --data-list All `
# --target-res 1000 `
# --gdal-path C:/OSGeo4W64/ `
# --skip-download True `
# --remove-na False `
# --use-hpc True

# ------------------------------------------------- Main code begins --------------------------------------------------


import argparse
from uvalibs.dataops import prepare_data
from uvalibs.sysops import boolean_string


def run_map_ml(args):
    """
    Driver function to prepare data.
    :param args: Namespace containing the following keys:

    Keys: Description
    ___________________________________________________________________________________________________________________
    input_shp: Input shapefile
    load_files: Set True to load existing CSV
    skip_download: Set True to load existing GEE data
    output_dir: Output directory
    start_year: Start year for downloading data sets
    end_year: End year for downloading data sets
    data_list: List of data sets to use/download. Valid names include 'GPM', 'MODIS_ET', 'SMOS_SMAP', 'MODIS_NDWI', etc.
    Set All to use all data sets defined in this project
    target_res: Target resolution (m) for rasters
    gdal_path: GDAL path
    remove_na: Set True to remove NA values from the final CSV
    use_hpc: Set False to run on local machine
    ___________________________________________________________________________________________________________________
    :return: None
    """

    prepare_data(
        args.input_shp,
        args.output_dir,
        data_list=args.data_list,
        data_start_year=args.start_year,
        data_end_year=args.end_year,
        already_prepared=args.load_files,
        target_res=args.target_res,
        gdal_path=args.gdal_path,
        remove_na=args.remove_na,
        skip_download=args.skip_download,
        use_hpc=args.use_hpc
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Flags to run HydroMAP_ML')
    parser.add_argument('--input-shp', type=str, required=True, help='Input shapefile')
    parser.add_argument('--load-files', type=boolean_string, default=False,
                        help='Set True to load existing CSV')
    parser.add_argument('--skip-download', type=boolean_string, default=True,
                        help='Set True to load existing GEE data')
    parser.add_argument('--output-dir', type=str, required=True, help='Output directory')
    parser.add_argument('--start-year', type=int, required=True, help='Start year for downloading data sets')
    parser.add_argument('--end-year', type=int, required=True, help='End year for downloading data sets')
    parser.add_argument('--data-list', type=str, nargs='+', required=True,
                        help="List of data sets to use/download. Valid names include 'GPM', 'MODIS_ET',"
                             " 'SMOS_SMAP', 'MODIS_NDWI', etc. Set All to use all data sets defined in this project")
    parser.add_argument('--target-res', type=int, default=1000, help='Target resolution (m) for rasters')
    parser.add_argument('--gdal-path', type=str, help='GDAL path')
    parser.add_argument('--remove-na', type=boolean_string, default=False,
                        help='Set True to remove NA values from the final CSV')
    parser.add_argument('--use-hpc', type=boolean_string, default=True, help='Set False to run on local machine')
    map_args = parser.parse_args()
    run_map_ml(map_args)
