# Author: Sayantan Majumdar
# Email: smxnv@mst.edu


import pandas as pd
import geopandas as gpd
import ee
import numpy as np
import os
import urllib.request
from .sysops import makedirs, make_proper_dir_name
from .rasterops import read_raster_as_arr
from .rasterops import reproject_rasters
from .vectorops import shp2raster


def get_gee_dict(get_key_list=False):
    """
    Get the available GEE data dictionary (global-scale data sets).
    :param get_key_list: Set True to get only the key list
    :return: GEE data, band, and scale dictionaries if get_key_list if False. Otherwise, only list of keys
    """

    gee_data_dict = {
        'SMOS_SMAP': ['NASA_USDA/HSL/soil_moisture', 'NASA_USDA/HSL/SMAP10KM_soil_moisture'],
        'LANDSAT_NDWI': ['LANDSAT/LE07/C01/T1_8DAY_NDWI', 'LANDSAT/LC08/C01/T1_8DAY_NDWI'],
        'LANDSAT_NDVI': ['LANDSAT/LE07/C01/T1_8DAY_NDVI', 'LANDSAT/LC08/C01/T1_8DAY_NDVI'],
        'GPM': 'NASA/GPM_L3/IMERG_MONTHLY_V06',
        'NASADEM': 'NASA/NASADEM_HGT/001',
        'DMSP_VIIRS': ['NOAA/DMSP-OLS/NIGHTTIME_LIGHTS', 'NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG'],
        'VIIRS_NDVI': 'NOAA/VIIRS/001/VNP13A1',
        'VIIRS_EVI': 'NOAA/VIIRS/001/VNP13A1',
        'VIIRS_EVI2': 'NOAA/VIIRS/001/VNP13A1',
        'MODIS_Day_LST': 'MODIS/006/MOD11A2',
        'MODIS_Night_LST': 'MODIS/006/MOD11A2',
        'MODIS_Google_NDVI': 'MODIS/MCD43A4_006_NDVI',
        'MODIS_Google_EVI': 'MODIS/MCD43A4_006_EVI',
        'MODIS_Terra_NDVI': 'MODIS/006/MOD13Q1',
        'MODIS_Terra_EVI': 'MODIS/006/MOD13Q1',
        'MODIS_NDWI': 'MODIS/MCD43A4_006_NDWI',
        'MODIS_BURNED_AREA': 'MODIS/006/MCD64A1',
        'MODIS_NPP': 'MODIS/006/MOD17A3HGF',
        'MODIS_FIRE': 'MODIS/006/MOD14A2',
        'MODIS_LAI': 'MODIS/006/MCD15A3H',
        'MODIS_ET': 'MODIS/006/MOD16A2',
        'MODIS_AEROSOL': 'MODIS/061/MOD08_M3',
        'MODIS_CIRRUS': 'MODIS/061/MOD08_M3',
        'MODIS_WV': 'MODIS/006/MCD19A2_GRANULES',
        'SW_Recurrence': 'JRC/GSW1_3/GlobalSurfaceWater',
        'SW_Occurence': 'JRC/GSW1_3/GlobalSurfaceWater',
        'SW_Seasonality': 'JRC/GSW1_3/GlobalSurfaceWater',
        'FAO_ACTUAL_ET': 'FAO/WAPOR/2/L1_AETI_D',
        'FAO_EVAPORATION': 'FAO/WAPOR/2/L1_E_D',
        'NASA_LANCE_FIRMS': 'FIRMS',
        'FLDAS_STORM_RO': 'NASA/FLDAS/NOAH01/C/GL/M/V001',
        'FLDAS_BF_GW_RO': 'NASA/FLDAS/NOAH01/C/GL/M/V001',
        'FLDAS_SM': 'NASA/FLDAS/NOAH01/C/GL/M/V001',
    }
    gee_band_dict = {
        'SMOS_SMAP': ['ssm', 'ssm'],
        'LANDSAT_NDWI': ['NDWI', 'NDWI'],
        'LANDSAT_NDVI': ['NDVI', 'NDVI'],
        'GPM': 'precipitation',
        'NASADEM': 'elevation',
        'DMSP_VIIRS': ['avg_vis', 'avg_rad'],
        'VIIRS_NDVI': 'NDVI',
        'VIIRS_EVI': 'EVI',
        'VIIRS_EVI2': 'EVI2',
        'MODIS_Day_LST': 'LST_Day_1km',
        'MODIS_Night_LST': 'LST_Night_1km',
        'MODIS_Google_NDVI': 'NDVI',
        'MODIS_Google_EVI': 'EVI',
        'MODIS_Terra_NDVI': 'NDVI',
        'MODIS_Terra_EVI': 'EVI',
        'MODIS_NDWI': 'NDWI',
        'MODIS_BURNED_AREA': 'BurnDate',
        'MODIS_NPP': 'Npp',
        'MODIS_FIRE': 'FireMask',
        'MODIS_LAI': 'Lai',
        'MODIS_ET': 'ET',
        'MODIS_AEROSOL': 'Aerosol_Optical_Depth_Land_Ocean_Mean_Mean',
        'MODIS_CIRRUS': 'Cirrus_Fraction_SWIR_FMean',
        'MODIS_WV': 'Column_WV',
        'SW_Recurrence': 'recurrence',
        'SW_Occurence': 'occurence',
        'SW_Seasonality': 'seasonality',
        'FAO_ACTUAL_ET': 'L1_AETI_D',
        'FAO_EVAPORATION': 'L1_E_D',
        'NASA_LANCE_FIRMS': 'T21',
        'FLDAS_STORM_RO': 'Qs_tavg',
        'FLDAS_BF_GW_RO': 'Qsb_tavg',
        'FLDAS_SM': 'SoilMoi00_10cm_tavg'
    }

    gee_scale_dict = {
        'SMOS_SMAP': [1, 1],
        'LANDSAT_NDWI': [1, 1],
        'LANDSAT_NDVI': [1, 1],
        'GPM': 1,
        'NASADEM': 1,
        'DMSP_VIIRS': [1, 1],
        'VIIRS_NDVI': 0.0001,
        'VIIRS_EVI': 0.0001,
        'VIIRS_EVI2': 0.0001,
        'MODIS_Day_LST': 0.02,
        'MODIS_Night_LST': 0.02,
        'MODIS_Google_NDVI': 1,
        'MODIS_Google_EVI': 1,
        'MODIS_Terra_NDVI': 0.0001,
        'MODIS_Terra_EVI': 0.0001,
        'MODIS_NDWI': 1,
        'MODIS_BURNED_AREA': 1,
        'MODIS_NPP': 0.0001,
        'MODIS_FIRE': 1,
        'MODIS_LAI': 0.1,
        'MODIS_ET': 0.1,
        'MODIS_AEROSOL': 0.001,
        'MODIS_CIRRUS': 0.0001,
        'MODIS_WV': 0.001,
        'SW_Recurrence': 1,
        'SW_Occurence': 1,
        'SW_Seasonality': 1,
        'FAO_ACTUAL_ET': 0.1,
        'FAO_EVAPORATION': 0.1,
        'NASA_LANCE_FIRMS': 1,
        'FLDAS_STORM_RO': 1,
        'FLDAS_BF_GW_RO': 1,
        'FLDAS_SM': 1
    }
    if not get_key_list:
        return gee_data_dict, gee_band_dict, gee_scale_dict
    return list(gee_data_dict.keys())


def save_gee_data(gee_data, gee_scale, gee_aoi, data_name, month, year, output_dir):
    """
    Save GEE data
    :param gee_data: EE Image object
    :param gee_scale: GEE scale in m
    :param gee_aoi: Area of Interest bounding box
    :param data_name: Name of the data set
    :param month: Month as int
    :param year: Year as int
    :param output_dir: Output directory
    :return: None
    """

    gee_url = gee_data.getDownloadUrl({
        'scale': gee_scale,
        'crs': 'EPSG:4326',
        'region': gee_aoi,
        'format': 'GEO_TIFF'
    })
    month_str = str(month)
    if month < 10:
        month_str = '0' + month_str
    local_file_name = output_dir + '{}_{}{}.tif'.format(data_name, month_str, year)
    print('Dowloading', local_file_name, '...')
    urllib.request.urlretrieve(gee_url, local_file_name)


def download_gee_data(year_list, start_month, end_month, outdir, data_extent, data='MODIS_ET', gee_scale=1000):
    """
    Download monthly GEE data. (Note: MOD16 has to be divided by 10 (line 38) as its original scale is 0.1 mm/8 days.)
    :param year_list: List of years in %Y format
    :param start_month: Start month in %m format
    :param end_month: End month in %m format
    :param outdir: Download directory
    :param data_extent: Data extent as a list in [minx, miny, maxx, maxy] format
    :param data: Name of the data set, MOD16 for MOD16 ET, SMOS_SMAP for SMOS/SMAP soil moisture, etc. The data name
    should be present in the GEE data dict keys
    :param gee_scale: GEE Data Scale in m
    :return: None
    """

    ee.Initialize()
    gee_data_dict, gee_band_dict, gee_scale_dict = get_gee_dict()
    gee_aoi = ee.Geometry.Rectangle(data_extent)
    multi_collection_dict = {
        'SMOS_SMAP': (2015, 4),
        'LANDSAT_NDWI': (2013, 4),
        'LANDSAT_NDVI': (2013, 4),
        'DMSP_VIIRS': (2014, 1)
    }
    multi_collection_data_list = list(multi_collection_dict.keys())
    if data in multi_collection_data_list:
        data_collection = ee.ImageCollection(gee_data_dict[data][0])
    elif data == 'NASADEM':
        data_collection = ee.Image(gee_data_dict[data])
        save_gee_data(data_collection, gee_scale, gee_aoi, data, 1, year_list[0], outdir)
        return
    else:
        data_collection = ee.ImageCollection(gee_data_dict[data])
    for year in year_list:
        band_name = gee_band_dict[data]
        band_scale = gee_scale_dict[data]
        for month in range(start_month, end_month + 1):
            try:
                start_date = ee.Date.fromYMD(year, month, 1)
                if month == 12:
                    end_date = ee.Date.fromYMD(year + 1, 1, 1)
                else:
                    end_date = ee.Date.fromYMD(year, month + 1, 1)
                if data in ['MODIS_ET', 'FAO_ACTUAL_ET', 'FAO_EVAPORATION', 'GPM']:
                    gee_data = data_collection.select(band_name).filterDate(start_date, end_date).sum()
                elif data in multi_collection_data_list:
                    year_check, month_check = multi_collection_dict[data]
                    if year >= year_check and month >= month_check:
                        data_collection = ee.ImageCollection(gee_data_dict[data][1])
                        band_name = gee_band_dict[data][1]
                        band_scale = gee_scale_dict[data][1]
                    else:
                        data_collection = ee.ImageCollection(gee_data_dict[data][0])
                        band_name = gee_band_dict[data][0]
                        band_scale = gee_scale_dict[data][0]
                    gee_data = data_collection.select(band_name).filterDate(start_date, end_date).mean()
                elif data in ['MODIS_Day_LST', 'MODIS_Night_LST']:
                    gee_data = data_collection.select(band_name).filterDate(start_date, end_date).median()
                else:
                    gee_data = data_collection.select(band_name).filterDate(start_date, end_date).mean()
                if band_scale != 1:
                    gee_data = gee_data.multiply(band_scale)
                save_gee_data(gee_data, gee_scale, gee_aoi, data, month, year, outdir)
            except ee.ee_exception.EEException:
                print('{} data not available for {}/{}'.format(data, month, year))


def prepare_data(input_shp, output_dir, data_list=('MODIS_ET', 'GPM'), data_start_year=2012, data_end_year=2021,
                 target_res=1000, gdal_path='/usr/bin', remove_na=False, already_prepared=False):
    """
    Prepare data for a region of interest
    :param input_shp: Input shapefile of the country containing administrative boundaries
    :param output_dir: Output directory to store files
    :param data_list: List of data sets to download, valid names include 'SM_IDAHO', 'MOD16', 'SMOS_SMAP',
    'DROUGHT', 'PRISM', 'TMIN', 'TMAX', 'WS', 'RO', 'NDWI', 'SPH', 'DEF', 'VPD', 'VPD_SMAP', 'MODIS'.
    :param data_start_year: Start year
    :param data_end_year: End year
    :param target_res: Target resolution for rasters in metres
    :param gdal_path: Path to gdal system path
    :param remove_na: Set True to remove NA values.
    :param already_prepared: Set True if subset csv already exists
    :return: Pandas data frame containing all the pixelwise data associated with each administrative boundary present
    in input_shp
    """

    output_csv = output_dir + 'NASA_ROSES_Data.csv'
    if not already_prepared:
        makedirs([make_proper_dir_name(output_dir)])
        input_gdf = gpd.read_file(input_shp)
        input_gdf['idx'] = input_gdf.index + 1
        modified_input_shp = output_dir + input_shp[input_shp.rfind('/') + 1:]
        input_gdf.to_file(modified_input_shp)
        xres = yres = target_res / 111e+3
        admin_raster = output_dir + input_shp[input_shp.rfind('/') + 1: input_shp.rfind('.')] + '.tif'
        shp2raster(
            modified_input_shp,
            admin_raster,
            value_field='idx',
            xres=xres, yres=yres,
            output_crs='EPSG:4326'
        )
        data_extent = input_gdf.total_bounds.tolist()
        year_list = range(data_start_year, data_end_year + 1)
        gee_file_dir = make_proper_dir_name(output_dir + 'GEE_Files')
        reproj_dir = make_proper_dir_name(gee_file_dir + 'Reprojected')
        makedirs([gee_file_dir, reproj_dir])
        if 'All' in data_list:
            data_list = get_gee_dict(get_key_list=True)
        for data in data_list:
            download_gee_data(year_list, 1, 12, gee_file_dir, data_extent, data, target_res)
        reproject_rasters(
            gee_file_dir,
            ref_raster=admin_raster,
            outdir=reproj_dir,
            gdal_path=gdal_path
        )
        output_df = pd.DataFrame()
        output_dict = {}
        admin_raster_arr = read_raster_as_arr(admin_raster, get_file=False)
        print('Creating CSV...')
        for data in data_list:
            if data != 'NASADEM':
                for year in year_list:
                    for month in range(1, 13):
                        month_str = str(month)
                        if month < 10:
                            month_str = '0' + month_str
                        raster_file = reproj_dir + '{}_{}{}.tif'.format(data, month_str, year)
                        if os.path.exists(raster_file):
                            raster_arr = read_raster_as_arr(raster_file, get_file=False)
                            if data == 'DMSP_VIIRS' and year >= 2014:
                                raster_arr = np.ceil(
                                    6.5 + 57.4 * (1 / (1 + np.exp(1.9 * (np.log(raster_arr + 1) - 10.8))))
                                )
                        else:
                            raster_arr = np.full_like(admin_raster_arr, fill_value=np.nan)
                        raster_val_list = raster_arr.ravel().tolist()
                        output_dict[data] = raster_val_list
                        raster_val_list = raster_arr.ravel().tolist()
                        output_dict['YEAR'] = [year] * len(raster_val_list)
                        output_dict['MONTH'] = [month] * len(raster_val_list)
                        if output_df.empty:
                            output_df = pd.DataFrame(data=output_dict)
                        else:
                            output_df = output_df.append(
                                pd.DataFrame(data=output_dict)
                            )
        num_periods = (year_list[-1] - year_list[0] + 1) * 12
        if 'NASADEM' in data_list:
            raster_file = reproj_dir + 'NASADEM_012012.tif'
            raster_arr = read_raster_as_arr(raster_file, get_file=False)
            output_df['NASADEM'] = raster_arr.ravel().tolist() * num_periods
        output_df['idx'] = admin_raster_arr.ravel().tolist() * num_periods
        nan_values = [np.inf, -np.inf]
        if remove_na:
            nan_values.append(np.nan)
        output_df = output_df[~output_df.isin(nan_values).any(1)]
        input_gdf = input_gdf.drop(columns=['geometry'])
        output_df = output_df.merge(input_gdf, on='idx', how='inner')
        output_df = reindex_df(output_df, ordering=True)
        output_df.to_csv(output_csv, index=False)
    else:
        output_df = pd.read_csv(output_csv)
    return output_df


def reindex_df(df, column_names=None, ordering=False):
    """
    Reindex dataframe columns
    :param df: Input dataframe
    :param column_names: Dataframe column names, these must be df headers
    :param ordering: Set True to apply ordering
    :return: Reindexed dataframe
    """
    if not column_names:
        column_names = df.columns
        ordering = True
    if ordering:
        column_names = sorted(column_names)
    return df.reindex(column_names, axis=1)
