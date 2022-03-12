# Author: Sayantan Majumdar
# Email: smxnv@mst.edu

import rasterio as rio
import numpy as np
import json
import os
import subprocess
import geopandas as gpd
from rasterio.mask import mask
from shapely.geometry import box
from fiona import transform
from glob import glob
from .sysops import make_gdal_sys_call_str


def read_raster_as_arr(raster_file, band=1, get_file=True, rasterio_obj=False, change_dtype=True):
    """
    Get raster array
    :param raster_file: Input raster file path
    :param band: Selected band to read (Default 1)
    :param get_file: Get rasterio object file if set to True
    :param rasterio_obj: Set true if raster_file is a rasterio object
    :param change_dtype: Change raster data type to float if true
    :return: Raster numpy array and rasterio object file (get_file=True and rasterio_obj=False)
    """

    if not rasterio_obj:
        raster_file = rio.open(raster_file)
    else:
        get_file = False
    raster_arr = raster_file.read(band)
    if change_dtype:
        raster_arr = raster_arr.astype(np.float32)
        if raster_file.nodata is not None:
            raster_arr[np.isclose(raster_arr, raster_file.nodata)] = np.nan
    if get_file:
        return raster_arr, raster_file
    return raster_arr


def write_raster(raster_data, raster_file, transform_, outfile_path, no_data_value, ref_file=None, out_crs=None):
    """
    Write raster file in GeoTIFF format
    :param raster_data: Raster data to be written
    :param raster_file: Original rasterio raster file containing geo-coordinates
    :param transform_: Affine transformation matrix
    :param outfile_path: Outfile file path
    :param no_data_value: No data value for raster (default float32 type is considered)
    :param ref_file: Write output raster considering parameters from reference raster file
    :param out_crs: Output crs
    :return: None
    """
    if ref_file:
        raster_file = rio.open(ref_file)
        transform_ = raster_file.transform
    crs = raster_file.crs
    if out_crs:
        crs = out_crs
    with rio.open(
            outfile_path,
            'w',
            driver='GTiff',
            height=raster_data.shape[0],
            width=raster_data.shape[1],
            dtype=raster_data.dtype,
            crs=crs,
            transform=transform_,
            count=raster_file.count,
            nodata=no_data_value
    ) as dst:
        dst.write(raster_data, raster_file.count)


def get_nodata():
    """
    Return fixed no data value for rasters
    :return: Pre-defined/hard-coded default no data value for the USGS MAP Project
    """

    return -32767


def crop_raster(input_raster_file, ref_file, output_raster_file):
    """
    Crop raster to bbox extents
    :param input_raster_file: Input raster file path
    :param ref_file: Reference raster or shape file to crop input_raster_file
    :param output_raster_file: Output raster file path
    :return: None
    """

    input_raster = rio.open(input_raster_file)
    if '.shp' in ref_file:
        ref_raster_ext = gpd.read_file(ref_file)
    else:
        ref_raster = rio.open(ref_file)
        minx, miny, maxx, maxy = get_raster_extent(ref_raster, is_rio_obj=True)
        ref_raster_ext = gpd.GeoDataFrame({'geometry': box(minx, miny, maxx, maxy)}, index=[0],
                                          crs=ref_raster.crs.to_string())
    ref_raster_ext = ref_raster_ext.to_crs(crs=input_raster.crs.data)
    ref_raster_ext.to_file('../Outputs/Test_Ref_Ext.shp')
    coords = [json.loads(ref_raster_ext.to_json())['features'][0]['geometry']]
    out_img, out_transform = mask(dataset=input_raster, shapes=coords, crop=True)
    out_img = out_img.squeeze()
    write_raster(out_img, input_raster, transform_=out_transform, outfile_path=output_raster_file,
                 no_data_value=input_raster.nodata)


def crop_rasters(input_raster_dir, input_mask_file, outdir, pattern='*.tif', prefix=''):
    """
    Crop multiple rasters in a directory
    :param input_raster_dir: Directory containing raster files which are named as *_<Year>.*
    :param input_mask_file: Mask file (shapefile) used for cropping
    :param outdir: Output directory for storing masked rasters
    :param pattern: Raster extension
    :param prefix: Output file prefix
    :return: None
    """

    for raster_file in glob(input_raster_dir + pattern):
        out_raster = outdir + prefix + raster_file[raster_file.rfind(os.sep) + 1:]
        crop_raster(raster_file, input_mask_file, out_raster)


def reproject_rasters(input_raster_dir, ref_raster, outdir, pattern='*.tif', gdal_path='/usr/bin/',
                      resampling_func='near', verbose=True):
    """
    Reproject rasters in a directory
    :param input_raster_dir: Directory containing raster files which are named as *_<Year>.*
    :param ref_raster: Reference raster file to consider while reprojecting
    :param outdir: Output directory for storing reprojected rasters
    :param pattern: Raster extension
    :param gdal_path: GDAL directory path, in Windows replace with OSGeo4W directory path, e.g. '/usr/bin/gdal/' on
    Linux or Mac and 'C:/OSGeo4W64/' on Windows, the '/' at the end is mandatory
    :param resampling_func: Resampling function ('near', 'bilinear', 'cubic', cubicspline', 'lanczos', 'average',
    'mode', 'max', 'min', 'med', 'q1', 'q3')
    :param verbose: Set True to print system call info
    :return: None
    """

    for raster_file in glob(input_raster_dir + pattern):
        out_raster = outdir + raster_file[raster_file.rfind(os.sep) + 1:]
        reproject_raster_gdal_syscall(
            raster_file,
            from_raster=ref_raster,
            outfile_path=out_raster,
            gdal_path=gdal_path,
            verbose=verbose,
            resampling_func=resampling_func
        )


def reproject_raster_gdal_syscall(input_raster_file, outfile_path, resampling_factor=1, resampling_func='near',
                                  downsampling=True, from_raster=None, keep_original=False, gdal_path='/usr/bin/',
                                  verbose=False, dst_xres=None, dst_yres=None):
    """
    Reproject raster using GDAL system call.
    :param input_raster_file: Input raster file
    :param outfile_path: Output file path
    :param resampling_factor: Resampling factor (default 3)
    :param resampling_func: Resampling function
    :param downsampling: Downsample raster (default True)
    :param from_raster: Reproject input raster considering another raster
    :param keep_original: Set True to only use the new projection system from 'from_raster'. The original raster extent
    is not changed
    :param gdal_path: GDAL directory path, in Windows replace with OSGeo4W directory path, e.g. '/usr/bin/gdal/' on
    Linux or Mac and 'C:/OSGeo4W64/' on Windows, the '/' at the end is mandatory
    :param verbose: Set True to print system call info
    :param dst_xres: Target xres in input_raster_file units. Set resampling_factor to None
    :param dst_yres: Target yres in input_raster_file units. Set resampling factor to None
    :return: None
    """

    src_raster_file = rio.open(input_raster_file)
    rfile = src_raster_file
    if from_raster and not keep_original:
        if isinstance(from_raster, str):
            rfile = rio.open(from_raster)
        else:
            rfile = from_raster
        resampling_factor = 1
    xres, yres = rfile.res
    extent = get_raster_extent(rfile, is_rio_obj=True)
    dst_proj = rfile.crs.to_string()
    no_data = src_raster_file.nodata
    if dst_xres and dst_yres:
        xres, yres = dst_xres, dst_yres
    elif resampling_factor:
        if not downsampling:
            resampling_factor = 1 / resampling_factor
        xres, yres = xres * resampling_factor, yres * resampling_factor
    args = ['-t_srs', dst_proj, '-te', str(extent[0]), str(extent[1]), str(extent[2]), str(extent[3]),
            '-dstnodata', str(no_data), '-r', str(resampling_func), '-tr', str(xres), str(yres), '-ot', 'Float32',
            '-overwrite', input_raster_file, outfile_path]
    sys_call = make_gdal_sys_call_str(gdal_path=gdal_path, gdal_command='gdalwarp', args=args, verbose=verbose)
    subprocess.call(sys_call)


def get_raster_extent(input_raster, new_crs=None, is_rio_obj=False):
    """
    Get raster extents using rasterio
    :param input_raster: Input raster file path
    :param new_crs: Specify a new crs to convert original extents
    :param is_rio_obj: Set True if input_raster is a Rasterio object and not a string
    :return: Raster extents in a list
    """

    if not is_rio_obj:
        input_raster = rio.open(input_raster)
    raster_extent = input_raster.bounds
    left = raster_extent.left
    bottom = raster_extent.bottom
    right = raster_extent.right
    top = raster_extent.top
    if new_crs:
        raster_crs = input_raster.crs.to_string()
        if raster_crs != new_crs:
            new_coords = reproject_coords(raster_crs, new_crs, [[left, bottom], [right, top]])
            left, bottom = new_coords[0]
            right, top = new_coords[1]
    return [left, bottom, right, top]


def reproject_coords(src_crs, dst_crs, coords):
    """
    Reproject coordinates. Copied from https://bit.ly/3mBtowB
    Author: user2856 (StackExchange user)
    :param src_crs: Source CRS
    :param dst_crs: Destination CRS
    :param coords: Coordinates as tuple of lists
    :return: Transformed coordinates as tuple of lists
    """

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    xs, ys = transform.transform(src_crs, dst_crs, xs, ys)
    return [[x, y] for x, y in zip(xs, ys)]
