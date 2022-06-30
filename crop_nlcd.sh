# Crop NLCD raster to the dimensions of Upper Oconee watershed

module load gdal
cd /project/qdr/spatial

# Get projection of NLCD raster, write it to a file, and transform the Upper Oconee boundary shapefile to that projection so they match
# gdalsrsinfo nlcd_2019_land_cover_l48_20210604.img -o wkt > nlcd_projection.wkt
# gdaltransform -t_srs nlcd_projection.wkt upper_oconee.gpkg upper_oconee_AEA.gpkg

gdalwarp -crop_to_cutline -cutline upper_oconee_AEA.gpkg -of GTiff nlcd_2019_land_cover_l48_20210604.img nlcd_land_cover_upper_oconee.tif
gdalwarp -crop_to_cutline -cutline upper_oconee_AEA.gpkg -of GTiff nlcd_2019_impervious_l48_20210604.img nlcd_impervious_upper_oconee.tif
