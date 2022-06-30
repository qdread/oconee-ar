# Download data and get ready to make maps


# Site locations ----------------------------------------------------------


library(sf)
library(readxl)
library(dplyr)
library(stars)

site_locs <- read_xlsx('project/AR paper_sampling sites_location.xlsx')
names(site_locs) <- c('site', 'lat', 'lon')

site_sp <- st_as_sf(site_locs, coords = c('lon', 'lat'), crs = st_crs('+proj=longlat')) %>%
  mutate(area = substr(site, 1, 4))

st_write(site_sp, 'project/sites.gpkg', driver = 'GPKG')



# Boundary of UO watershed ------------------------------------------------

# This will be used as a template to get any other needed data.
# https://www.usgs.gov/national-hydrography/access-national-hydrography-products
# It is code 03070101 see https://water.usgs.gov/lookup/getwatershed?03070101
# Download page https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/WBD/HU2/GPKG/
# File https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/HU2/GPKG/WBD_03_HU2_GPKG.zip

# Read in the watershed boundaries dataset for HUC2 unit 03 (southeastern US)
st_layers('~/spatial_data/WBD_03_HU2_GPKG.gpkg')
huc8s <- st_read(dsn = '~/spatial_data/WBD_03_HU2_GPKG.gpkg', layer = 'WBDHU8')
upper_oconee <- huc8s %>% filter(huc8 == '03070101')

st_write(upper_oconee, 'project/upper_oconee.gpkg', driver = 'GPKG')


# NLCD data for UO watershed ----------------------------------------------

# Page https://www.mrlc.gov/data?f%5B0%5D=year%3A2019
# NLCD 2019 landcover data at https://s3-us-west-2.amazonaws.com/mrlc/nlcd_2019_land_cover_l48_20210604.zip
# NLCD 2019 percent imperviousness data at https://s3-us-west-2.amazonaws.com/mrlc/nlcd_2019_impervious_l48_20210604.zip

nlcd <- read_stars('~/spatial_data/nlcd_2019_land_cover_l48_20210604.img')
upper_oconee <- st_read('project/upper_oconee.gpkg')
upper_oconee_AEA <- st_transform(upper_oconee, crs = st_crs(nlcd))
st_write(upper_oconee_AEA, 'project/upper_oconee_AEA.gpkg', driver = 'GPKG')

gdalwarp(srcfile = '~/spatial_data/nlcd_2019_impervious_l48_20210604.img', dstfile = 'project/nlcd_upper_oconee.tif', of = 'GTiff', 
         crop_to_cutline = TRUE, cutline = 'project/upper_oconee_AEA.gpkg', dryrun=T)

nlcd_upper_oconee <- st_crop(nlcd, upper_oconee_AEA)

nlcd_imperv <- read_stars('~/spatial_data/nlcd_2019_impervious_l48_20210604.img')
nlcd_imperv_upper_oconee <- st_crop(nlcd_imperv, upper_oconee_AEA)

write_stars(nlcd_upper_oconee, 'project/nlcd_2019_upper_oconee.tif', driver = 'GTiff')
write_stars(nlcd_imperv_upper_oconee, 'project/nlcd_imperv_2019_upper_oconee.tif', driver = 'GTiff')

# River and stream paths in UO watershed ----------------------------------

# https://www.usgs.gov/national-hydrography/access-national-hydrography-products
# Downloaded from https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/NHD/HU8/GPKG/
# File https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/GPKG/NHD_H_03070101_HU8_GPKG.zip