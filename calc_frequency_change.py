"""
Script to calculate TC frequency changes and merge with mangrove centroids

Author: Sarah Hülsen

Note: The calculations in this script are computationally expensive and were originally run on ETH Zürich computation cluster Euler.
"""
# Calculate TC hazard frequency changes:
## - calculate TC frequency per defined intensity category for historical and future period
## - calculate changemaps
## - merge with mangrove data


from pathlib import Path
import numpy as np
import geopandas as gpd 
import copy
from shapely.geometry import box
import pandas as pd
from climada.hazard import Hazard
from climada.engine import ImpactCalc

data_dir = Path('../data/')

# Eckert IV equal area projection
proj_eck4 = '+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'  
# define hazard categories
# Low: Saffir-Simpson Cat1, Cat2 (33-49)
# Medium: Cat3, Cat4 (50-70)
# High: Cat5 and above (70 and above)
intensity_cat = {
    'tc': [0, 33, 50, 70, 300],
    }

resolution = '0150'
models = ['miroc6', 'cesm2', 'ecearth6']
basins = ['AP', 'IO', 'SH', 'WP']
ssps = ['245', '370', '585']
ecosystem_file = 'Bunting_mangroves_2020_global_rslr_all_ssps_2020-2100.csv'
source = 'Bunting'


# define functions to create frequency maps
def haz_to_cathaz(haz, cat):
    """Maps windspeed intensity to defined categories"""
    haz_cat = copy.deepcopy(haz)
    haz_cat.intensity.data = np.digitize(haz.intensity.data, cat) - 1
    haz_cat.intensity.eliminate_zeros()
    return haz_cat

def _make_freq_gdf(freqs_dict, centroids):
    """Creates geodataframe from calculated frequencies per given category and centroids"""
    gdf_freqs = gpd.GeoDataFrame(freqs_dict)
    if len(centroids.geometry) < centroids.size:
        centroids.set_geometry_points()
    gdf_freqs['geometry'] = centroids.geometry
    gdf_freqs.set_geometry('geometry')
    return gdf_freqs

def freqs_per_cat_gdf(haz_cat):
    """Calculates frequency per given category"""
    cat = np.unique(haz_cat.intensity.data).astype(int)
    haz_type = haz_cat.haz_type
    haz_freqs_dict = {}
    for c in cat:
        cat_freq = ImpactCalc.eai_exp_from_mat((haz_cat.intensity == c).astype(int), haz_cat.frequency)
        cat_freq = np.squeeze(np.asarray(cat_freq))
        haz_freqs_dict[haz_type + str(c)] = cat_freq
    return _make_freq_gdf(haz_freqs_dict, haz_cat.centroids)

# load ecosystem data
exp = pd.read_csv(data_dir.joinpath('mangroves', ecosystem_file))
exp_gdf = gpd.GeoDataFrame(data=exp, geometry=gpd.points_from_xy(exp["longitude"], exp["latitude"]), crs="EPSG:4326")
exp_bounds = exp_gdf.total_bounds
# clip hazard data data to ecosystem bounds
bounds_buffered = box(*exp_bounds).buffer(2) # add 2 degree buffer to the bounds
buffered_gdf = gpd.GeoDataFrame(geometry=[bounds_buffered], crs=exp_gdf.crs)
# reproject to equal area (Eckert IV)
exp_gdf = exp_gdf.to_crs(proj_eck4)


# create frequency maps for current climate & climate change scenarios
for model in models:
    # current climate
    for basin in basins:
            # create current climate frequency map per basin
            haz = Hazard.from_hdf5(data_dir.joinpath( f"TC_{basin}_{resolution}as_MIT_{model}_20thcal_H08.hdf5"))
            tc_cat = haz_to_cathaz(haz, intensity_cat['tc'])
            tc_freqs = freqs_per_cat_gdf(tc_cat)
            tc_freqs['lon_cent'] = tc_freqs.geometry.x
            tc_freqs['lat_cent'] = tc_freqs.geometry.y
            basin_poly = gpd.read_file(data_dir.joinpath('basin_polys', f'{basin}_poly.shp'))
            # clip frequency map to given basin
            tc_freqs_clip = gpd.clip(tc_freqs, basin_poly)
            tc_freqs_clip.to_csv(data_dir.joinpath("frequency_maps", "current", f"TC_{basin}_{resolution}as_MIT_{model}_20thcal_freqmap.csv"), index=False)
    # merge basins to create global frequency map
    freq_AP = pd.read_csv(data_dir.joinpath("frequency_maps", "current", f"TC_AP_{resolution}as_MIT_{model}_20thcal_freqmap.csv"))
    freq_IO = pd.read_csv(data_dir.joinpath("frequency_maps", "current", f"TC_IO_{resolution}as_MIT_{model}_20thcal_freqmap.csv"))
    freq_SH = pd.read_csv(data_dir.joinpath("frequency_maps", "current", f"TC_SH_{resolution}as_MIT_{model}_20thcal_freqmap.csv"))
    freq_WP = pd.read_csv(data_dir.joinpath("frequency_maps", "current", f"TC_WP_{resolution}as_MIT_{model}_20thcal_freqmap.csv"))
    freq_global = pd.concat([freq_AP, freq_IO, freq_SH, freq_WP])
    freq_global.to_csv(data_dir.joinpath("frequency_maps", "current", f"TC_global_{resolution}as_MIT_{model}_20thcal_freqmap.csv"), index=False)
    
    #future climate
    for ssp in ssps:
        for basin in basins:
            # create climate change frequency map per basin
            haz = Hazard.from_hdf5(haz_data_dir.joinpath(f"TC_{basin}_{resolution}as_MIT_{model}_ssp{ssp}_2cal_H08.hdf5"))
            tc_cat = haz_to_cathaz(haz, intensity_cat['tc'])
            tc_freqs = freqs_per_cat_gdf(tc_cat)
            tc_freqs['lon_cent'] = tc_freqs.geometry.x
            tc_freqs['lat_cent'] = tc_freqs.geometry.y
            basin_poly = gpd.read_file(data_dir.joinpath('basin_polys', f'{basin}_poly.shp'))
            # clip frequency map to given basin
            tc_freqs_clip = gpd.clip(tc_freqs, basin_poly)
            tc_freqs_clip.to_csv(data_dir.joinpath("frequency_maps", "future", f"TC_{basin}_{resolution}as_MIT_{model}_ssp{ssp}_2cal_freqmap.csv"), index=False)

        # merge basins to create global frequency map
        freq_AP = pd.read_csv(data_dir.joinpath('frequency_maps', 'future', f"TC_AP_{resolution}as_MIT_{model}_ssp{ssp}_2cal_freqmap.csv"))
        freq_IO = pd.read_csv(data_dir.joinpath('frequency_maps', 'future', f"TC_IO_{resolution}as_MIT_{model}_ssp{ssp}_2cal_freqmap.csv"))
        freq_SH = pd.read_csv(data_dir.joinpath('frequency_maps', 'future', f"TC_SH_{resolution}as_MIT_{model}_ssp{ssp}_2cal_freqmap.csv"))
        freq_WP = pd.read_csv(data_dir.joinpath('frequency_maps', 'future', f"TC_WP_{resolution}as_MIT_{model}_ssp{ssp}_2cal_freqmap.csv"))
        freq_global = pd.concat([freq_AP, freq_IO, freq_SH, freq_WP])
        freq_global.to_csv(data_dir.joinpath('frequency_maps', 'future', f"TC_global_{resolution}as_MIT_{model}_ssp{ssp}_2cal_freqmap.csv"), index=False)


# change maps per scenario and model
for model in models:
    current = pd.read_csv(data_dir.joinpath("frequency_maps", "current", f"TC_global_{resolution}as_MIT_{model}_20thcal_freqmap.csv"))

    for ssp in ssps:
        future = pd.read_csv(data_dir.joinpath('frequency_maps', 'future', f'TC_global_{resolution}as_MIT_{model}_ssp{ssp}_2cal_freqmap.csv'))
        change_df = pd.merge(current, future, on=['lat_cent', 'lon_cent', 'geometry'], suffixes=('_current', '_future'))
        # Calculate the change for each TC category column
        for cat in range(1, 4):
            current_col = f'TC{cat}_current'
            future_col = f'TC{cat}_future'
            change_col = f'TC{cat}_change'
            ratio_col = f'TC{cat}_ratio'
            change_df[change_col] = change_df[future_col] - change_df[current_col]
            change_df[ratio_col] = (change_df[future_col] - change_df[current_col]) / change_df[current_col]
        change_df.to_csv(data_dir.joinpath('frequency_maps', 'change_maps', f'TC_global_{resolution}as_MIT_{model}_ssp{ssp}_2cal_changemap.csv'), index=False)
        
        # combine with mangrove data
        change_df = pd.read_csv(data_dir.joinpath('frequency_maps', 'change_maps', f'TC_global_{resolution}as_MIT_{model}_ssp{ssp}_2cal_changemap.csv'))
        change_df['geometry'] = gpd.GeoSeries.from_wkt(change_df['geometry'])
        change_gdf = gpd.GeoDataFrame(change_df, geometry='geometry', crs="EPSG:4326")
        change_clip = gpd.clip(change_gdf, buffered_gdf)
        # reproject to equal area (Eckert IV)
        change_clip = change_clip.to_crs(proj_eck4)
        # match points & save
        join = gpd.sjoin_nearest(exp_gdf, change_clip, distance_col="TC_distance")
        join = join.to_crs("EPSG:4326")
        join.to_csv(data_dir.joinpath('mangroves', 'frequency_maps', f'{source}_mangroves_TC_global_{resolution}as_MIT_{model}_ssp{ssp}_2cal_changemap.csv'), index=False)


# calculate model median frequency & median change
for ssp in ssps:
    miroc6 = pd.read_csv(data_dir.joinpath('frequency_maps', 'change_maps', f"TC_global_{resolution}as_MIT_miroc6_ssp{ssp}_2cal_changemap.csv"))
    cesm2 = pd.read_csv(data_dir.joinpath('frequency_maps', 'change_maps', f"TC_global_{resolution}as_MIT_cesm2_ssp{ssp}_2cal_changemap.csv"))
    ecearth6 = pd.read_csv(data_dir.joinpath('frequency_maps', 'change_maps', f"TC_global_{resolution}as_MIT_ecearth6_ssp{ssp}_2cal_changemap.csv"))
    
    # merge model frequency maps
    dfs = [miroc6, cesm2, ecearth6]
    df_all = pd.concat(dfs)
    df_median = df_all.groupby(['lon_cent', 'lat_cent', 'geometry'],
                                as_index=False)[['TC1_current',
                                                'TC2_current',
                                                'TC3_current',
                                                'TC1_future',
                                                'TC2_future',
                                                'TC3_future',
                                                'TC1_change',
                                                'TC2_change',
                                                'TC3_change',
                                                'TC1_ratio',
                                                'TC2_ratio',
                                                'TC3_ratio']].median()
    df_median.to_csv(data_dir.joinpath('frequency_maps', 'change_maps', f"TC_global_{resolution}as_MIT_median_ssp{ssp}_2cal_changemap.csv"), index=False)

# combine with mangrove data
for ssp in ssps:
    df_median = pd.read_csv(data_dir.joinpath('frequency_maps', 'change_maps', f"TC_global_{resolution}as_MIT_median_ssp{ssp}_2cal_changemap.csv"))
    df_median['geometry'] = gpd.GeoSeries.from_wkt(df_median['geometry'])
    median_gdf = gpd.GeoDataFrame(df_median, geometry='geometry', crs="EPSG:4326")
    median_clip = gpd.clip(median_gdf, buffered_gdf)
    # reproject to equal area (Eckert IV)
    median_clip = median_clip.to_crs(proj_eck4)
    # match points & save
    join = gpd.sjoin_nearest(exp_gdf, median_clip, distance_col="TC_distance")
    join = join.to_crs("EPSG:4326")
    join.to_csv(data_dir.joinpath('mangroves', 'frequency_maps', f'{source}_mangroves_TC_global_{resolution}as_MIT_median_ssp{ssp}_2cal_changemap.csv'), index=False)

# Analyse changemaps
## - first filter too low frequencies
## - potentially filter return period 1 in 250: f=0.004 (in future climate)
## - inf values mean newly affected, -1 values mean not affected anymore
## - store in new column and set to NaN in original columns

models = ['miroc6', 'cesm2', 'ecearth6', 'median']
# Frequency threshold (1/250years = 0.004)
threshold = 0.004

def newly_affected(ratio):
  if ratio == np.inf:
    return 1
  else:
    return 0
  
def newly_unaffected(ratio):
  if ratio == -1:
    return 1
  else:
    return 0

for ssp in ssps:
  for mod in models:
    df = pd.read_csv(data_dir.joinpath(f'mangroves', 'frequency_maps', f'{source}_mangroves_TC_global_{resolution}as_MIT_{mod}_ssp{ssp}_2cal_changemap.csv'))
    for cat in range(1, 4):
        current_col = f'TC{cat}_current'
        future_col = f'TC{cat}_future'
        change_col = f'TC{cat}_change'
        ratio_col = f'TC{cat}_ratio'
        affected_col = f'TC{cat}_new_affected' # column for newly affected
        unaffected_col = f'TC{cat}_new_unaffected' # column for no longer affected
        df[current_col] = np.where(df[current_col]<threshold, float(0), df[current_col])
        df[future_col] = np.where(df[future_col]<threshold, float(0), df[future_col]) # filter rows by frequency threshold
        df[ratio_col] = (df[future_col] - df[current_col]) / df[current_col] # recalculate ratio
        df[affected_col] = df[ratio_col].apply(newly_affected)
        df[affected_col] = np.where(df[change_col]==0, float(0), df[affected_col]) # set newly affected to 0 where there was no change
        df[unaffected_col] = df[ratio_col].apply(newly_unaffected)
        # set NaN values
        df[ratio_col] = df[ratio_col].replace(np.inf, pd.NA)
        df[ratio_col] = df[ratio_col].replace(-1, pd.NA)
        # save new dataframe
        df.to_csv(data_dir.joinpath('mangroves', 'frequency_maps', 'analysed', f'{source}_mangroves_TC_global_{resolution}as_MIT_{mod}_ssp{ssp}_2cal_changemap_analysed.csv'), index=False)
    