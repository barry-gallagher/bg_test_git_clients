# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 12:20:37 2022

@author: Casiano.Koprowski.LA
"""

import configparser as cp
from datetime import datetime, date
from glob import glob
import logging
import os
import re
import sys

from shapely.geometry import shape
from shapely.ops import unary_union
import numpy as np
from osgeo import osr, gdal, ogr
from pyproj import CRS
from pyproj.transformer import AreaOfInterest
from pyproj.exceptions import CRSError
from tqdm import tqdm
from shapely.geometry import MultiPolygon, Polygon, MultiPoint, Point, GeometryCollection, LineString, LinearRing

from fuse_dev.fuse.fuse_processor import FuseProcessor
from fuse_dev.fuse.meta_review.meta_review import flatten_geometry, MetadataDatabase, connect_with_retries, database_has_table, split_URL_port
from fuse_dev.fuse.meta_review.summary_table import MetricsDatabaseAggregate, SurveySummaryDatabase, SurveySummaryTable
from fuse_dev.fuse.raw_read.noaa.enc import write_geopackage, open_geopackage
from nbs import configs
from nbs_utils.points_utils import from_npz
from nbs_utils.gdal_utils import spatial_reference_from_metadata
from transform_dev.datum_transform.transform import reproject_horizontal_coordinates
from nbs.bruty.nbs_postgres import start_integer

LOGGER = logging.getLogger('proc.enc_sup')

PLOT_DIRECTORY = os.path.join(os.path.expanduser('~'), 'Downloads')

FILE_EXTENSIONS = ['.a93', '.bag', '.tif', '.xyz']


class Area:
    def __init__(self, o: shape):
        self.o = o['data_coverage']

    def __lt__(self, other: shape) -> bool:
        # within = self.o.within(other.o)
        size = self.o.area < other.o.area
        return size # or within


def run_configs(config_filenames: [str], ignore_filename: str, create_summary_tables: bool = True) -> [SurveySummaryTable]:
    """Assess and update data as needed."""
    start_time = datetime.now()
    log_name = f'_fuse_enc_{datetime.now():%Y%m%d_%H%M%S}.log'

    summary_tables = []
    logging_errors = {}

    for config_index, (config_filename, config_file) in enumerate(configs.iter_configs(config_path)):
        process_errors = {}
        processor = FuseProcessor(config_file, log_name)

        utm_zone = int(processor._config['to_horiz_key'])
        in_crs = CRS.from_epsg(4326)
        out_crs = CRS.from_dict({'datum': 'NAD83', 'proj': 'utm', 'zone': utm_zone, 'south': False})
        try:
            out_crs = CRS.from_epsg(out_crs.to_epsg())
        except CRSError:
            LOGGER.warning(f'No EPSG for NAD83 in utm {utm_zone}.  Proceeding with approximate definition.')
        hostname, database, table_name = processor._config['enc_table'].split('/')

        with open(processor._config['credentials_file']) as database_credentials_file:
            username, password = [line.strip() for line in database_credentials_file][:2]
        direct_hostname, port = split_URL_port(hostname)
        processor._meta_obj = MetadataDatabase(hostname, database, username, password, table_name, fields=processor._cols, crs=out_crs)
        processor._meta_obj.crs = out_crs
        LOGGER.info(f"Working in {database}.{table_name}.")

        try:
            LOGGER.info("Transforming geometry associated with ENCs")
            # Transforming and storing the original ENC MCOVR geometries
            records = processor.records_for_read_type()
            # def only_process_searched(records_list):
            #     # this is for debugging bad names/sorind records
            #     return
            #     for i in range(len(records_list)-1, -1, -1):
            #         record = records_list[i]
            #         if "US4DE11M" not in record['survey_id']:
            #             records_list.pop(i)
            # only_process_searched(records)
            for record in tqdm(records, mininterval=.7):
                tr_poly = []
                geom = record['provided_coverage']
                if isinstance(geom, MultiPolygon):
                    provided_coverage = record['provided_coverage']
                    for geom in provided_coverage.geoms:
                        geom_ring = list(zip(*geom.exterior.xy))
                        try:
                            tr_points = reproject_horizontal_coordinates(geom_ring, in_crs, out_crs, strict_xy_axis_order=True)
                        except ValueError as error:
                            LOGGER.warning(f"{error.__class__.__name__}: {error} --> attempting to pass in an Area of Interest to resolve conflict.  No further action needed.")
                            geom_ring = np.array(geom_ring)
                            aoi_bounds = AreaOfInterest(np.max(geom_ring[:, 0]), np.min(geom_ring[:, 1]),
                                                        np.min(geom_ring[:, 0]), np.max(geom_ring[:, 1]))
                            tr_points = reproject_horizontal_coordinates(geom_ring, in_crs, out_crs,
                                                                         strict_xy_axis_order=True,
                                                                         area_of_interest=aoi_bounds)

                        tr_poly.append(Polygon(tr_points))
                else:
                    geom_ring = list(zip(*geom.exterior.xy))
                    try:
                        tr_points = reproject_horizontal_coordinates(geom_ring, in_crs, out_crs, strict_xy_axis_order=True)
                    except ValueError as error:
                        LOGGER.warning(f"{error.__class__.__name__}: {error} --> attempting to pass in an Area of Interest to resolve conflict.  No further action needed.")
                        geom_ring = np.array(geom_ring)
                        aoi_bounds = AreaOfInterest(np.max(geom_ring[:, 0]), np.min(geom_ring[:, 1]),
                                                    np.min(geom_ring[:, 0]), np.max(geom_ring[:, 1]))
                        tr_points = reproject_horizontal_coordinates(geom_ring, in_crs, out_crs,
                                                                     strict_xy_axis_order=True,
                                                                     area_of_interest=aoi_bounds)
                    tr_poly.append(Polygon(tr_points))
                record['data_coverage'] = flatten_geometry(MultiPolygon(tr_poly))
                processor._meta_obj[record['from_filename']] = record

            for ridx, record in enumerate(records):
                LOGGER.info(f"Working on {record['from_filename']}: {ridx+1} of {len(records)}")
                if 'to_filename' not in record or record['to_filename'] == '':
                    msg = "No to_filename found.  Marked as never_post."
                    LOGGER.warning(msg)
                    if 'notes' in record and msg not in record['notes']:
                        record['notes'] = f"{record['notes']}; {msg}"
                    else:
                        record['notes'] = msg
                    record['never_post'] = True
                    record['interpolation_mask'] = None
                    processor._meta_obj[record['from_filename']] = record
                    continue
                else:
                    subdir_path = os.path.split(record['from_path'])[0]
                    original_path = os.path.split(subdir_path)[0]
                    file_marker = os.path.join(original_path, 'THIS CHART IS RETIRED')
                    depricated = False
                    if not os.path.exists(original_path) or os.path.exists(file_marker):
                        depricated = True

                    if depricated:
                        msg = 'Tile marked as depricated. Marked as never_post'
                        LOGGER.warning(msg)
                        if 'notes' in record and msg not in record['notes']:
                            record['notes'] = f"{record['notes']}; {msg}"
                        else:
                            record['notes'] = msg
                        record['never_post'] = True
                        record['interpolation_mask'] = None
                        processor._meta_obj[record['from_filename']] = record
                        continue
            # This is here to make sure that surveys with never_post = True flags are not included in the list.
            records_no_retired = [record for record in records if not record['never_post'] == True]
            # Sort polygons by largest to smallest (area)
            sorted_areas = sorted(records_no_retired, key=Area, reverse=True)
            modified_records = []

            for ridx, record in enumerate(sorted_areas):
                # Ensure the geometry is valid
                geom = record['data_coverage']
                if not geom.is_valid:
                    print(f"Invalid geometry for {record['from_filename']}")

                # Loop over smaller surveys
                for item in reversed(sorted_areas[ridx + 1:]):
                    if geom.intersects(item['data_coverage']):
                        other_data_coverage = item['data_coverage']

                        # Check and clean geometries
                        if not geom.is_valid:
                            geom = geom.buffer(0)
                        if not other_data_coverage.is_valid:
                            other_data_coverage = other_data_coverage.buffer(0)
                        geom = filter_geometries_with_area(geom)
                        if geom is None:  # only geometries were bad slivers (line segments) from intersections
                            break
                        other_data_coverage = filter_geometries_with_area(other_data_coverage)

                        if record['point_spacing'] > item['point_spacing']:
                            geom = geom.difference(other_data_coverage)
                        elif record['point_spacing'] == item['point_spacing']:
                            intersection = geom.intersection(other_data_coverage)
                            if intersection.area > 0:
                                LOGGER.warning(f"Same chart scale supersession was performed for overlapping tiles {record['from_filename']} and {item['from_filename']}")
                                geom = geom.difference(other_data_coverage)
                        elif record['point_spacing'] < item['point_spacing']:
                            item['data_coverage'] = other_data_coverage.difference(geom)
                            LOGGER.warning(f"Removing {item['from_filename']} overlap from {record['from_filename']} due to point spacing.")
                        else:
                            raise NotImplementedError(f'Found a chart comparison for {record["from_filename"]} that does not have a point spacing that matches the current comparison cases.')
                    else:
                        # if a geometry does not intersect another geometry there isn't anything to do here.
                        pass

                # Update interpolation_mask and append the modified record
                record['interpolation_mask'] = filter_geometries_with_area(geom)
                modified_records.append(record)
            # We will use modified_records for further processing if needed

            # Use intersection to trim points from outside of the resultant areas
            for record in modified_records:
                if record['interpolation_mask'] is not None:
                    cookie_cutter = record['interpolation_mask']
                    output = record['to_filename']
                    output_folder = os.path.split(output)[0]
                    original_output = os.path.join(output_folder, f"{record['from_filename']}.gpkg")
                    data = open_geopackage(original_output)

                    # Original ENC Z values are swapped for an index
                    pindex = np.array(list(range(len(data))))
                    points = data.copy()
                    points[:, 2] = pindex
                    points = points[:, [0, 1, 2]].astype(float)

                    # Trim points via intersection between points and geometry
                    mp_shape = MultiPoint(points)
                    trimmed_shape = mp_shape.intersection(cookie_cutter)
                    if not trimmed_shape.is_empty and trimmed_shape.is_valid:
                        if isinstance(trimmed_shape, Point):
                            trimmed_shape = MultiPoint([trimmed_shape])

                        # Use indexes of remaining points to pull from the original data
                        trimmed_idx = np.array([(*point.xy[0], *point.xy[1], point.z) for point in trimmed_shape.geoms])
                        trimmed_points = data[trimmed_idx[:, 2].astype(int)]

                        # Write out new files and update 'to_filename'
                        npy_file = os.path.join(output_folder, f"{record['from_filename']}_trimmed.npy")
                        np.save(npy_file, trimmed_points)

                        gpkg_file = os.path.join(output_folder, f"{record['from_filename']}_trimmed.gpkg")
                        record['read_type'] = 'enc'
                        processor._point_writer.write(npy_file, gpkg_file, record)
                        record['to_filename'] = gpkg_file
                        processor._meta_obj[record['from_filename']] = record
                        if os.path.exists(npy_file):
                            os.remove(npy_file)
                    else:
                        msg = f"No proper ENC soundings remain after supersession for {record['from_filename']}. Keeping the interpolation mask data in enc table to potentially represent land area."
                        LOGGER.warning(msg)
                        if 'notes' in record and msg not in record['notes']:
                            record['notes'] = f"{record['notes']}; {msg}"
                        else:
                            record['notes'] = msg
                            
                        record['never_post'] = True
                        processor._meta_obj[record['from_filename']] = record
                else:
                    msg = "ENC superseded entirely. Original output kept as to_filename and marked as never_post."
                    LOGGER.warning(msg)
                    if 'notes' in record and msg not in record['notes']:
                        record['notes'] = f"{record['notes']}; {msg}"
                    else:
                        record['notes'] = msg

                    record['never_post'] = True
                    record['interpolation_mask'] = None
                    processor._meta_obj[record['from_filename']] = record

            records = processor.records_for_read_type()
            # only_process_searched(records)
            # Create new table for rehashed data
            hostname, database, table_name_ = processor._config['enc_table'].split('/')
            table_name = f"{table_name_}_{date.today():%Y%m%d}"

            with open(processor._config['credentials_file']) as database_credentials_file:
                username, password = [line.strip() for line in database_credentials_file][:2]
            # FIXME there is a smallgotcha here, if you run this one the same day to fix a sorind that needs correction,
            #   the previous records will still be in the MetadataDatabase we are about to connect to.  
            #   Basically it assumes you aren't running supersession twice in one day.
            new_meta = MetadataDatabase(hostname, database, username, password, table_name, fields=processor._cols, crs=out_crs)
            LOGGER.info(f"Working in {database}.{table_name}.")
            remapped_in_region = {}

            # Rehash ENC data and store by SORIND values
            for ridx, record in enumerate(records):
                try:
                    LOGGER.info(f"Splitting {record['from_filename']}: {ridx+1} of {len(records)}")
                    output = record['to_filename']
                    enc_name = record['survey_id']
                    enc_date = record['end_date']
                    output_folder = os.path.split(output)[0]
                    trimmed_output = os.path.join(output_folder, f"{record['from_filename']}_trimmed.gpkg")
                    data = open_geopackage(trimmed_output)

                    # Replacing None SORIND and SORDAT with 'Blank'
                    if any(data[:, 4] == None):
                        new_sorind = "Blank"
                        LOGGER.warning(f"Encountered soundings with no SORIND, given: {new_sorind}")
                        data[data[:, 4] == None, 4] = f'{new_sorind}'
                    if any(data[:, 5] == None):
                        new_sordat = "Blank"
                        LOGGER.warning(f"Encountered soundings with no SORDAT, given: {new_sordat}")
                        data[data[:, 5] == None, 5] = new_sordat

                    # Correcting bad SORIND
                    LOGGER.info("Resolving SORIND names and correcting as needed.")
                    sorind_list = np.unique(data[:, 4])
                    ident_breakdown = {'g': 'graph', 'r': 'reprt'}
                    format_string = f"US,U[1-4S],(?:{('|').join(ident_breakdown.values())}),"
                    ident_chars = f"[{('').join(ident_breakdown.keys())}][a-z]" + '{4}'
                    naming_standard = re.compile(format_string)  #, re.IGNORECASE)
                    index_breakdown = {
                        0: (2, re.compile("US", re.IGNORECASE)),
                        3: (2, re.compile("U[1-4S]", re.IGNORECASE)),
                        6: (5, re.compile(ident_chars, re.IGNORECASE)),
                        }
                    remapped_sorinds = {}
                    for sorind in sorind_list:
                        if sorind == "Blank":
                            continue

                        matched = True
                        matches_standard = re.match(naming_standard, sorind)
                        subset_sorind, source_name = sorind[:12], sorind[12:]
                        corrected = ''

                        if matches_standard is None:
                            matched = False
                            for stind, descr in index_breakdown.items():
                                enind, frmt = descr
                                section = subset_sorind[stind:stind+enind]
                                section_match = re.match(frmt, section)
                                if section_match is None:
                                    if stind == 6 and section.lower() not in ident_breakdown.values():
                                        # 'nsurf' and 'survy' should be 'graph'
                                        corrected += f"{ident_breakdown['g']},"
                                elif stind in (0, 3):
                                    section = section.upper()
                                    corrected += f"{section},"
                                elif stind == 6:
                                    letter = section[0].lower()
                                    corrected += f"{ident_breakdown[letter]},"
                            if corrected != 'Blank':
                                corrected += source_name
                        # There is a duplicate sorind in US4DE11M for Chart 12214
                        if sorind == 'US,U1,graph,Chart 12214':
                            corrected = 'US,US,graph,Chart 12214'
                            matched = False
                        # duplicate "graph" or "reprt" in the name -- found a graph,graph and graph,reprt
                        # If not handled Bruty will report "Two IDs reference the same file ID"
                        if re.search(r",\s*(graph|reprt),\s*(graph|reprt)", sorind, re.IGNORECASE):
                            if corrected:
                                pass  # I think this corrects the leading part of the sorind
                            else:
                                # remove all ",graph" or ",reprt" that are followed another ",graph" or ",reprt"
                                # FIXME We have seen cases where the second of the two graph,reprt is correct so remove the first
                                #   if the first occurence is correct we need better logic to find the matches
                                corrected = re.sub(r"(,\s*(graph|reprt))(?=,\s*(graph|reprt))", "", sorind, flags=re.IGNORECASE)
                                matched = False
                        # double comma where only one should be -- US5ILMEC / DD-38023 had this
                        # FIXME if sorind has two errors then only one will get corrected currently 
                        #   (i.e. corrected is not checked but is replaced)
                        if re.search(r"\s*(graph|reprt),\s*,", sorind, re.IGNORECASE):
                            if corrected:
                                pass  # I think this corrects the leading part of the sorind
                            else:
                                # remove all ",graph" or ",reprt" that are followed another ",graph" or ",reprt"
                                # FIXME We have seen cases where the second of the two graph,reprt is correct so remove the first
                                #   if the first occurence is correct we need better logic to find the matches
                                corrected = re.sub(r",\s*,", ",", sorind, flags=re.IGNORECASE)
                                matched = False
                        split_sorind = sorind.split(',')
                        if bool(re.search(r'^\s+|\s+$', split_sorind[-1])):
                            matched = False
                            if corrected != '':
                                split_sorind = corrected.split(',')
                                corrected = (',').join([section.strip() for section in split_sorind])
                            else:
                                corrected = (',').join([section.strip() for section in split_sorind])

                        if 'chart' in source_name:
                            matched = False
                            if corrected != '':
                                corrected = corrected.replace('chart', 'Chart')
                            else:
                                corrected = sorind.replace('chart', 'Chart')

                        if not matched:
                            if corrected == 'Blank':
                                LOGGER.warning(f"SORIND '{sorind}' could not be resolved to naming standards, replacing with: '{corrected}'")
                            else:
                                LOGGER.warning(f"SORIND '{sorind}' does not match naming standards, replacing with: '{corrected}'")
                            if corrected not in remapped_sorinds:
                                remapped_sorinds[corrected] = [sorind]
                            else:
                                remapped_sorinds[corrected].append(sorind)
                            data[data[:, 4] == sorind, 4] = corrected

                    if len(list(remapped_sorinds.keys())) != 0:
                        if enc_name not in remapped_in_region:
                            remapped_in_region[enc_name] = [remapped_sorinds]
                        else:
                            remapped_in_region[enc_name].append(remapped_sorinds)

                    LOGGER.info("Resolving SORIND/SORDAT pairings for output; SORDAT is corrected as needed.")
                    sorind_list = np.unique(data[:, 4])
                    for sorind in sorind_list:

                        data_idx = np.where(data[:, 4] == sorind)
                        sorind_subsets = data[data_idx]
                        sordat_list = np.unique(sorind_subsets[:, 5])
                        for sordat in sordat_list:
                            LOGGER.info(f"Resolving new entry for primary SORIND {sorind} and secondary SORDAT {sordat}")
                            data_idx = np.where(sorind_subsets[:, 5] == sordat)
                            constituants = sorind_subsets[data_idx]
                            subset_meta = record.copy()

                            # Accounting for bad SORDAT
                            try:
                                date_format = "%Y%m%d"
                                if sordat == 'Blank':
                                    raise ValueError("No SORDAT provided; Blank is used as SORDAT")
                                else:
                                    augmented = False
                                    if len(sordat) >= 4:
                                        y = int(sordat[:4])
                                        yr = sordat[:4]
                                        if len(sordat) >= 6:
                                            m = int(sordat[4:6])
                                            mr = sordat[4:6]
                                            if len(sordat) >= 8:
                                                d = int(sordat[6:])
                                                dr = sordat[6:]
                                            else:
                                                d = 1
                                                dr = '01'
                                                augmented = True
                                        else:
                                            m = 1
                                            mr = '01'
                                            d = 1
                                            dr = '01'
                                            augmented = True

                                    if not (1 < m < 12):
                                        augmented = True
                                        mr = '01'

                                    if m == 2:
                                        mmax = 29
                                    elif m in (1, 3, 5, 7, 8, 10, 12):
                                        mmax = 31
                                    else:
                                        mmax = 30
                                    if not (1 < d < mmax):
                                        augmented = True
                                        dr = '01'

                                    exp_sordat = yr + mr + dr
                                    if augmented:
                                        LOGGER.info(f"SORDAT modified for 'end_date' and 'start_date'. -> New: {exp_sordat}; Old: {sordat}")

                                    subset_meta['end_date'] = datetime.strptime(exp_sordat, date_format).date()
                                    subset_meta['start_date'] = subset_meta['end_date']
                            except ValueError as error:
                                LOGGER.warning(f"Unable to set 'end_date' for SORIND/SORDAT pairing {sorind}/{sordat}; 'end_date' is set to None. This error will not break processing -> {error.__class__.__name__}: {error}")
                                subset_meta['end_date'] = None
                                subset_meta['start_date'] = None

                            # Filling out metadata
                            subset_meta['from_filename'] = f"{sorind}_{sordat}_{enc_name}_{enc_date:%Y%m%d}"

                            if enc_name.lower() not in sorind.lower():
                                if 'chart' not in sorind.lower():
                                    subset_meta['source_indicator'] = f"{sorind}_{enc_name}"
                                elif sorind == 'Blank':
                                    subset_meta['source_indicator'] = sorind
                                else:
                                    subset_meta['source_indicator'] = f"{sorind}"

                            if sorind != 'Blank':
                                source_type = subset_meta['from_filename'].split(',')[-2].strip()
                                if source_type == 'reprt':
                                    first_name = f"{sorind}-{source_type}_{sordat}_{enc_name}_{enc_date:%Y%m%d}"
                                    sub_name = first_name.split(',')[-1].strip()
                                    subset_meta['survey_id'] = subset_meta['from_filename'].split(',')[-1].strip()
                                else:
                                    sub_name = subset_meta['from_filename'].split(',')[-1].strip()
                                    subset_meta['survey_id'] = sub_name
                            else:
                                sub_name = subset_meta['from_filename'].split(',')[-1].strip()
                                subset_meta['survey_id'] = sub_name

                            subset_meta['from_path'] = record['to_filename']
                            subset_meta['provided_coverage'] = subset_meta['interpolation_mask']
                            subset_meta['catzoc'] = 5
                            subset_meta['decay_score'] = 1
                            subset_meta['supersession_score'] = 20
                            subset_meta['data_coverage'] = None
                            subset_meta['interpolation_mask'] = None
                            subset_meta['to_filename'] = None

                            if sorind in remapped_sorinds:
                                subset_meta['notes'] = f"Contains corrected SORIND values from {', '.join(remapped_sorinds[sorind])}"
                            new_meta[subset_meta['from_filename']] = subset_meta.flatten()
                            LOGGER.info(f"SORIND/SORDAT pair resolved as new 'from_filename': {subset_meta['from_filename']}; saving output")

                            try:
                                constituants = constituants[np.lexsort((constituants[:, 1], constituants[:, 0]))]
                                is_same = False
                                
                                # cleaning up the survey_id and source_indicator if needed - see Issue: TASK #127494
                                sub_name_sanitized = sanitize_filename(sub_name, option = 'geopackage_file_name')
                                if sub_name_sanitized != sub_name: 
                                    original_survey_id = subset_meta['survey_id'] 
                                    original_source_indicator = subset_meta['source_indicator']
                                    subset_meta['survey_id'] = sanitize_filename(original_survey_id, option = 'others') 
                                    subset_meta['source_indicator'] = sanitize_filename(original_source_indicator, option = 'others')

                                    LOGGER.info(f"Survey ID changed from {original_survey_id} to {subset_meta['survey_id']} and Source Indicator changed from {original_source_indicator} to {subset_meta['source_indicator']}")
                                # Check if the geopackage of the constituent on disk holds the same data, if so then don't write again.
                                # This speeds up the combine process later in the NBS workflow
                                new_gpkg = os.path.join(output_folder, f"{sub_name_sanitized}.gpkg")
                                # This would be more convenient and faster with geopandas
                                if os.path.exists(new_gpkg):
                                    try:
                                        gpkg = gdal.OpenEx(new_gpkg)
                                    except:
                                        pass  # failed to open, just write a new file
                                    else:
                                        if gpkg and gpkg.GetLayerCount() == 1:
                                            lyr = gpkg.GetLayer(0)  # There should only be one layer in these files
                                            if lyr.GetGeomType() == ogr.wkbPoint:
                                                comp_old_srs = lyr.GetSpatialRef().ExportToWkt()
                                                # The geopackage writer ir changing the CRS from out_crs using the spatial_reference_from_metadata
                                                # which loses the authority code and changes the name as well, but we have to match it anyway
                                                comp_new_crs = spatial_reference_from_metadata(subset_meta).ExportToWkt()
                                                gpkg_data = open_geopackage(new_gpkg)
                                                # Sort how constituants was sorted
                                                gpkg_data = gpkg_data[np.lexsort((gpkg_data[:, 1], gpkg_data[:, 0]))]
                                                # The open_geopackage() reads everything as strings, so have to switch to numbers for numpy comparison
                                                try:
                                                    is_same_pts = np.all(np.equal(constituants[:, [0, 1, 2, 3]].astype(np.float64),
                                                                                  gpkg_data[:, [0, 1, 2, 3]].astype(np.float64)))
                                                except ValueError:  # happens if point count doesn't match
                                                    is_same = False
                                                else:
                                                    is_same_srs = comp_old_srs == comp_new_crs
                                                    is_same = is_same_pts and is_same_srs
                                            del lyr
                                    finally:
                                        if gpkg:
                                            del gpkg
                                subset_meta['read_type'] = 'enc'
                                subset_meta['to_filename'] = new_gpkg
                                new_meta[subset_meta['from_filename']] = subset_meta
                                if not is_same:
                                    # Saves the data for gpkg, npz
                                    msg = f"Supersession results for {subset_meta['from_filename']} have new/updated results. New files will be saved."
                                    # Write out new files and update 'to_filename'
                                    new_npy = os.path.join(output_folder, f"{sub_name_sanitized}.npy")
                                    np.save(new_npy, constituants)

                                    processor._point_writer.write(new_npy, new_gpkg, subset_meta)
                                    if os.path.exists(new_npy):
                                        os.remove(new_npy)
                                else:
                                    msg = f"Supersession results for {subset_meta['from_filename']} matched current results. Did not replace current products."
                                LOGGER.info(msg)
                            except Exception as error:
                                error_text = f"{subset_meta['from_filename']} encountered issue while writing out data; set to never_post --> {error.__class__.__name__}: {error}"
                                subset_meta['never_post'] = True
                                new_meta[subset_meta['from_filename']] = subset_meta
                                LOGGER.exception(error_text)
                                if config_filename not in process_errors:
                                    process_errors[config_filename] = [error_text]
                                else:
                                    process_errors[config_filename].append(error_text)
                                    del error_text

                except RuntimeError as error:
                    error_text = f"RuntimeError: Error reading ENC geopackage for {record['from_filename']}"
                    LOGGER.exception(error_text)

            formatted = pretty(remapped_in_region)
            LOGGER.info(f"Bad SORIND assignents corrected:{formatted}")

            # Move current finalized records to final production table
            records = new_meta.records_where({'from_filename': '!None'}, never_post_check=True)
            # Stores the new records to be written into the truncated metadata table
            flat_pack = [record.flatten() for record in records]
            # Create new table in metadata databse for rehashed data
            hostname, database, table_name = processor._config['metatable'].split('/')

            with open(processor._config['credentials_file']) as database_credentials_file:
                username, password = [line.strip() for line in database_credentials_file][:2]
            direct_hostname, port = split_URL_port(hostname)
            LOGGER.info(f"Working in {database}.{table_name}.")
            meta_table = MetadataDatabase(hostname, database, username, password, table_name, fields=processor._cols, crs=out_crs)
            meta_table.crs = out_crs

            # Check if the corresponding meta table exists
            check_table_query = f"SELECT 1 FROM information_schema.tables WHERE table_schema = 'public' AND table_name = '{table_name}' LIMIT 1;"
            cur = meta_table.connection.cursor()
            cur.execute(check_table_query)
            table_exists = cur.fetchone() is not None
            cur.close()

            # Truncate (remove all the entries) in the table if it exists - cleaning up the table
            # NOTE: This script only deletes entries in the metadata table and not in the enc dated table - {table_name}_{date}.
            #       As the metadata table is populated using entries from the enc dated table, you should not
            #       expect metadata table entries to be removed for potentail entry removal that occur in the corresponding enc table
            #       if the changes occur on the same day that the original enc table was manipulated.
            #       If you want the removal of the entries from metadata table to happen on the same day,
            #       you need to delete the dated enc table manually. Otherwise, wait for one day and then run supersession!
            start_val_for_nbs_id = start_integer(table_name)
            if table_exists:
                # store the old records to match nbs_id values
                records_with_nbs_id = {record['from_filename']: record for record in meta_table.records}
                max_nbs_id = max([record['nbs_id'] for record in records_with_nbs_id.values()]) + 1

                LOGGER.info(f"Truncating entries from the metadata table {table_name}.")
                truncate_query = f"TRUNCATE TABLE {table_name};"
                try:
                    cur = meta_table.connection.cursor()
                    cur.execute(truncate_query)
                    meta_table.connection.commit()
                    cur.close()
                except Exception as e:
                    LOGGER.error(f"Error occurred while  truncating table {table_name}: {e}")
            else:
                LOGGER.warning(f"The table {table_name} does not exist. This is a problem!")
                records_with_nbs_id = {}
                max_nbs_id = start_val_for_nbs_id

            # Reset the numbering for nbs_id column if the column already exists,
            # if the column does not exists, create it and set the correct start value
            handle_nbs_id_column(meta_table, table_name, start_val_for_nbs_id)

            LOGGER.info(f"{len(flat_pack)} records confirmed for final table")
            # Writes the stored records into the metadata table
            meta_table.insert_records(flat_pack)
            # Look at all the records in the new metadata table and use the old nbs_id if it exists.
            # Otherwise use the max_nbs_id we found and increase it for future use.
            # We have to use direct SQL as there is protection against setting nbs_id in the fuse meta_review .insert_records() code.
            cur = meta_table.connection.cursor()
            nbs_id_name = "nbs_id"
            prim_key_name = 'from_filename'
            for record in flat_pack:
                prim_key = record[prim_key_name]
                if prim_key in records_with_nbs_id:
                    record[nbs_id_name] = records_with_nbs_id[prim_key][nbs_id_name]
                else:
                    record[nbs_id_name] = max_nbs_id
                    max_nbs_id += 1
                # TODO make this an executemany
                cur.execute(f"""UPDATE {table_name} SET {nbs_id_name} = {record[nbs_id_name]} WHERE {prim_key_name} = '{prim_key}'""")
            cur.close()
            # We just updated the nbs_id values, so we need to reset the counter for the next record that gets added automatically
            handle_nbs_id_column(meta_table, table_name, max_nbs_id)

        except Exception as error:
            error_text = f"update supersession error for {config_filename} - {error.__class__.__name__}: {error}"
            LOGGER.exception(error_text)
            if config_filename not in process_errors:
                process_errors[config_filename] = [error_text]
            else:
                process_errors[config_filename].append(error_text)
                del error_text

        try:
            if create_summary_tables:
                summary_table = update_summary_table(processor)
                summary_tables.append(summary_table)
        except Exception as error:
            error_text = f"Did not create summary table for {config_filename} - {error.__class__.__name__}: {error}"
            LOGGER.exception(error_text)
            if 'summary_tables' not in process_errors:
                process_errors['summary_tables'] = [error_text]
            else:
                process_errors['summary_tables'].append(error_text)
            del error_text

        if len(process_errors) != 0:
            formatted_errors = []
            for survey_id, errors in process_errors.items():
                formatted_errors.append(f"{config_filename}")
                for text in errors:
                    formatted_errors.append(f"\t{text}")
            formatted_text = '\n'.join(formatted_errors)
            LOGGER.warning(f"These errors were reported during this run:\n{formatted_text}")
            logging_errors.update(process_errors)

    end_time = datetime.now()
    LOGGER.info(f'Completed supersession at {end_time}, took {end_time - start_time}')
    if len(logging_errors) != 0:
        formatted_errors = []
        for config_filename, errors in logging_errors.items():
            formatted_errors.append(f"{config_filename}")
            for text in errors:
                formatted_errors.append(f"\t{text}")
        formatted_text = '\n'.join(formatted_errors)
        print(f"These errors were reported during this run:\n{formatted_text}")
    processor.close_process_log()

    return summary_tables


def pretty(d: dict, indent: int = 0) -> str:
    """
    https://stackoverflow.com/a/3229493.

    Parameters
    ----------
    d : dict
        dictionary.
    indent : int, optional
        indent of nested values. The default is 0.

    Returns
    -------
    str
        string representation of dict.

    """
    text = ''
    for key, value in d.items():
        text += '\n' + ('\t' * indent) + str(key)
        if isinstance(value, dict):
            text += pretty(value, indent+1)
        elif isinstance(value, (tuple, list)):
            for item in value:
                text += '\n' + ('\t' * (indent+1)) + str(item)
        else:
            text += '\n' + ('\t' * (indent+1) + str(value))
    return text


def load_ignored(ignore_filename: str):
    """
    Load a list of strings that represent surveys to be ignored upon read.
    """

    ignore_list = []
    if os.path.isfile(ignore_filename):
        with open(ignore_filename) as ignore_file:
            ignore_list.extend(line.strip() for line in ignore_file if line != '')
    return ignore_list


def new_from_record_survey_ids(records: [dict], survey_directories: [str]) -> {str: str}:
    """
    Indentify new surveys from existing metadata records.

    Parameters
    ----------
    records : [dict]
        List of metadata records.
    survey_directories : [str]
        Input survey folder paths.

    Returns
    -------
    {str: str}
        A dict of 'survey_id': 'survey_directory' pairs of new surveys

    """
    surveys = []
    new_survey_folders = {}

    for record in records:
        if 'survey_id' in record:
            if record['survey_id'] not in surveys:
                surveys.append(record['survey_id'])

    for survey_folder in survey_directories:
        survey_id = os.path.basename(survey_folder)
        if survey_id not in surveys:
            new_survey_folders[survey_id] = survey_folder

    return new_survey_folders


def resort_directory_list(survey_directories: [str]) -> [str]:
    """
    Moves surveys with VR BAGs to the bottom of the list.

    Parameters
    ----------
    survey_directories : [str]
        Input survey folder paths.

    Returns
    -------
    [str]
        Sorted survey folder paths.

    """
    holding = []
    survey_holding = survey_directories.copy()

    for survey in survey_directories:
        survey_files = os.listdir(survey)
        if any(['VR' not in survey_file for survey_file in survey_files]):
            survey_holding.remove(survey)
            holding.append(survey)

    survey_holding.extend(holding)
    return survey_holding


def vr_only(survey_directories: [str]) -> [str]:
    """
    Moves surveys with VR BAGs to the bottom of the list

    Parameters
    ----------
    survey_directories : [str]
        Input survey folder paths.

    Returns
    -------
    [str]
        Surveys considered to contain VR BAGs.

    """
    survey_holding = survey_directories.copy()

    for survey in survey_directories:
        survey_files = os.listdir(survey)
        has_vr = False
        for item in survey_files:
            if '_VR_' in item.upper():
                has_vr = True
        if not has_vr:
            survey_holding.remove(survey)

    return survey_holding


def update_summary_table(processor: FuseProcessor) -> SurveySummaryTable:
    config = processor._config
    region = str(os.path.abspath(config['rawpath']).split('\\')[2])
    source_type = config['raw_reader_type']

    metadata_object = processor._meta_obj

    username = metadata_object.username
    password = metadata_object.password
    hostname, database, table = config['metatable'].split('/')

    summary_table = SurveySummaryDatabase(hostname, username, password, region=region, source_type=source_type, metadata_database=database, metadata_table=table)
    summary_table.update()
    return summary_table


def plot(regions: [str], source_types: [str], hostname: str, username: str, password: str, metadata_database: str):
    now = datetime.now()
    golden_ratio = (1 + 5 ** 0.5) / 2
    figure_width_inches = 16

    run_date = datetime(2020, 1, 17)

    metrics_aggregate = MetricsDatabaseAggregate(hostname, username, password, 'summary', run_date=run_date, regions=regions, source_types=source_types, metadata_database=metadata_database)

    pyplot.figure(figsize=(figure_width_inches, figure_width_inches / golden_ratio))
    metrics_aggregate.pie_chart(show=True)
    pyplot.savefig(os.path.join(PLOT_DIRECTORY, f'metrics_pie_{"_".join(regions)}_{now:%Y%m%d%H%M%S}.png'), bbox_inches='tight')

    pyplot.figure(figsize=(figure_width_inches, figure_width_inches / golden_ratio))
    metrics_aggregate.bar_graph(show=True)
    pyplot.savefig(os.path.join(PLOT_DIRECTORY, f'metrics_bar_{"_".join(regions)}_{now:%Y%m%d%H%M%S}.png'), bbox_inches='tight')


def create_summary(config_filename: str, regions: dict, source_types: dict, database: str = 'metadata_develop'):
    """
    Create and save summary tables for current run.

    Parameters
    ----------
    config_filename : str
        DESCRIPTION.
    regions : dict
        DESCRIPTION.
    source_types : dict
        DESCRIPTION.
    database : str, optional
        DESCRIPTION. The default is 'metadata_develop'.

    Returns
    -------
    None.

    """
    config = cp.ConfigParser(config_filename)
    hostname = config['Default']['metatable'].split('/')[0]
    with open(config['Default']['credentials']) as database_credentials_file:
        lines = [line.strip() for line in database_credentials_file.readlines()]
        username, password = lines[:2]

    plot(regions, source_types, hostname, username, password, database)


def filter_geometries_with_area(geometries, logger_output = False):
    """
    Filter geometries.

    Parameters:
    - geometries: Shapely geometry or list of geometries
    - logger_output: Whether or not the output of the function is recorded. default is False

    Returns:
    - Filtered geometries as Polygon or MultiPolygon
    """

    def filter_single_geometry(geom):
        if isinstance(geom, (Point, LineString, Polygon, MultiPolygon, LinearRing)):
            if isinstance(geom, (Polygon, MultiPolygon)) and geom.area >= 1:
                return flatten_geometry(geom)

            elif isinstance(geom, (LineString, LinearRing)) and geom.is_closed:
                polygon_geom = Polygon(geom)
                if polygon_geom.area >= 1:  # FIXME - this will fail if coordinates are degrees
                    return flatten_geometry(polygon_geom)

        if logger_output == True:
            LOGGER.warning(f"Geometry {type(geom)} is not valid or has area < 1. Number of geometries before: {num_before} and after: {num_after}")
        return None

    num_before = 0
    num_after = 0
    
    def count_geometries(geom):
        nonlocal num_before
        nonlocal num_after
        if isinstance(geom, (Point, LineString, LinearRing, Polygon)):
            num_before += 1
            updated_geom = filter_single_geometry(geom)
            if updated_geom:
                num_after += 1
            return updated_geom
        elif isinstance(geom, (MultiPolygon, GeometryCollection)):
            num_before += len(geom.geoms)
            filtered_geometries = [filter_single_geometry(sub_geom) for sub_geom in geom.geoms if filter_single_geometry(sub_geom) is not None]
            if filtered_geometries:
                num_after += len(filtered_geometries)
                return unary_union(filtered_geometries)  # Combine into a MultiPolygon
        return None
    
    result = count_geometries(geometries)
    
    return result


def sanitize_filename(input_name, option = 'others'):

    """
    Sanitize a given filename or file path by replacing invalid characters.

    Parameters:
    - input_name (str): The original filename or file path to be sanitized.
    - option (str): The sanitization option. Default is 'others'.
                   Options:
                   - 'geopackage_file_name': Sanitize for GeoPackage file names, replacing specific characters
                     with unique patterns.
                   - 'others': Replace characters not in [a-zA-Z0-9-,_ \n\.] with '-'.

    Returns:
    str: The sanitized filename or file path based on the specified option.
    """

    if option == 'geopackage_file_name':
        # Define a regular expression pattern to match characters not allowed in file names
        invalid_chars_pattern = re.compile(r'[<>:"=/\\|?*\n\r\t]')

        # Create a mapping of non-allowed characters to unique replacements
        replacement_mapping = {'<': 'Lt', '>': 'Gt', ':': 'Colon', '"': 'Quote', '/':'-', '\\': 'Backslash',
                            '|': 'Pipe', '?': 'Question', '*': 'Asterisk', '\n': 'Newline', '\r': 'Return', '\t': 'Tab',
                            '=':'Equal'}

        # Check if the input name contains invalid characters
        if invalid_chars_pattern.search(input_name):
            # Replace invalid characters with a pattern indicating the original character and its replacement
            def replace_callback(match):
                char = match.group()
                replacement = replacement_mapping.get(char, char)
                return f"{replacement}"

            # Generate the updated name by applying replacements and adding underscores
            sanitized_name = invalid_chars_pattern.sub(replace_callback, input_name)
            
            LOGGER.info(f"Geopackage file name changed from {input_name} to {sanitized_name}")

            return sanitized_name
        else:
            # If no invalid characters are found, return the original name
            return input_name

    elif option == 'others':
        return re.sub('[^a-zA-Z0-9-,_ \n\.]', '-', input_name)

def get_sequence_name(cur, table_name):
    cur.execute(f"SELECT relname FROM pg_class WHERE relname LIKE '{table_name}_nbs_id_seq%';")
    sequence_row = cur.fetchone()
    sequence_name = sequence_row[0] if sequence_row else None
    return sequence_name

def get_current_start_value(cur, sequence_name):
    cur.execute(f"SELECT last_value FROM {sequence_name};")
    current_start_val = cur.fetchone()[0]
    return current_start_val

def handle_nbs_id_column(meta_table, table_name, start_val_for_nbs_id):
    try:
        # Check if the column exists and if it's an identity column
        cur = meta_table.connection.cursor()
        cur.execute(f"SELECT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = '{table_name}' AND column_name = 'nbs_id');")
        column_exists = cur.fetchone()[0]

        if column_exists:
            # Get the corresponding sequence name dynamically
            sequence_name = get_sequence_name(cur, table_name)
            if not sequence_name:
                LOGGER.error(f"Sequence for column nbs_id not found.")
                return
            
            current_start_val = get_current_start_value(cur, sequence_name)
 
            if current_start_val != start_val_for_nbs_id:
                LOGGER.info(f"Column nbs_id is an identity column but its starting value is {current_start_val}. Setting the starting value to {start_val_for_nbs_id}...")
                cur.execute(f"ALTER SEQUENCE {sequence_name} RESTART WITH {start_val_for_nbs_id};")
                LOGGER.info(f"Starting value for column nbs_id set to {start_val_for_nbs_id} successfully.")
                meta_table.connection.commit()
            else:
                LOGGER.info("Column nbs_id is already an identity column with the specified starting value.")

        else:
            # Create the column and its sequence if do not exist
            LOGGER.info("Column nbs_id does not exist. Creating the column...")
            cur.execute(f"ALTER TABLE {table_name} ADD COLUMN nbs_id BIGINT GENERATED BY DEFAULT AS IDENTITY (START WITH {start_val_for_nbs_id});")
            LOGGER.info("Column nbs_id created successfully.")
            meta_table.connection.commit()

    except Exception as e:
        LOGGER.error(f"Failed to handle the scenarios for column nbs_id. Error: {e}")
        meta_table.connection.rollback()

if __name__ == '__main__':
    """ Run all local configs. """
    if len(sys.argv) == 1:
        config_path_file = 'fuse_config_path.config'
    else:
        config_path_file = sys.argv[2]

    config = {}
    config_file = cp.ConfigParser()
    config_file.read(config_path_file)
    sections = config_file.sections()
    for section in sections:
        config_file_section = config_file[section]
        for key in config_file_section:
            config[key] = config_file_section[key]
    config_path = config['config_path']

    summary_tables = run_configs(config_path, 'ignore_list.txt', create_summary_tables=False)

    # if len(summary_tables) != 0:
    #     regions = {summary_table.region for summary_table in summary_tables}
    #     source_types = {summary_table.source_type for summary_table in summary_tables}
    #     create_summary(config_filenames[0], regions, source_types, 'metadata_develop')
