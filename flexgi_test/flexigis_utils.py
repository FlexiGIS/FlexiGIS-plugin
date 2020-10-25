import glob
import logging
import subprocess
import os
from pathlib import Path
import pandas as pd
import numpy as np
from qgis.core import *
from qgis.PyQt.QtCore import QVariant
import qgis


# if Path("log").exists():
#     pass
# else:
#     os.mkdir("log")
# log_file = os.path.join("log", "flexigis_plugin.log")
#     # create a log file
# logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s',
#                     filename=log_file,
#                     level=logging.DEBUG)


def filter_pbf_with_poly(pbffile, polyfile, output_filename):
    try:
        result = subprocess.run(["osmosis", "--read-pbf", "file="+str(pbffile),
            "--tag-filter", "accept-ways", "building=*", "--used-node",
                        "--bounding-polygon", "file=" +
                            str(polyfile), "--buffer", "outPipe.0=building",
                        "--read-pbf", "file=" +
                            str(pbffile), "--tag-filter", "accept-ways", "highway=*", "--used-node",
                        "--bounding-polygon", "file=" +
                            str(polyfile), "--buffer", "outPipe.0=highway", "--read-pbf", "file="+str(pbffile),
                        "--tag-filter", "accept-ways", "landuse=*", "--used-node", "--bounding-polygon", "file=" +
                            str(polyfile),
                        "--buffer", "outPipe.0=landuse_1", "--read-pbf", "file=" +
                            str(pbffile), "--tag-filter", "accept-relations",
                        "landuse=*", "--used-node", "--bounding-polygon", "file=" +
                            str(polyfile), "--buffer", "outPipe.0=landuse_2",
                        "--merge", "inPipe.0=landuse_1", "inPipe.1=landuse_2", "--buffer", "outPipe.0=landuse_all",
                        "--merge", "inPipe.0=landuse_all", "inPipe.1=building", "--buffer", "outPipe.0=landuse_building",
                        "--merge", "inPipe.0=landuse_building", "inPipe.1=highway", "--write-pbf", "file=" + str(output_filename)],
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
# logging.info(result)
    except RuntimeError:
        pass
# logging.error("Osmosis RuntimeError!")


def osm_convert(in_file, out_file):
    result = subprocess.run(["osmconvert", str(in_file), "-o=" + str(out_file) + ".o5m"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
# logging.info(result)


def osm_filter(input_o5m, osm_tag, output_osm):
    result = subprocess.run(["osmfilter", str(input_o5m)+".o5m",
                             "--keep=" + str(osm_tag), "-o=" + str(output_osm) + ".osm"],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # logging.info(result)


def osm_shapefile(output_osm):
    result = subprocess.run(["ogr2ogr", "-skipfailures", "-f",
                             "ESRI Shapefile", str(output_osm), str(output_osm) + ".osm"],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # logging.info(result)


def compute_area(dataset, width):
    """Compute area for each line feature and return a dataframe object.

    dataset: Dataframe object
    width: Dictionary, unique highway category as key and the width in meters
    as value.

    """
    Area = []
    dataset = dataset.set_index(["highway"])
    for key, value in width.items():
        area = dataset.loc[key]["length"]*value
        Area.append(area)
    Area = pd.concat(Area)
    dataset["area"] = Area.values
    dataset_new = dataset.reset_index()
    return dataset_new.sort_values("highway")


def filter_lines(shp_file, out_dir):
	layer = QgsVectorLayer(shp_file, "", 'ogr')
	highway_lines = {'living_street', 'motorway', 'pedestrian', 'primary', 'secondary', 'service', 'tertiary', 'trunk',
                  'motorway_link', 'primary_link', 'secondary_link',
                  'tertiary_link', 'trunk_link'}
	width = [7.5, 15.50, 7.5, 10.5, 9.5, 7.5, 9.5, 9.5, 6.5, 6.5, 6.5, 6.5, 6.5]

	highway_width = dict(zip(highway_lines, width))

	osm_id = [
	    feature["osm_id"]
	    for feature in layer.getFeatures()
	    if feature["highway"] in highway_lines
	]

	highway = [
	    feature["highway"]
	    for feature in layer.getFeatures()
	    if feature["highway"] in highway_lines
	]

	crsSrc = QgsCoordinateReferenceSystem(4326)
	crsDrc = QgsCoordinateReferenceSystem(3857)
	tr = QgsCoordinateTransform(crsSrc, crsDrc, QgsProject.instance())
	geom2 = []
	length = []

	for feature in layer.getFeatures():
		if feature['highway'] in highway_lines:
			geom = feature.geometry()
			geom.transform(tr)
			your_string = geom.asWkt()
			length_ = round(geom.length(), 2)
			geom2.append(your_string)
			length.append(length_)

	df_line = pd.DataFrame(
	    {
		'osm_id': osm_id,
		'highway': highway,
		'length': length,
		'geometry': geom2
	    }
	)

	df_line = compute_area(df_line, highway_width)
	df_line = df_line.reset_index()
	df_line = df_line[[
            "highway", "osm_id", "length", "area", "geometry"]]
	df_line.to_csv(os.path.join(out_dir, "highway_lines.csv"))

#*************** Highway square ***********************


def filter_squares(shp_file, out_dir):
	layer = QgsVectorLayer(shp_file, "", 'ogr')
	square_feature = {'crossing', 'footway', 'living_street',
                   'pedestrian', 'platform', 'residential',
                   'service', 'traffic_island'
                   }
	osm_id = [
            feature["osm_id"]
            for feature in layer.getFeatures()
        ]

	osm_way_id = [
		feature["osm_way_id"]
		for feature in layer.getFeatures()
	]

	# get available line features from tag
	other_tags = [
		feature["other_tags"]
		for feature in layer.getFeatures()
	]

	crsSrc = QgsCoordinateReferenceSystem(4326)
	crsDrc = QgsCoordinateReferenceSystem(3857)
	tr = QgsCoordinateTransform(crsSrc, crsDrc, QgsProject.instance())
	geom2 = []
	Area = []

	for feature in layer.getFeatures():
		geom = feature.geometry()
		geom.transform(tr)
		your_string = geom.asWkt()
		area_ = round(geom.area(), 2)
		geom2.append(your_string)
		Area.append(area_)

	df_square = pd.DataFrame(
	    {
		'osm_id': osm_id,
		'osm_way_id': osm_way_id,
		'other_tags': other_tags,
		'Area': Area,
		'geometry': geom2
	    }
	)

	# filter other_tags attribute to extract highway square features
	highway_shapes = df_square[df_square["other_tags"] !=NULL]
    # create a clean highway tag column
	tag_val = []
	tag = "highway"
	for i in highway_shapes.other_tags.values:
		if tag in i:
			tag_true = 1
			tag_val.append(tag_true)
		else:
			tag_false = 0
			tag_val.append(tag_false)
	highway_shapes["bool"] = tag_val
	highway_shapes = highway_shapes[highway_shapes["bool"] == 1]
	filter_tag = []
	for n in highway_shapes["other_tags"].values:
		ll = n.split(",")
		for lll in ll:
			if lll.startswith('"highway"'):
				filter_tag.append(lll)
	feature_list = []
	for i in filter_tag:
		m = i.split('"')
		feature_list.append(m[3])
	highway_shapes["highway"] = feature_list
	highway_shapes_df = highway_shapes.loc[highway_shapes["highway"].isin(
	    square_feature)]
	highway_shapes_df = highway_shapes_df.sort_values("highway")
	highway_shapes_df = highway_shapes_df.reset_index()
	highway_shapes_df = highway_shapes_df[[
	    "osm_id", "osm_way_id", "highway",  "Area", "geometry"]]

	highway_shapes_df.to_csv(os.path.join(out_dir, "highway_squares.csv"))

def shape_to_csv(shapefile_dir):
    dir_path = shapefile_dir
    dir_names = glob.glob(os.path.join(dir_path, "*"))

    for i in dir_names:
        j = gpd.read_file(os.path.join(i, str(os.path.basename(i))+".shp"))
        j.to_csv(os.path.join(
            i, str(os.path.basename(i))+".csv"), encoding="utf-8")

def csvLayerNames(layer_path):
	layer_names = glob.glob(os.path.join(layer_path, "*.csv"))
	xx = [os.path.basename(i) for i in layer_names]
	return xx
