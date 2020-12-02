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
		'area': Area,
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
	    "osm_id", "osm_way_id", "highway",  "area", "geometry"]]

	highway_shapes_df.to_csv(os.path.join(out_dir, "highway_squares.csv"))

# convert csv to shapefile
def shape_to_csv(dir_path, layer_name):
	in_csvfile = os.path.join(dir_path, layer_name)
	out_shpfile = os.path.join(dir_path, os.path.splitext(in_csvfile)[0])
	result = subprocess.run(["ogr2ogr", "-oo", "GEOM_POSSIBLE_NAMES=geometry*", "-f", "ESRI Shapefile", str(out_shpfile)+".shp", str(in_csvfile)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def csvLayerNames(layer_path):
	layer_names = glob.glob(os.path.join(layer_path, "*.csv"))
	xx = [os.path.basename(i) for i in layer_names]
	return xx


# simulate street light electrcity demand
def streetLightDemnd(input_path, input_path2, output_path):
    """Get standard load profile."""
    standardLoad = pd.read_csv(input_path)
    # SL1: All urban lights are operated as
    # - evening (16:15) - midnight (00:00) => 'On'
    # - midnight (00:15) - early morning  (05:15) => 'Off'
    #  - Early morning (05:30 )- morning (09:00) => 'On'
    SL1 = standardLoad['SB1'] / 1000
    # All urban roads are illuminated all night
    SL2 = standardLoad['SB2'] / 1000
    x_sl = 0.01464
    timestamp = standardLoad['Zeitstempel']
    # get planet OSM data for highway (line and polygon)
    osmLines = pd.read_csv(os.path.join(input_path2, 'highway_lines.csv'))
    osmSquares = pd.read_csv(os.path.join(input_path2, 'highway_squares.csv'))
    osmData = pd.concat([osmLines[["highway", "area"]], osmSquares[["highway", "area"]]], ignore_index=True)
    # for secenario 1 & 2
    osmArea = osmData["area"].sum()
    squaresArea = osmSquares["area"].sum()
    # scenario 3 => main roads only for all night
    mainRoads = ['living_street', 'motorway', 'pedestrian', 'primary',
                    'secondary', 'service', 'tertiary', 'trunk']
    mainRoadsArea = osmLines[osmLines["highway"].
                                        isin(mainRoads)]["area"].sum()
    # scenario 4 => minus the selected main roads
    # operated using SL1 (on and off)
    secondaryRoadsArea = osmLines[~osmLines["highway"].
                                            isin(mainRoads)]["area"].sum()
    mainRoadandsquareArea = mainRoadsArea + squaresArea

    _streetLoad_ = []
    fname = open(os.path.join(output_path,
                                "streetlight_load.csv"), 'w')
    fname.write('TIME;Mode2;Mode1;Mode3\n')
    for i, row in standardLoad.iterrows():
        row = str(timestamp[i]) + ';' + \
            str(osmArea*SL1[i]*x_sl) + ';' +\
            str(osmArea*SL2[i]*x_sl) + ';' +\
            str((mainRoadandsquareArea*SL2[i]*x_sl) + (secondaryRoadsArea*SL1[i]*x_sl)) + '\n'
        _streetLoad_.append(row)
    fname.writelines(_streetLoad_)
    fname.close()


def optimizationCommodities(input_path, input_path2, output_path):
    """Get solar data."""
    # demand and supply
    pv = pd.read_csv(input_path, parse_dates=True)
    wind = pd.read_csv(input_path2, parse_dates=True)

    wind['pv'] = pv['pv']
    feedin = wind
    df6 = pd.read_csv(os.path.join(
        output_path, "streetlight_load.csv"),
        delimiter=";", index_col='TIME', parse_dates=True)

    df6 = df6.iloc[0:len(wind.index)]  # TODO: fix index
    feedin['demand_Mode2'] = df6['Mode2'].values
    feedin['demand_Mode1'] = df6['Mode1'].values
    feedin['demand_Mode3'] = df6['Mode3'].values
    #feedin['demand_SB2_SB1'] = df6['SB2/SB1[kWh]'].values
    feedin.to_csv(os.path.join(
        output_path, 'optimization-commodities.csv'))
