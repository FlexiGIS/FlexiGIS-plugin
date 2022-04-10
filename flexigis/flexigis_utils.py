import glob
import logging
import subprocess
import os
from pathlib import Path
import pandas as pd
import numpy as np
from qgis.core import *

# from qgis.PyQt.QtCore import QVariant
import qgis
from PyQt5.QtGui import *
from qgis.PyQt import QtGui
from random import randrange

# from qgis.PyQt.QtGui import QPolygonF, QColor

log_path = Path(__file__).parent
if log_path.exists():
    pass
else:
    log_path.mkdir()
log_file = log_path / "flexigis_plugin.log"
# create a log file
logging.basicConfig(
    format="%(asctime)s: %(levelname)s: %(message)s",
    filename=log_file,
    level=logging.DEBUG,
)


def filter_pbf_with_poly(pbffile, polyfile, output_filename):
    """
    This function uses Osmosis tool to extract an area of interest from a pbf
    file, provided the poly file (polygon) of the destrict or area of interest is specified.

    input: .osm.pbf of a city e.g. Bremen_latest.osm.pbf
    input: .poly file of a destrict in Bremen e.g Neustadt.ploy
    output: this returns .osm.pbf of Neustadt
    """
    try:
        result = subprocess.run(
            [
                "osmosis",
                "--read-pbf",
                "file=" + str(pbffile),
                "--tag-filter",
                "accept-ways",
                "building=*",
                "--used-node",
                "--bounding-polygon",
                "file=" + str(polyfile),
                "--buffer",
                "outPipe.0=building",
                "--read-pbf",
                "file=" + str(pbffile),
                "--tag-filter",
                "accept-ways",
                "highway=*",
                "--used-node",
                "--bounding-polygon",
                "file=" + str(polyfile),
                "--buffer",
                "outPipe.0=highway",
                "--read-pbf",
                "file=" + str(pbffile),
                "--tag-filter",
                "accept-ways",
                "landuse=*",
                "--used-node",
                "--bounding-polygon",
                "file=" + str(polyfile),
                "--buffer",
                "outPipe.0=landuse_1",
                "--read-pbf",
                "file=" + str(pbffile),
                "--tag-filter",
                "accept-relations",
                "landuse=*",
                "--used-node",
                "--bounding-polygon",
                "file=" + str(polyfile),
                "--buffer",
                "outPipe.0=landuse_2",
                "--merge",
                "inPipe.0=landuse_1",
                "inPipe.1=landuse_2",
                "--buffer",
                "outPipe.0=landuse_all",
                "--merge",
                "inPipe.0=landuse_all",
                "inPipe.1=building",
                "--buffer",
                "outPipe.0=landuse_building",
                "--merge",
                "inPipe.0=landuse_building",
                "inPipe.1=highway",
                "--write-pbf",
                "file=" + str(output_filename),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        logging.info(result)
    except RuntimeError as e:
        logging.error("Osmosis RuntimeError!: {}".format(e))


def osm_convert(in_file, out_file):
    """
    This function converts .osm.pbf file into .o5m file.
    input: .osm.pdf
    output: .o5m
    """
    result = subprocess.run(
        ["osmconvert", str(in_file), "-o=" + str(out_file) + ".o5m"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    logging.info(result)


def osm_filter(input_o5m, osm_tag, output_osm):
    """
    The o5m file is filtered using OSMfilter tool, which extracts features informations
    of specified Urban-infrastructure Tag e.g higway building, landuse
    """
    result = subprocess.run(
        [
            "osmfilter",
            str(input_o5m) + ".o5m",
            "--keep=" + str(osm_tag),
            "-o=" + str(output_osm) + ".osm",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    logging.info(result)


def osm_shapefile(output_osm):
    """
    Export the filtered osm Tag file to shape files of polygon, lines and points.

    """
    result = subprocess.run(
        [
            "ogr2ogr",
            "-skipfailures",
            "-f",
            "ESRI Shapefile",
            str(output_osm),
            str(output_osm) + ".osm",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    logging.info(result)


def compute_area(dataset, width):
    """Compute area for each line feature and return a dataframe object.

    dataset: Dataframe object
    width: Dictionary, unique highway category as key and the width in meters
    as value.

    """
    dataset["area"] = np.nan
    for key, value in width.items():
        dataset.loc[key]["area"] = dataset.loc[key]["length"] * value
        # Area.append(area)
    # Area = pd.concat(Area)
    # dataset["area"] = Area.values
    return dataset


def filter_lines(shp_file, out_dir):
    """ "
    Filter, coordinate geoprocess, calculate width (meters) and area (meters squared) of urban highway lines
    and export result to csv file

    input: shape file e.g Lines.shp
    output: csv file of geoprrocessed highway information
    """
    layer = QgsVectorLayer(shp_file, "", "ogr")
    highway_lines = {
        "living_street",
        "motorway",
        "residential",
        "pedestrian",
        "primary",
        "secondary",
        "service",
        "tertiary",
        "trunk",
        "unclassified",
        "motorway_link",
        "primary_link",
        "secondary_link",
        "tertiary_link",
        "trunk_link",
    }
    width = [
        7.5,
        15.50,
        6.5,
        7.5,
        10.5,
        9.5,
        7.5,
        9.5,
        9.5,
        7.5,
        6.5,
        6.5,
        6.5,
        6.5,
        6.5,
    ]

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
        if feature["highway"] in highway_lines:
            geom = feature.geometry()
            geom.transform(tr)
            your_string = geom.asWkt()
            length_ = round(geom.length(), 2)
            geom2.append(your_string)
            length.append(length_)

    df_line = pd.DataFrame(
        {
            "osm_id": osm_id,
            "highway": highway,
            "length": length,
            "geometry": geom2,
        }
    )

    check_features = set(df_line.highway.unique()).intersection(highway_lines)
    width_dict = {k: highway_width[k] for k in check_features}
    df_line = df_line.set_index("highway")
    df_line = compute_area(df_line, width_dict)
    df_line = df_line.reset_index()
    df_line = df_line[["highway", "osm_id", "length", "area", "geometry"]]
    df_line = df_line.sort_values("highway")
    df_line.to_csv(os.path.join(out_dir, "highway_lines.csv"))


# *************** Highway square ***********************


def filter_squares(shp_file, out_dir):
    """ "
    Filter, coordinate geoprocess, calculate area (meters squared) of urban highway squares
    and export result to csv file

    input: shape file e.g Polygons.shp
    output: csv file of geoprrocessed highway information
    """
    layer = QgsVectorLayer(shp_file, "", "ogr")
    square_feature = {
        "crossing",
        "footway",
        "living_street",
        "pedestrian",
        "platform",
        "residential",
        "service",
        "traffic_island",
    }
    osm_id = [feature["osm_id"] for feature in layer.getFeatures()]
    osm_way_id = [feature["osm_way_id"] for feature in layer.getFeatures()]

    # get available line features from tag
    other_tags = [feature["other_tags"] for feature in layer.getFeatures()]

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
            "osm_id": osm_id,
            "osm_way_id": osm_way_id,
            "other_tags": other_tags,
            "area": Area,
            "geometry": geom2,
        }
    )

    # filter other_tags attribute to extract highway square features
    highway_shapes = df_square[df_square["other_tags"] != NULL]
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
    highway_shapes_df = highway_shapes.loc[
        highway_shapes["highway"].isin(square_feature)
    ]
    highway_shapes_df = highway_shapes_df.sort_values("highway")
    highway_shapes_df = highway_shapes_df.reset_index()
    highway_shapes_df = highway_shapes_df[
        ["osm_id", "osm_way_id", "highway", "area", "geometry"]
    ]
    highway_shapes_df.to_csv(os.path.join(out_dir, "highway_squares.csv"))


# convert csv to shapefile


def shape_to_csv(dir_path, layer_name):
    """
    Export result stored as csv layers to shape files
    input: csv file
    output: shapefile
    """
    in_csvfile = os.path.join(dir_path, layer_name)
    out_shpfile = os.path.join(dir_path, os.path.splitext(in_csvfile)[0])
    result = subprocess.run(
        [
            "ogr2ogr",
            "-oo",
            "GEOM_POSSIBLE_NAMES=geometry*",
            "-f",
            "ESRI Shapefile",
            str(out_shpfile) + ".shp",
            str(in_csvfile),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def csvLayerNames(layer_path):
    """
    Get names of layers stored as csv files.
    """
    layer_names = glob.glob(os.path.join(layer_path, "*.csv"))
    xx = [os.path.basename(i) for i in layer_names]
    return xx


# simulate street light electrcity demand
def streetLightDemnd(input_path, input_path2, output_path):
    """Calculate's street light electricity demand, using highway information such as:
    Highway Lines area (metersquared),
    Highway squares area (metersquared),
    street light electricity usage index,
    street light standard load profile.

    also the renewable generation can be simulated by parsing
    wind and solar PV time series data.

    input: all input data are in csv file format
    output: timeseries of electricity demand in csv.

    """
    standardLoad = pd.read_csv(input_path)
    # SL1: All urban lights are operated as
    # - evening (16:15) - midnight (00:00) => 'On'
    # - midnight (00:15) - early morning  (05:15) => 'Off'
    #  - Early morning (05:30 )- morning (09:00) => 'On'
    SL1 = standardLoad["SB1"] / 1000
    # All urban roads are illuminated all night
    SL2 = standardLoad["SB2"] / 1000
    x_sl = 0.01464
    timestamp = standardLoad["Zeitstempel"]
    # get planet OSM data for highway (line and polygon)
    osmLines = pd.read_csv(os.path.join(input_path2, "highway_lines.csv"))
    osmSquares = pd.read_csv(os.path.join(input_path2, "highway_squares.csv"))
    osmData = pd.concat(
        [osmLines[["highway", "area"]], osmSquares[["highway", "area"]]],
        ignore_index=True,
    )
    # for secenario 1 & 2
    osmArea = osmData["area"].sum()
    squaresArea = osmSquares["area"].sum()
    # scenario 3 => main roads only for all night
    mainRoads = [
        "living_street",
        "motorway",
        "pedestrian",
        "primary",
        "secondary",
        "service",
        "tertiary",
        "trunk",
    ]
    mainRoadsArea = osmLines[osmLines["highway"].isin(mainRoads)]["area"].sum()
    # scenario 4 => minus the selected main roads
    # operated using SL1 (on and off)
    secondaryRoadsArea = osmLines[~osmLines["highway"].isin(mainRoads)][
        "area"
    ].sum()
    mainRoadandsquareArea = mainRoadsArea + squaresArea

    _streetLoad_ = []
    fname = open(os.path.join(output_path, "streetlight_load.csv"), "w")
    fname.write("TIME;Mode2;Mode1;Mode3\n")
    for i, row in standardLoad.iterrows():
        row = (
            str(timestamp[i])
            + ";"
            + str(osmArea * SL1[i] * x_sl)
            + ";"
            + str(osmArea * SL2[i] * x_sl)
            + ";"
            + str(
                (mainRoadandsquareArea * SL2[i] * x_sl)
                + (secondaryRoadsArea * SL1[i] * x_sl)
            )
            + "\n"
        )
        _streetLoad_.append(row)
    fname.writelines(_streetLoad_)
    fname.close()


def optimizationCommodities(input_path, input_path2, output_path):
    """Get renewable generation data."""
    # demand and supply
    pv = pd.read_csv(input_path, parse_dates=True)
    wind = pd.read_csv(input_path2, parse_dates=True)

    wind["pv"] = pv["pv"]
    feedin = wind
    df6 = pd.read_csv(
        os.path.join(output_path, "streetlight_load.csv"),
        delimiter=";",
        index_col="TIME",
        parse_dates=True,
    )

    df6 = df6.iloc[0 : len(wind.index)]  # TODO: fix index
    feedin["demand_Mode2"] = df6["Mode2"].values
    feedin["demand_Mode1"] = df6["Mode1"].values
    feedin["demand_Mode3"] = df6["Mode3"].values
    # feedin['demand_SB2_SB1'] = df6['SB2/SB1[kWh]'].values
    feedin.to_csv(os.path.join(output_path, "optimization-commodities.csv"))


def simulate_urban_demand(
    input_destination1, input_destination2, output_destination
):
    """Simulate urban building infrastructure electricity demand for different infrastructure
    clusters. The urban buildings are clustered into 5 different categories
    1. Agricultural
    2. Commercial
    3. Educational
    4. Industrial
    5. Residential
    output: csv file of electricity demand.
    """

    dfa = pd.read_csv(os.path.join(input_destination2, "agricultural.csv"))
    dfc = pd.read_csv(os.path.join(input_destination2, "commercial.csv"))
    dfe = pd.read_csv(os.path.join(input_destination2, "educational.csv"))
    dfi = pd.read_csv(os.path.join(input_destination2, "industrial.csv"))
    dfr = pd.read_csv(os.path.join(input_destination2, "residential.csv"))
    # get standaerd load profile
    df_SLP = pd.read_csv(input_destination1)
    norm_ = 1000

    # SLP * area * electricity usage index
    load_a = (df_SLP["L0"] / norm_) * dfa.area.sum() * (120 / norm_)
    load_i = dfi.area.sum() * (645 / norm_) * df_SLP["G3"] / norm_
    load_c = dfc.area.sum() * (201 / norm_) * df_SLP["G0"] / norm_
    load_e = dfe.area.sum() * (142 / norm_) * (df_SLP["G1"] / norm_)
    load_r = dfr.area.sum() * (df_SLP["H0"] / norm_) * (146 / norm_)

    demand_ = pd.DataFrame(
        {
            "time": df_SLP.Zeitstempel,
            "agricultural": load_a,
            "commercial": load_c,
            "educational": load_e,
            "industrial": load_i,
            "residential": load_r,
        }
    )
    demand_.to_csv(
        os.path.join(output_destination, "urban_building_elect_demand.csv")
    )


def pv_feedin_generation(
    input_destination1,
    input_destination2,
    input_destination3,
    output_destination,
):
    """
    Calculate renewable generation of the building. The solar PV roof top generation is calculated.
    """
    """Get wind and pv feed in data."""

    dfa = pd.read_csv(os.path.join(input_destination3, "agricultural.csv"))
    dfc = pd.read_csv(os.path.join(input_destination3, "commercial.csv"))
    dfe = pd.read_csv(os.path.join(input_destination3, "educational.csv"))
    dfi = pd.read_csv(os.path.join(input_destination3, "industrial.csv"))
    dfr = pd.read_csv(os.path.join(input_destination3, "residential.csv"))

    area_a, area_c = dfa["area"].sum(), dfc["area"].sum()
    area_e, area_i = dfe["area"].sum(), dfi["area"].sum()
    area_r = dfr["area"].sum()

    pv_peak_power = 206.7989
    pv_module_area = 1.646
    norm_ = 1000

    pv_feedin = pd.read_csv(input_destination1, parse_dates=True)
    wind_feedin = pd.read_csv(input_destination2, parse_dates=True)

    pv_feedin = pv_feedin.reindex()
    time_stamp = pv_feedin[["time"]]
    pv_data = pv_feedin["pv"]
    wind_power_data = wind_feedin["wind"]

    # roof top area for pv installation
    areaRoof_a, areaRoof_c = area_a * 0.267, area_c * 0.267
    areaRoof_e, areaRoof_i = area_e * 0.267, area_i * 0.267
    areaRoof_r = area_r * 0.578

    # calculate aggregate peak power at rooftop
    total_roof_to_area = (
        areaRoof_a + areaRoof_c + areaRoof_e + areaRoof_i + areaRoof_r
    )
    speak_power_agg = (
        pv_peak_power * (total_roof_to_area / pv_module_area)
    ) / norm_

    # """Simulate Agricultural building type electricity demand and PV feedin supply."""
    pv_roofTop_generation = time_stamp
    pv_roofTop_generation["agricultural"] = (
        pv_data * pv_peak_power * (areaRoof_a / pv_module_area)
    ) / norm_
    pv_roofTop_generation["commercial"] = (
        pv_data * pv_peak_power * (areaRoof_c / pv_module_area)
    ) / norm_
    pv_roofTop_generation["educational"] = (
        pv_data * pv_peak_power * (areaRoof_e / pv_module_area)
    ) / norm_
    pv_roofTop_generation["industrial"] = (
        pv_data * pv_peak_power * (areaRoof_i / pv_module_area)
    ) / norm_
    pv_roofTop_generation["residential"] = (
        pv_data * pv_peak_power * (areaRoof_r / pv_module_area)
    ) / norm_
    pv_roofTop_generation["wind generation"] = wind_power_data

    pv_roofTop_generation.to_csv(
        os.path.join(output_destination, "pv_rooftop_generation.csv")
    )


# ********************export layer to map**********************


def symbolize_layer(dir_path, layerfile_name):
    """Categorize and symbolize map in QGIS desktop.

    When the layer is exported as map, the urban infrastructure features
    are symbolized with different color code.
    """
    layer_path = os.path.join(dir_path, layerfile_name)
    file_name = os.path.splitext(layer_path)
    layer_name = os.path.basename(file_name[0])
    uri = (
        "file://"
        + layer_path
        + "?type=csv&detectTypes=yes&wktField=geometry&crs=epsg:3857&spatialIndex=no&subsetIndex=no&watchFile=no"
    )
    vLayer = QgsVectorLayer(uri, str(layer_name), "delimitedtext")

    # display layer on qgis desktop
    QgsProject.instance().addMapLayer(vLayer)

    # get unique layer
    layer = qgis.utils.iface.activeLayer()
    field_names = [field.name() for field in layer.fields()]

    if "highway" in field_names:
        unique_layer = layer.fields().indexFromName("highway")
    elif "landuse" in field_names:
        unique_layer = layer.fields().indexFromName("landuse")

    elif "category" in field_names:
        unique_layer = layer.fields().indexFromName("category")

    unique_values = layer.uniqueValues(unique_layer)
    # fill categories
    categories = []
    for unique_value in unique_values:
        # initialize the default symbol for this geometry
        symbol = QgsSymbol.defaultSymbol(layer.geometryType())

        # configure symbol layer
        layer_style = {}
        layer_style["color"] = "%d, %d, %d" % (
            randrange(0, 256),
            randrange(0, 256),
            randrange(0, 256),
        )
        layer_style["outline"] = "#000000"
        symbol_layer = QgsSimpleFillSymbolLayer.create(layer_style)

        # replace default symbol with the configured one
        if symbol_layer is not None:
            symbol.changeSymbolLayer(0, symbol_layer)

        # create renderer object
        category = QgsRendererCategory(unique_value, symbol, str(unique_value))
        # entry for the list of category items
        categories.append(category)

    # create renderer object
    if "highway" in field_names:
        renderer = QgsCategorizedSymbolRenderer("highway", categories)
    elif "landuse" in field_names:
        renderer = QgsCategorizedSymbolRenderer("landuse", categories)

    elif "category" in field_names:
        renderer = QgsCategorizedSymbolRenderer("category", categories)

    # assign the created renderer to the layer
    if renderer is not None:
        layer.setRenderer(renderer)

    layer.triggerRepaint()


# ++++++ filter building and landuse


def transformGeo(x):
    # helper function, transforms coordinate reference system from 4326 to 3857
    crsSrc = QgsCoordinateReferenceSystem(4326)
    crsDrc = QgsCoordinateReferenceSystem(3857)
    tr = QgsCoordinateTransform(crsSrc, crsDrc, QgsProject.instance())
    x.transform(tr)
    return x


def transformLanduseGeo1(geom):
    # for landuse subroutine
    crsSrc = QgsCoordinateReferenceSystem(25832)
    crsDrc = QgsCoordinateReferenceSystem(3857)
    tr = QgsCoordinateTransform(crsSrc, crsDrc, QgsProject.instance())
    geom.transform(tr)
    return geom

def transformLanduseGeo2(geom):
    # for building subroutine
    crsSrc = QgsCoordinateReferenceSystem(25832)  # bei CLC epsg:25832
    crsDrc = QgsCoordinateReferenceSystem(4326)
    tr = QgsCoordinateTransform(crsSrc, crsDrc, QgsProject.instance())
    geom.transform(tr)
    return geom


def computeArea(x):
    # helper function, computes area of a QGIS feature coordinate
    area = round(x.area(), 2)
    return area


def geoToAswkt(x):
    # helper function, converts georeference to Aswkt
    x = x.asWkt()
    return x


def landuseLayers(landuseShapefile, out_dir):
    """geoprocess landuse layer Tag and export as csv file."""
    landuseLayers = get_landuseLayers(landuseShapefile)
    landuseLayers["geometry"] = landuseLayers.loc[:, "geometry"].map(
        transformGeo
    )
    landuseLayers["area"] = landuseLayers.loc[:, "geometry"].map(computeArea)
    landuseLayers["geometry"] = landuseLayers.loc[:, "geometry"].map(geoToAswkt)
    landuseLayers.to_csv(os.path.join(out_dir, "landuse.csv"))


def get_landuseLayers(landuseShapefile):
    """Extract relevant landuse Tag information from shapefile."""
    landuseLayer = QgsVectorLayer(landuseShapefile, "", "ogr")
    # extract landuse features from layers
    osm_id_ = [feature["osm_id"] for feature in landuseLayer.getFeatures()]
    osm_way_id = [
        feature["osm_way_id"] for feature in landuseLayer.getFeatures()
    ]
    landuse = [feature["landuse"] for feature in landuseLayer.getFeatures()]

    # replace NULL osm_id with the corresponding non NULL values from osm_way_id list
    osm_id = []
    for i in range(len(osm_id_)):
        if osm_id_[i] == NULL:
            osm_id_[i] = osm_way_id[i]
        osm_id.append(osm_id_[i])

    # calculate geometry area
    geometry_ = []
    area = []
    for feature in landuseLayer.getFeatures():
        geom = feature.geometry()
        geometry_.append(geom)

    # create a dataframe of cleaned landuse layers
    landuseLayers = pd.DataFrame(
        {"osm_id": osm_id, "landuse": landuse, "geometry": geometry_}
    )
    landuseLayers = landuseLayers.sort_values(by="landuse")
    cluster = [
        "farmyard",
        "greenhouse_horticulture",
        "farmland",
        "commercial",
        "retail",
        "industrial",
        "residential",
    ]
    landuseLayers = landuseLayers[landuseLayers["landuse"].isin(cluster)]
    return landuseLayers


def refactor_landuse(landuseShapefile, out_dir, flag='landuse'):
    '''Extract relevant landuse Tag information from shapefile.'''
    landuseLayer = QgsVectorLayer(landuseShapefile, "", 'ogr')
    # extract landuse features from layers
    osm_id_ = [feature['id'] for feature in landuseLayer.getFeatures()]
    landuse = [feature['landcover'] for feature in landuseLayer.getFeatures()]

    # replace NULL osm_id with the corresponding non NULL values from osm_way_id list
    #osm_id = []
    #for i in range(len(osm_id_)):
    #    if osm_id_[i] == NULL:
    #        osm_id_[i] = 1
    #    osm_id.append(osm_id_[i])

    # append geometry data to dataframe
    geometry_ = []
    for feature in landuseLayer.getFeatures():
        geom = feature.geometry()
        geometry_.append(geom)

    # create a dataframe of cleaned landuse layers
    landuseLayers = pd.DataFrame(
        {
            'osm_id': osm_id_,
            'landuse': landuse,
            'geometry': geometry_
        }
    )

    landuse_map = {
        'Pasture, meadows and other perma\nnent grasslands under agricultural\nuse': 'agricultural',
        'Peatbogs': 'agricultural',
        'Fruit tree and berry plantations': 'agricultural',
        'Agricultural farms': 'agricultural',
        'Broad-leaved forest': 'agricultural',
        'Mixed forest': 'agricultural',
        'Non-irrigated arable land': 'agricultural',
        'Coniferous forest': 'agricultural',
        'Water bodies': 'agricultural',
        'Water courses': 'agricultural',
        'Airports': 'commercial',
        'Sport and leisure facilities': 'commercial',
        'Port area': 'industrial',
        'Road and rail networks and associ\nated land': 'industrial',
        'Mineral extraction sites': 'industrial',
        'Dump sites': 'industrial',
        'Green urban area': 'residential',
        'Continuous urban fabric': 'residential',
        'Discontinuous urban fabric': 'residential'
    }

    #landuseLayers = landuseLayers.sort_values(by='landuse')
    landuseLayers['landuse'] = pd.Series([landuse_map.get(n, n) for n in landuseLayers['landuse']])

    cluster = ['commercial', 'agricultural', 'industrial', 'residential']
    landuseLayers = landuseLayers[landuseLayers['landuse'].isin(cluster)]

    if flag == 'landuse':
        landuseLayers['geometry'] = landuseLayers.loc[:,
                            'geometry'].map(transformLanduseGeo1)
    if flag == 'building':
        landuseLayers['geometry'] = landuseLayers.loc[:,
                                    'geometry'].map(transformLanduseGeo2)
    landuseLayers['area'] = landuseLayers.loc[:, 'geometry'].map(computeArea)
    landuseLayers['geometry'] = landuseLayers.loc[:,
                                'geometry'].map(geoToAswkt)
    landuseLayers.to_csv(os.path.join(out_dir, 'landuse.csv'))


def get_buildingLayers(buildingShapefile):
    """Preprocess building data. The building infrastructure are clustered into 5 different categories
    1. Agricultural
    2. Commercial
    3. Educational
    4. Industrial
    5. Residential

    """
    buildingLayer = QgsVectorLayer(buildingShapefile, "", "ogr")
    # extract landuse features from layers

    osm_id_ = [feature["osm_id"] for feature in buildingLayer.getFeatures()]
    osm_way_id = [
        feature["osm_way_id"] for feature in buildingLayer.getFeatures()
    ]
    building = [feature["building"] for feature in buildingLayer.getFeatures()]

    # replace NULL osm_id with the corresponding non NULL values from osm_way_id list
    osm_id = []
    for i in range(len(osm_id_)):
        if osm_id_[i] == NULL:
            osm_id_[i] = osm_way_id[i]
        osm_id.append(osm_id_[i])

    # change coordinate refrence system
    # calculate geometry area
    geometry_ = []
    area = []
    for feature in buildingLayer.getFeatures():
        geom = feature.geometry()
        geometry_.append(geom)

    # create a dataframe of cleaned landuse layers
    buildingLayers = pd.DataFrame(
        {"osm_id": osm_id, "building": building, "geometry": geometry_}
    )
    buildingLayers = buildingLayers.sort_values(by="building")
    return buildingLayers


def getIntersections(df_landuse, df_building):
    """Helper function for extracting intersect geometries between polygons,
    e.g a building type inside a it's landuse polygons.
    """
    geoms = []
    for lg in df_landuse.geometry:
        for bg in df_building.geometry:
            if bg.intersects(lg):
                geoms.append(bg)

    resolution = df_building[df_building.geometry.isin(geoms)]
    return resolution


def layersBuildings(landuse, building, type="residential"):
    # find intersection for different buildiing features with their respective landuse type

    # residential, commercial, industrial and agricultural
    if type == "residential":
        residentials = getIntersections(landuse, building)
        return residentials
    elif type == "commercial":
        commercials = getIntersections(landuse, building)
        return commercials
    elif type == "agricultural":
        agriculturals = getIntersections(landuse, building)
        return agriculturals
    elif type == "industrial":
        industrials = getIntersections(landuse, building)
        return industrials


def education(buildings_, out_dir):
    # Educational

    # buildings_ = get_buildingLayers(buildingShapefile)
    educational = buildings_[
        buildings_.building.isin(["kindergarten", "school", "university"])
    ]
    educational["geometry"] = educational["geometry"].map(transformGeo)
    educational["area"] = educational["geometry"].map(computeArea)
    educational["geometry"] = educational["geometry"].map(geoToAswkt)
    educational["category"] = ["educational" for i in range(len(educational))]
    educational.to_csv(os.path.join(out_dir, "educational.csv"))


def residential_building(buildings_, landuseShapefile, out_dir):
    # Residential buildings

    # buildings_ = get_buildingLayers(buildingShapefile)
    buildings = buildings_[
        ~buildings_.building.isin(["kindergarten", "school", "university"])
    ]
    landuses = get_landuseLayers(landuseShapefile)
    landuseResidential = landuses[landuses.landuse == "residential"]
    residential = layersBuildings(
        landuseResidential, buildings, type="residential"
    )
    residential["geometry"] = residential["geometry"].map(transformGeo)
    residential["area"] = residential["geometry"].map(computeArea)
    residential["geometry"] = residential["geometry"].map(geoToAswkt)
    residential["category"] = ["residential" for i in range(len(residential))]
    residential.to_csv(os.path.join(out_dir, "residential.csv"))


def commercial_building(buildings_, landuseShapefile, out_dir):
    # Commercial buildings

    # buildings_ = get_buildingLayers(buildingShapefile)
    buildings = buildings_[
        ~buildings_.building.isin(["kindergarten", "school", "university"])
    ]
    landuses = get_landuseLayers(landuseShapefile)
    landusecomm = landuses[landuses.landuse.isin({"commercial", "retail"})]
    commercial = layersBuildings(landusecomm, buildings, type="commercial")
    commercial["geometry"] = commercial["geometry"].map(transformGeo)
    commercial["area"] = commercial["geometry"].map(computeArea)
    commercial["geometry"] = commercial["geometry"].map(geoToAswkt)
    commercial["category"] = ["commercial" for i in range(len(commercial))]
    commercial.to_csv(os.path.join(out_dir, "commercial.csv"))


def agricultural_building(buildings_, landuseShapefile, out_dir):
    # Agricultural buildings

    # buildings_ = get_buildingLayers(buildingShapefile)
    buildings = buildings_[
        ~buildings_.building.isin(["kindergarten", "school", "university"])
    ]
    landuses = get_landuseLayers(landuseShapefile)
    land_agric = landuses[
        landuses.landuse.isin(
            {"farmyard", "farmland", "greenhouse_horticulture"}
        )
    ]
    agricultural = layersBuildings(land_agric, buildings, type="agricultural")
    agricultural["geometry"] = agricultural["geometry"].map(transformGeo)
    agricultural["area"] = agricultural["geometry"].map(computeArea)
    agricultural["geometry"] = agricultural["geometry"].map(geoToAswkt)
    agricultural["category"] = [
        "agricultural" for i in range(len(agricultural))
    ]
    agricultural.to_csv(os.path.join(out_dir, "agricultural.csv"))


def industrial_building(buildings_, landuseShapefile, out_dir):
    # Industrial buildings

    # buildings_ = get_buildingLayers(buildingShapefile)
    buildings = buildings_[
        ~buildings_.building.isin(["kindergarten", "school", "university"])
    ]
    landuses = get_landuseLayers(landuseShapefile)
    landuseIndustrial = landuses[landuses.landuse == "industrial"]
    industrial = layersBuildings(
        landuseIndustrial, buildings, type="industrial"
    )
    industrial["geometry"] = industrial["geometry"].map(transformGeo)
    industrial["area"] = industrial["geometry"].map(computeArea)
    industrial["geometry"] = industrial["geometry"].map(geoToAswkt)
    industrial["category"] = ["industrial" for i in range(len(industrial))]
    industrial.to_csv(os.path.join(out_dir, "industrial.csv"))


def building_layers(buildingShapefile, landuseShapefile, out_dir):
    # generate clusters of different features, out put are csv files
    try:
        buildings_ = get_buildingLayers(buildingShapefile)
        education(buildings_, out_dir)
        residential_building(buildings_, landuseShapefile, out_dir)
        commercial_building(buildings_, landuseShapefile, out_dir)
        agricultural_building(buildings_, landuseShapefile, out_dir)
        industrial_building(buildings_, landuseShapefile, out_dir)

        logging.info("Buildings categorization complete")
    except RuntimeError as e:
        logging.error(e)


# def create_mapLayout(out_path):
# 	layers = {layer for layer in QgsProject.instance().mapLayers().values()}
# 	project = QgsProject.instance()
# 	manager = project.layoutManager()
# 	layoutName = 'layout1'
# 	layouts_list = manager.printLayouts()

# 	# remove duplicate layouts
# 	for layout in layouts_list:
# 		if layout.name() == layoutName:
# 			manager.removeLayout(layout)

# 	layout = QgsPrintLayout(project)
# 	layout.initializeDefaults()
# 	layout.setName(layoutName)
# 	manager.addLayout(layout)

# 	# create map item in the layout
# 	map = QgsLayoutItemMap(layout)
# 	map.setRect(20, 20, 20, 20)

# 	# set map extent
# 	ms = QgsMapSettings()
# 	for i in layers:
# 		ms.setLayers([i])  # layers to be mapped
# 	rect = QgsRectangle(ms.fullExtent())
# 	rect.scale(1.0)
# 	ms.setExtent(rect)
# 	map.setExtent(rect)
# 	map.setBackgroundColor(QColor(255, 255, 255, 0))
# 	layout.addLayoutItem(map)
# 	map.attemptMove(QgsLayoutPoint(23.626, 7.091, QgsUnitTypes.LayoutMillimeters))
# 	map.attemptResize(QgsLayoutSize(180, 180, QgsUnitTypes.LayoutMillimeters))

# 	# add legend and scale bar
# 	legend = QgsLayoutItemLegend(layout)
# 	legend.setTitle('')
# 	layerTree = QgsLayerTree()
# 	for i in layers:
# 		layerTree.addLayer(i)
# 	legend.model().setRootGroup(layerTree)
# 	layout.addLayoutItem(legend)
# 	legend.attemptMove(QgsLayoutPoint(230, 15, QgsUnitTypes.LayoutMillimeters))

# 	scalebar = QgsLayoutItemScaleBar(layout)
# 	scalebar.setStyle('Double Box')
# 	scalebar.setUnits(QgsUnitTypes.DistanceKilometers)
# 	scalebar.setNumberOfSegments(5)
# 	scalebar.setNumberOfSegmentsLeft(0)
# 	scalebar.setUnitsPerSegment(5)
# 	scalebar.setLinkedMap(map)
# 	scalebar.setUnitLabel('km')
# 	scalebar.setFont(QFont('Arial', 14))
# 	scalebar.update()
# 	layout.addLayoutItem(scalebar)
# 	scalebar.attemptMove(QgsLayoutPoint(
# 		53.143, 187.091, QgsUnitTypes.LayoutMillimeters))
# 	scalebar.attemptResize(QgsLayoutSize(
# 		134.753, 12.500, QgsUnitTypes.LayoutMillimeters))

# export layout to png or pdf
# layout = manager.layoutByName(layoutName)
# exporter = QgsLayoutExporter(layout)
# map_name = os.path.join(out_path,'map.png')
# exporter.exportToImage(map_name, QgsLayoutExporter.ImageExportSettings())
# exporter.exportToPdf(map_name, QgsLayoutExporter.PdfExportSettings())
