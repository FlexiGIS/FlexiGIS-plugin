# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Flexigis
                                 A QGIS plugin
 This is a test plugin for flexigis
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2020-09-16
        git sha              : $Format:%H$
        copyright            : (C) 2020 by chinonso
        email                : chinonso.unaichi@dlr.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
import os.path
import subprocess
import os
from pathlib import Path

from qgis.PyQt.QtCore import QSettings, QTranslator, QCoreApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction, QFileDialog, QMessageBox
from qgis.core import Qgis

import webbrowser
from PyQt5.Qt import QApplication, QUrl, QDesktopServices
import sys
from PyQt5.QtGui import *
from qgis.PyQt import QtGui

# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
from .flexigis_dialog import Filter_Dialog, Geoprocess_Dialog, Simulate_Dialog
from .flexigis_utils import (filter_pbf_with_poly, osm_convert, osm_filter,
                             osm_shapefile, csvLayerNames, symbolize_layer, simulate_urban_demand)
from .flexigis_utils import (filter_lines, filter_squares, shape_to_csv,
                             streetLightDemnd, optimizationCommodities, landuseLayers, building_layers, pv_feedin_generation)
from .flexigis_utils import refactor_landuse

class flexigis:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'flexigis{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)
            QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&flexigis')

        # Check if plugin was started the first time in current QGIS session
        # Must be set in initGui() to survive plugin reloads
        self.first_start = None

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('flexigis', message)

    def add_action(
            self,
            icon_path,
            text,
            callback,
            enabled_flag=True,
            add_to_menu=True,
            add_to_toolbar=True,
            status_tip=None,
            whats_this=None,
            parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            # Adds plugin icon to Plugins toolbar
            self.iface.addToolBarIcon(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/flexigis/flexigis.png'
        self.add_action(
            icon_path,
            text=self.tr(u''),
            callback=self.run,
            parent=self.iface.mainWindow())

        # will be set False in run()
        self.first_start = True

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&flexigis'),
                action)
            self.iface.removeToolBarIcon(action)

# ++++++++++++++++++++custom function begins++++++++++++++++++++++++++++++++++++
    # Block 1 => OSM File Filter
    def help_page(self):
        webbrowser.open('https://github.com/FlexiGIS/FlexiGIS-plugin/blob/master/flexigis/help/source/index.rst')

    def on_text_changed_b3(self):
        self.dlg1.b3.setEnabled(bool(self.dlg1.lineEdit1.text())
                                and bool(self.dlg1.lineEdit2.text()))

    def on_text_changed_b4(self):
        self.dlg1.b4.setEnabled(bool(self.dlg1.lineEdit3.text())
                                and bool(self.dlg1.lineEdit2.text()) and bool(self.dlg1.lineEdit1.text()))

    def select_pbf_block1(self):
        file_extentsion = "pbf(*.osm.pbf)"
        input_file = QFileDialog.getOpenFileName(
            self.dlg1, "Select input file", "", file_extentsion)
        filename_pbf = input_file[0]
        self.dlg1.lineEdit1.setText(filename_pbf)

    def select_polygon_block1(self):
        file_extentsion = "polygon(*.poly)"
        input_file = QFileDialog.getOpenFileName(
            self.dlg1, "Select input file", "", file_extentsion)
        filename_poly = input_file[0]
        self.dlg1.lineEdit2.setText(filename_poly)

    def on_b3_click(self):
        self.pbf_file = self.dlg1.lineEdit1.text()
        self.poly_file = self.dlg1.lineEdit2.text()
        output_name = "urban-infrastructure.osm.pbf"

        if os.path.isfile(self.pbf_file) and os.path.isfile(self.poly_file):
            _infile_ = os.path.splitext(self.pbf_file)[0]
            output_file_path = os.path.dirname(_infile_)
            self.dlg1.lineEdit3.setText(
                os.path.join(output_file_path, output_name))
        else:
            self.iface.messageBar().pushMessage(
                "Ensure input files (osm.pbf/poly) are valid", level=Qgis.Warning, duration=4)

    def popup_button(self, i):
        self.dlg1.lineEdit1.setText("")
        self.dlg1.lineEdit2.setText("")
        self.dlg1.lineEdit3.setText("")

    def on_b4_click(self):
        msg = QMessageBox()
        msg.setWindowTitle("OSM file filter")
        self.pbf_filtered = self.dlg1.lineEdit3.text()
        if Path(self.pbf_filtered).suffixes == Path(self.pbf_file).suffixes:
            if os.path.isfile(self.poly_file):
                filter_pbf_with_poly(
                    self.pbf_file, self.poly_file, self.pbf_filtered)
                self.iface.messageBar().pushMessage(
                    "osm.pbf data filter using polygon done! see output file {}".format(self.pbf_filtered), level=Qgis.Success, duration=3)

                # msg.setText("Done!")
                # msg.setIcon(QMessageBox.Information)
                # msg.setInformativeText(
                #     "osm.pbf file filter using bounding polygons completed.")
                
                # _ = msg.exec_()
            else:
                self.iface.messageBar().pushMessage(
                    "Ensure input poly file exits", level=Qgis.Warning, duration=4)
                msg.setText("File type error")
                msg.setIcon(QMessageBox.Warning)
                msg.setInformativeText(
                    "Use the browse button to choose a valid .poly file")
                msg.buttonClicked.connect(self.popup_button)
                _ = msg.exec_()
        else:
            self.iface.messageBar().pushMessage(
                "Output file extension error", level=Qgis.Warning, duration=4)
            msg.setText("Output file error")
            msg.setIcon(QMessageBox.Critical)
            msg.setInformativeText("Set output file extension to osm.pbf")
            msg.setDetailedText(
                "To properly set output file name, choose valid osm.pbf file and poly file using the browse file buttons and click Ok to set default Outputfile name.")
            msg.buttonClicked.connect(self.popup_button)
            _ = msg.exec_()

    # Block 2

    def popup_button_block2(self, i):
        self.dlg2.lineEdit1_2.setText("")
        self.dlg2.lineEdit2_2.setText("")

    def select_pbf_block2(self):
        file_extension = "pbf(*.osm.pbf)"
        input_file = QFileDialog.getOpenFileName(
            self.dlg2, "Select poly file", "", file_extension)
        self.urban_pbf = input_file[0]
        self.dlg2.lineEdit1_2.setText(self.urban_pbf)

    def select_shp_block2(self):
        file_extension = "shp(*.shp)"
        input_file = QFileDialog.getOpenFileName(
            self.dlg2, "Select shapefile", "", file_extension)
        self.shp_file = input_file[0]
        self.dlg2.lineEdit4_2.setText(self.shp_file)

    def select_shp2_block2(self):
        file_extension = "shp(*.shp)"
        input_file = QFileDialog.getOpenFileName(
            self.dlg2, "Select shapefile", "", file_extension)
        self.shp_file = input_file[0]
        self.dlg2.lineEdit5_2.setText(self.shp_file)

    def on_b6_click(self):
        msg = QMessageBox()
        msg.setWindowTitle("Info")
        if os.path.isfile(self.dlg2.lineEdit1_2.text()):
            self.out_dir = os.path.splitext(self.dlg2.lineEdit1_2.text())[0]
            self.out_dir = os.path.splitext(self.out_dir)[0]
            if Path(self.out_dir).exists():
                if self.dlg2.comboBox2_2.currentText() == "":
                    self.dlg2.comboBox2_2.addItems(
                        [layer for layer in csvLayerNames(self.out_dir)])
            else:
                os.mkdir(self.out_dir)
            self.dlg2.lineEdit2_2.setText(self.out_dir)
            self.dlg2.lineEdit3_2.setText(self.out_dir)
        else:
            self.iface.messageBar().pushMessage(
                "Input osm.pbf file error", level=Qgis.Warning, duration=4)
            msg.setText("Input file error")
            msg.setIcon(QMessageBox.Warning)
            msg.setInformativeText(
                "Use the browse button to choose a valid .osm.pbf file and click Ok")
            msg.buttonClicked.connect(self.popup_button_block2)
            _ = msg.exec_()

    def on_text_changed_b6(self):
        self.dlg2.b2_2.setEnabled(bool(self.dlg2.lineEdit1_2.text()))

    def on_text_changed_b7(self):
        self.dlg2.b3_2.setEnabled(bool(self.dlg2.lineEdit2_2.text()))


    def on_b3_2_click(self):
        msg = QMessageBox()
        msg.setWindowTitle("OSM geoprocessing")
        input_filename = self.dlg2.lineEdit1_2.text()  # Filename of PBF-File
        input_file_dirname = os.path.split(input_filename)[0]
        # out_file = os.path.basename(self.dlg.lineEdit_5.text())
        out_file_dirname = self.dlg2.lineEdit2_2.text()  # Directory of Outputs
        osm_tag = self.dlg2.comboBox1_2.currentText()  # -> "building", "highway", "landuse"
        outfile_tag = os.path.join(input_file_dirname, osm_tag)
        # landuse_file_tag = os.path.join(input_file_dirname, 'landuse')
        landuse_out_tag = os.path.join(input_file_dirname, 'landuse')
        landuse_input_file_name = self.dlg2.lineEdit4_2.text()  # Input of EO landuse data
        buildings_input_file_name = self.dlg2.lineEdit5_2.text()  # Input of EO building data

        #TODO: check out_file and Out_dir ==> make the naming less redundant!
        if osm_tag == "highway":
            highway_cv_filenames = ['highway_squares.csv', 'highway_lines.csv']
            osm_convert(input_filename, out_file_dirname)
            osm_filter(out_file_dirname, osm_tag, outfile_tag)
            osm_shapefile(outfile_tag)
            filter_lines(os.path.join(
                outfile_tag, "lines.shp"), out_file_dirname)
            filter_squares(os.path.join(
                outfile_tag, "multipolygons.shp"), out_file_dirname)
    #         # filter_highway_points(os.path.join(
    #         #     osm_tag, "points.shp"), os.path.join(out_file, "highway_points"))
            if self.dlg2.comboBox2_2.currentText() == "":
                self.dlg2.comboBox2_2.addItems(
                    [layer for layer in csvLayerNames(out_file_dirname)])
            elif all(item in [self.dlg2.comboBox2_2.itemText(i) for i in range(self.dlg2.comboBox2_2.count())] for item in highway_cv_filenames):
                pass
            else:
                self.dlg2.comboBox2_2.addItems(
                    [layer for layer in csvLayerNames(out_file_dirname) if layer in highway_cv_filenames])
            self.iface.messageBar().pushMessage(
                "Highway data geoprocessing complete.", level=Qgis.Success, duration=4)

        elif osm_tag == "landuse":
            landuse_cv_filenames = ['landuse.csv']
            if self.dlg2.checkBox_3.isChecked():  # Only process OSM data?
                osm_convert(input_filename, out_file_dirname)
                osm_filter(out_file_dirname, osm_tag, outfile_tag)
                osm_shapefile(outfile_tag)
                landuseLayers(os.path.join(
                    outfile_tag, "multipolygons.shp"), out_file_dirname)
            else:  # Process EO data
                refactor_landuse(os.path.join(landuse_input_file_name), out_file_dirname, osm_tag)

            # Add generated files to comboBox for later post-processing to shapefiles etc
            if self.dlg2.comboBox2_2.currentText() == "":  # if no elements available -> add them
                self.dlg2.comboBox2_2.addItems(
                    [layer for layer in csvLayerNames(out_file_dirname)])  # only 'landuse.csv' should appear here
            elif "landuse.csv" in [self.dlg2.comboBox2_2.itemText(i) for i in range(self.dlg2.comboBox2_2.count())]:
                pass
            else:
                self.dlg2.comboBox2_2.addItems(
                    [layer for layer in csvLayerNames(out_file_dirname) if layer in landuse_cv_filenames])
            self.iface.messageBar().pushMessage(
                "Landuse data geoprocessing complete.", level=Qgis.Success, duration=4)

        elif osm_tag == "building":
            building_cv_filenames = ['agricultural.csv', 'commercial.csv',
                                        'educational.csv', 'industrial.csv', 'residential.csv']
            if self.dlg2.checkBox_3.isChecked():  # Only process OSM data?
                osm_convert(input_filename, out_file_dirname)
                osm_filter(out_file_dirname, "building", outfile_tag)
                osm_shapefile(outfile_tag)

                # landuse
                osm_convert(input_filename, out_file_dirname)
                osm_filter(out_file_dirname, "landuse", landuse_file_tag)
                osm_shapefile(landuse_file_tag)

                building_layers(os.path.join(outfile_tag, "multipolygons.shp"),
                                os.path.join(landuse_file_tag, "multipolygons.shp"), out_file_dirname)
            else:  # Process EO data only
                # Generate landuse-data from input and export to shapefile
                refactor_landuse(os.path.join(landuse_input_file_name), out_file_dirname, osm_tag)
                shape_to_csv(out_file_dirname, 'landuse.csv')

                # intersect landuse and building layers
                building_layers(buildings_input_file_name,
                                os.path.join(out_file_dirname, "landuse.shp"), out_file_dirname)


            # Add generated files to comboBox for later post-processing to shapefiles etc
            if self.dlg2.comboBox2_2.currentText() == "":
                self.dlg2.comboBox2_2.addItems(
                    [layer for layer in csvLayerNames(out_file_dirname)])
            elif all(item in [self.dlg2.comboBox2_2.itemText(i) for i in range(self.dlg2.comboBox2_2.count())] for item in building_cv_filenames):
                pass
            else:
                self.dlg2.comboBox2_2.addItems(
                    [layer for layer in csvLayerNames(out_file_dirname) if layer in building_cv_filenames])
            self.iface.messageBar().pushMessage(
                "Building categories geoprocessing complete.", level=Qgis.Success, duration=4)
        else:
            self.iface.messageBar().pushMessage(
                "Filter/Geoprocessing for this tag is under developement", level=Qgis.Info, duration=4)


    def checkBox_click(self):
        if self.dlg2.checkBox.isChecked():
            self.dlg2.checkBox_2.setChecked(False)
        if self.dlg2.comboBox2_2.currentText() != "":
            self.dlg2.b4_2.setEnabled(bool(self.dlg2.lineEdit3_2.text()))

    def checkBox2_click(self):
        if self.dlg2.checkBox_2.isChecked():
            self.dlg2.checkBox.setChecked(False)
        if self.dlg2.comboBox2_2.currentText() != "":
            self.dlg2.b4_2.setEnabled(bool(self.dlg2.lineEdit3_2.text()))

    def on_click_b4_2(self):
        if self.dlg2.checkBox.isChecked():
            dir_name = self.dlg2.lineEdit3_2.text()
            layer_name = self.dlg2.comboBox2_2.currentText()
            shape_to_csv(dir_name, layer_name)
            self.iface.messageBar().pushMessage("shapefile successfully generated! see output in {}".format(
                str(dir_name)), level=Qgis.Success, duration=4)
        else:
            dir_name = self.dlg2.lineEdit3_2.text()
            layer_name = self.dlg2.comboBox2_2.currentText()
            symbolize_layer(dir_name, layer_name)
            # create_mapLayout(dir_name)
            self.iface.messageBar().pushMessage("Map layers successfully generated! see output in {}".format(
                str(dir_name)), level=Qgis.Success, duration=4)

    # Window 3 >>

    def selectSLP(self):
        file_extentsion = "csv(*.csv)"
        input_file = QFileDialog.getOpenFileName(
            self.dlg3, "Select input file", "", file_extentsion)
        stl_csv = input_file[0]
        self.dlg3.lineEdit1_3.setText(stl_csv)

    def selectPV(self):
        file_extentsion = "csv(*.csv)"
        input_file = QFileDialog.getOpenFileName(
            self.dlg3, "Select input file", "", file_extentsion)
        stl_csv = input_file[0]
        self.dlg3.lineEdit2_3.setText(stl_csv)

    def selectWind(self):
        file_extentsion = "csv(*.csv)"
        input_file = QFileDialog.getOpenFileName(
            self.dlg3, "Select input file", "", file_extentsion)
        stl_csv = input_file[0]
        self.dlg3.lineEdit3_3.setText(stl_csv)

    def selectLayer_dir(self):
        input_dir = str(QFileDialog.getExistingDirectory(
            self.dlg3, "Select directory"))
        #stl_csv = input_file
        self.dlg3.lineEdit5_3.setText(input_dir)

    def on_textchanged_l1_3(self):
        self.dlg3.b4_3.setEnabled(bool(self.dlg3.lineEdit1_3.text()))

    def on_textchanged_l4_3(self):
        self.dlg3.b5_3.setEnabled(bool(self.dlg3.lineEdit4_3.text()))

    def on_click_b4_3(self):
        slp_csv = self.dlg3.lineEdit1_3.text()
        if os.path.isfile(slp_csv):
            dir_name = os.path.split(slp_csv)[0]
            self.dlg3.lineEdit4_3.setText(dir_name)

        else:
            self.iface.messageBar().pushMessage(
                "Input field not a valid file path. Select a valid standard load profile csv file", level=Qgis.Critical, duration=4)

    # simulate energy demand for street light
    def simulateDemand_streetLight(self):
        slp_path = self.dlg3.lineEdit1_3.text()
        wind_path = self.dlg3.lineEdit2_3.text()
        pv_path = self.dlg3.lineEdit3_3.text()
        in_dir = self.dlg3.lineEdit4_3.text()
        layer_dir = self.dlg3.lineEdit5_3.text()
        folder_name = "demand_profile"
        new_dir = os.path.join(in_dir, folder_name)

        if Path(new_dir).exists():
            pass
        else:
            os.mkdir(new_dir)
        if self.dlg3.comboBox1_3.currentText() == "Street light elect. demand":
            if os.path.isfile(slp_path) and os.path.isdir(layer_dir):
                if pv_path == "" or wind_path == "":
                    streetLightDemnd(slp_path, layer_dir, new_dir)
                    self.iface.messageBar().pushMessage(
                        "Street light elect. demand simulation done!", level=Qgis.Success, duration=4)
                elif os.path.isfile(pv_path) and os.path.isfile(wind_path):
                    streetLightDemnd(slp_path, layer_dir, new_dir)
                    optimizationCommodities(pv_path, wind_path, new_dir)
                    self.iface.messageBar().pushMessage(
                        "Street light elect. demand simulation done!", level=Qgis.Success, duration=4)
                else:
                    self.iface.messageBar().pushMessage(
                        "Ensure selected file paths exit!", level=Qgis.Critical, duration=4)
            else:
                self.iface.messageBar().pushMessage(
                    "Ensure selected file paths and Layer directory path exit!", level=Qgis.Critical, duration=4)

        elif self.dlg3.comboBox1_3.currentText() == "Urban infrastructure elect. demand":
            if os.path.isfile(slp_path) and os.path.isdir(layer_dir):
                if pv_path == "" or wind_path == "":
                    simulate_urban_demand(slp_path, layer_dir, new_dir)
                    self.iface.messageBar().pushMessage(
                        "Urban building electricity demand simulation done!", level=Qgis.Success, duration=4)
                elif os.path.isfile(pv_path) and os.path.isfile(wind_path):
                    simulate_urban_demand(slp_path, layer_dir, new_dir)
                    pv_feedin_generation(
                        pv_path, wind_path, layer_dir, new_dir)
                    self.iface.messageBar().pushMessage(
                        "Urban building electricity demand and renewable generation simulation done!", level=Qgis.Success, duration=4)
                else:
                    self.iface.messageBar().pushMessage(
                        "Ensure selected file paths exit!", level=Qgis.Critical, duration=4)
            else:
                self.iface.messageBar().pushMessage(
                    "Ensure selected file paths and Layer directory path exit!", level=Qgis.Critical, duration=4)

    # window Navigations << >>
    def close_widget(self):
        self.dlg1.close()
        self.dlg2.close()
        self.dlg3.close()

    def open_widget2(self):
        self.dlg1.close()
        self.dlg2.show()

    def open_widget3(self):
        self.dlg2.close()
        self.dlg3.show()

    def previous_widget1(self):
        self.dlg2.close()
        self.dlg1.show()

    def previous_widget2(self):
        self.dlg3.close()
        self.dlg2.show()

    def run(self):
        """Run method that performs all the real work"""
        # Create the dialog with elements (after translation) and keep reference
        # Only create GUI ONCE in callback, so that it will only load when the plugin is started
        # window 1
        if self.first_start == True:
            self.first_start = False
            self.dlg1 = Filter_Dialog()

            self.dlg1.b3.setEnabled(False)
            self.dlg1.b4.setEnabled(False)
            self.dlg1.b5.setEnabled(False)

            self.dlg1.b1.clicked.connect(self.select_pbf_block1)
            self.dlg1.b2.clicked.connect(self.select_polygon_block1)
            self.dlg1.lineEdit3.textChanged.connect(
                self.on_text_changed_b4)
            self.dlg1.lineEdit1.textChanged.connect(
                self.on_text_changed_b3)
            self.dlg1.lineEdit2.textChanged.connect(
                self.on_text_changed_b3)

            self.dlg1.b3.clicked.connect(self.on_b3_click)
            self.dlg1.b4.clicked.connect(self.on_b4_click)
            self.dlg1.b8.clicked.connect(self.help_page)

        # window 2
        self.dlg2 = Geoprocess_Dialog()
        self.dlg2.b2_2.setEnabled(False)
        self.dlg2.b3_2.setEnabled(False)
        self.dlg2.b4_2.setEnabled(False)
        self.dlg2.lineEdit1_2.textChanged.connect(self.on_text_changed_b6)
        self.dlg2.lineEdit2_2.textChanged.connect(self.on_text_changed_b7)
        self.dlg2.b1_2.clicked.connect(self.select_pbf_block2)
        self.dlg2.b2_2.clicked.connect(self.on_b6_click)
        self.dlg2.b3_2.clicked.connect(self.on_b3_2_click)
        self.dlg2.b9_2.clicked.connect(self.select_shp_block2)
        self.dlg2.b10_2.clicked.connect(self.select_shp2_block2)
        self.dlg2.checkBox.stateChanged.connect(self.checkBox_click)
        self.dlg2.checkBox_2.stateChanged.connect(self.checkBox2_click)
        self.dlg2.b4_2.clicked.connect(self.on_click_b4_2)
        self.dlg2.b5_2.clicked.connect(self.help_page)

        # window 3
        self.dlg3 = Simulate_Dialog()
        self.dlg3.b4_3.setEnabled(False)
        self.dlg3.b5_3.setEnabled(False)
        self.dlg3.b8_3.setEnabled(False)
        self.dlg3.b1_3.clicked.connect(self.selectSLP)
        self.dlg3.b2_3.clicked.connect(self.selectPV)
        self.dlg3.b3_3.clicked.connect(self.selectWind)
        self.dlg3.b10_3.clicked.connect(self.selectLayer_dir)
        self.dlg3.lineEdit1_3.textChanged.connect(self.on_textchanged_l1_3)
        self.dlg3.b4_3.clicked.connect(self.on_click_b4_3)
        self.dlg3.lineEdit4_3.textChanged.connect(self.on_textchanged_l4_3)
        self.dlg3.b5_3.clicked.connect(self.simulateDemand_streetLight)
        self.dlg3.b6_3.clicked.connect(self.help_page)

        # widgets navigations
        self.dlg1.b7.clicked.connect(self.close_widget)
        self.dlg1.b6.clicked.connect(self.open_widget2)
        self.dlg2.b7_2.clicked.connect(self.open_widget3)
        self.dlg2.b6_2.clicked.connect(self.previous_widget1)
        self.dlg2.b8_2.clicked.connect(self.close_widget)
        self.dlg3.b9_3.clicked.connect(self.close_widget)
        self.dlg3.b7_3.clicked.connect(self.previous_widget2)

        # show the dialog
        self.dlg1.show()
        # Run the dialog event loop
        result = self.dlg1.exec_()
        # See if OK was pressed

        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass
