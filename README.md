# FlexiGIS-plugin: A QGIS GUI interface for FlexigIS-light.


## Installation: testing the plugin on Qgis

- Clone the plugin repository (since the repo is not public yet, you will be ask for the flexigis github account username and password)

```console
terminal:~$ https://github.com/FlexiGIS/FlexiGIS-plugin.git
```

- cd into the FlexiGIS-plugin directory and copy the  *flexgi_test* folder to yours QGIS3 python path.

``` 
terminal/FlexiGIS-plugin:~$ cp -r flexigi_test /home/user/.local/share/QGIS/QGIS3/profiles/default/python/plugins
```

After following the above steps, start your QGIS3 deskop environment. First install plugin reloader, this is used to reload any installed plugin in QGIS. To reload the flexgi_test plugin, ensure you can see it on the list of installed plugins in your QGIS3 (Plugin management section).

In your QGIS menu bar go to: plugins > Manage and install plugins > installed (to check installed plugin). To ensure that the newly created plugin is present in the QGIS plugin path.

To reload plugin go to the plugin reloader icon, click on the drop down menu, select configuration, then select the plugin you want to reload. In our case we select the 'flexgi_test'. The final hit the reload icon of the plugin reloader. This is done any time the plugin backed or front end codes are edited.

## Dependencies

- python >=3.5
- pandas
- gdal 
- osmfilter
- osmosis

*NOTE: you can check your QGIS3 python version by running the below commands on your QGIS python terminal

````
import sys

sys.version
````
