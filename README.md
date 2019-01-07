# BlackGEM_tiling
Sky-tiling code for the BlackGEM array

Installation:
============
The simplest way to install the package is by using virtual environment. It has been tested using conda for Python3. 
It should also work in other virtual environments. 
In conda, create an environment:

>> conda create -n blackgem_tiling

Activate the environment:

>> source activate blackgem_tiling

(Install python.app to make sure that you can use pythonw otherwise matplotlib import will complain in conda).

Install this package:

>> python setup.py install

This should install the tiling code in your system. To run the code you will beed a sky-map fits file. One example
fits file will be includes in the repo. 

Running the module:
==================

In the python interpreter window, import the module:

>> from BlackGEM_tiling import ranked_tiling

Create a tile object:

>> tileObj = ranked_tiling.RankedTileGenerator()

Create the set of ranked tiles for a given sky-map:

>> ranked_tiles = tileObj.getRankedTiles('bayestar.fits.gz')

The return is a pandas DataFrame of the following form:

rank |  tile_index | tile_prob

- - - - - - - - - - - - - 

0   |  13921.0 |  0.058268

1   |    620.0 |  0.054954

2   |  14064.0 |  0.048635

3   |    536.0 |  0.038242

4   |    621.0 |  0.026141

...



