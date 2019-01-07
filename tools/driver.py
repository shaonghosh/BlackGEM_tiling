import numpy as np
import ranked_tiling as rt

tileObj = rt.RankedTileGenerator()
bayestarFile = '/Users/ghosh4/Downloads/2016_fits/320494/bayestar.fits.gz'
ranked_tiles = tileObj.getRankedTiles(fitsfilename=bayestarFile)
include90 = np.cumsum(ranked_tiles.tile_prob.values) < 0.9
include90[np.sum(include90)] = True
tileObj.plotTiles()
