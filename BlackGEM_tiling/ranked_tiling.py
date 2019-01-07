# Copyright (C) 2019 Shaon Ghosh
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import os
import sys
import pickle
from pkg_resources import resource_filename

import numpy as np
import pandas as pd
import healpy as hp


def getSkymapData(skymapFilename, res=256):
    healpix_skymap = hp.read_map(skymapFilename, verbose=False)
    skymapUD = hp.ud_grade(healpix_skymap, res, power=-2)
    healpix_skymap = skymapUD
    npix = len(healpix_skymap)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, np.arange(0, npix))
    ra_pixels = np.rad2deg(phi)
    dec_pixels = np.rad2deg(0.5*np.pi - theta)
    data_pixels = healpix_skymap[np.arange(0, npix)]
    skymapdata = [np.arange(npix), ra_pixels, dec_pixels, data_pixels]
    return skymapdata

class RankedTileGenerator:
    def __init__(self, precomputedTileFile=None,
                 tilefile=None):
        if precomputedTileFile is None:
            precomputedTileFile = resource_filename('BlackGEM_tiling',
                                                    'preComputed_256_NN.dat')
        File = open(precomputedTileFile, 'rb')

        if tilefile is None:
            tilefile = resource_filename('BlackGEM_tiling',
                                         'BlackGEMtilesPG_AV.dat')
        self.tileData = pickle.load(File)
        self.IDs = np.array(list(self.tileData.keys()))
        ID, RA_cent, Dec_cent, gal_AV = np.loadtxt(tilefile, unpack=True)
        ID = ID.astype(int)
        tile_values = np.vstack((RA_cent, Dec_cent, gal_AV)).T
        self.tileDict = dict(zip(ID, tile_values))


    def getRankedTiles(self, verbose=False, fitsfilename=None,
                       skymapdata=False, res=256):
        if fitsfilename:
            skymap = getSkymapData(fitsfilename, res=res)
        elif skymapdata:
            skymap = [skymapdata[0], skymapdata[1],
                      skymapdata[2], skymapdata[3]]

        self.pixel_id_all = skymap[0]
        self.point_ra_all = skymap[1]
        self.point_dec_all = skymap[2]
        self.point_pVal_all = skymap[3]
        pTile = []
        TileProbSum = 0.0
        if verbose:
            print('Computing Ranked-Tiles...')

        tileList = list(self.tileData.keys())

        for tile in tileList:
            pixelsInThisTile = self.tileData[tile]
            pvalThesePix = self.point_pVal_all[pixelsInThisTile]
            pTile.append(np.sum(pvalThesePix))

        pTile = np.array(pTile)
        sorted_indices = np.argsort(-1*pTile)
        output = np.vstack((self.IDs[sorted_indices],
                            pTile[sorted_indices])).T

        self.df = pd.DataFrame(output, columns=['tile_index', 'tile_prob'])
        return self.df


    def plotTiles(self, CI=0.9):
        from .allskymap_basic import AllSkyMap
        import pylab as pl

        ranked_tiles = self.df
        include_CI = np.cumsum(ranked_tiles.tile_prob.values) < CI
        include_CI[np.sum(include_CI)] = True
        ranked_tiles_CI = ranked_tiles[include_CI]
        tileIDs_CI = ranked_tiles_CI.tile_index
        ra_tiles_CI = []
        dec_tiles_CI = []
        gal_av_tiles_CI = []
        for id in tileIDs_CI:
            [ra_tile, dec_tile, gal_av] = self.tileDict[id]
            ra_tiles_CI.append(ra_tile)
            dec_tiles_CI.append(dec_tile)
            gal_av_tiles_CI.append(gal_av)

        ra_tiles_CI = np.array(ra_tiles_CI)
        dec_tiles_CI = np.array(dec_tiles_CI)

        order = np.argsort(-self.point_pVal_all)
        ra = self.point_ra_all[order]
        dec = self.point_dec_all[order]
        pVal = self.point_pVal_all [order]
        include = np.cumsum(pVal) < CI
        include[np.sum(include)] = True
        ra_CI = ra[include]
        dec_CI = dec[include]
        pVal_CI = pVal[include]

        pl.rcParams.update({'font.size': 16})
        m = AllSkyMap(projection='hammer')
        RAP_map, DecP_map = m(ra_CI, dec_CI)
        RAP_center, DecP_center = m(ra_tiles_CI, dec_tiles_CI)
        m.drawparallels(np.arange(-90.,120.,20.), color='grey',
                        labels=[False,True,True,False],
                        labelstyle='+/-')

        m.drawmeridians(np.arange(0.,420.,30.), color='grey')
        m.drawmapboundary(fill_color='white')
        lons = np.arange(-150,151,30)
        m.label_meridians(lons, fontsize=16, vnudge=1, halign='left', hnudge=-1)
        m.scatter(RAP_map, DecP_map, c='yellow', s=3, alpha=0.8)
        m.scatter(RAP_center, DecP_center, c='red', s=10)

        pl.savefig('skymap_with_tiles.png')
