# Copyright (C) 2018 Shaon Ghosh
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

from __future__ import print_function

import os
import sys
import time
import pickle
import argparse
from collections import defaultdict

import numpy as np
import pylab as pl
import healpy as hp

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import constants as c, units as u

parser = argparse.ArgumentParser(description='Pre-compute the tiles')
parser.add_argument('--nside', type=int, help='Value of nside')
parser.add_argument('--tilefile', type=str, help='Name of the tile file')
parser.add_argument('--outpath', type=str, default='.',
                    help='Path to output file')
parser.add_argument('--method', type=str, default='NT',
                    help="'NT' (=nearest tile) or 'area'")
parser.add_argument('--FOV', type=float, default=2.7, help="BlackGEM's FOV")
parser.add_argument('--protect', action='store_true',
                    help='Set True to write-protect the output file')

args = parser.parse_args()

nside = args.nside
npix = hp.nside2npix(nside)
pixel_index = np.arange(npix)
RAs, Decs = hp.pix2ang(nside, pixel_index, lonlat=True, nest=False)

insideorNot = np.zeros(npix, dtype='bool')

filename = '_'.join(['preComputed', str(nside), args.method]) + '.dat'
filename = os.path.join(args.outpath, filename)
File = open(filename, 'wb')
ID, RA_cent, Dec_cent, gal_AV = np.loadtxt(args.tilefile, unpack=True)
tileCents = SkyCoord(RA_cent, Dec_cent, frame='icrs', unit='deg')

coveredRAs = np.array([])
coveredDecs = np.array([])
coveredPixels = np.array([])
TimeRecord = []
numPoints = []
output = {}

if args.method == 'NT':
    start = time.time()
    pixtiledict = {}
    for pix_id, pixRa, pixDec in zip(pixel_index, RAs, Decs):
        pixCoord = SkyCoord(pixRa, pixDec, unit='deg')
        sep = tileCents.separation(pixCoord)
        minsep = np.argmin(sep)
        nearestTile = tileCents[minsep]
        nearestTileID = ID.astype(int)[minsep]
        pixtiledict[pix_id] = nearestTileID

    # Inverting the pixel tile association dictionary #
    tilepixdict = defaultdict(list)
    for key,value in pixtiledict.items():
        tilepixdict[value].append(key)
    tilepixdict = dict(tilepixdict)

    end = time.time()
    print('Time taken to run loop = {}'.format(end - start))

pickle.dump(tilepixdict, File)
File.close()


# dec_interval = np.sqrt(FOV)
# for thisID, thisRA, thisDec, in zip(ID, RA_cent, Dec_cent):
#     print('Declination angle = {}'.format(thisDec))
