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
import multiprocessing
from collections import defaultdict

import ray
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
parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='Set true for stdout statements')
args = parser.parse_args()


@ray.remote
def loop_through_pixels(pixel_index, RAs, Decs):
    pixtiledict = {}
    for pix_id, pixRa, pixDec in zip(pixel_index, RAs, Decs):
        pixCoord = SkyCoord(pixRa, pixDec, unit='deg')
        sep = tileCents.separation(pixCoord)
        minsep = np.argmin(sep)
        nearestTile = tileCents[minsep]
        nearestTileID = ID.astype(int)[minsep]
        pixtiledict[pix_id] = nearestTileID
    return pixtiledict

ID, RA_cent, Dec_cent, gal_AV = np.loadtxt(args.tilefile, unpack=True)
tileCents = SkyCoord(RA_cent, Dec_cent, frame='icrs', unit='deg')

coveredRAs = np.array([])
coveredDecs = np.array([])
coveredPixels = np.array([])
TimeRecord = []
numPoints = []
output = {}

if args.verbose:
    ray.init(logging_level=1)
else:
    ray.init(logging_level=40)
cores = multiprocessing.cpu_count()
if args.verbose:
    print("Total number of cores in this machine: {}".format(cores))

workers = cores
nside = args.nside
npix = hp.nside2npix(nside)
split = np.array_split(np.arange(npix), workers)
futures = []
start = time.time()
for this_pix_set in split:
    RAs, Decs = hp.pix2ang(nside, this_pix_set, lonlat=True, nest=False)
    futures.append(loop_through_pixels.remote(this_pix_set, RAs, Decs))

ray.get(futures)

pixtiledict = {}
for future in futures:
    pixtiledict.update(ray.get(future))

# Inverting the pixel tile association dictionary #
tilepixdict = defaultdict(list)
for key,value in pixtiledict.items():
    tilepixdict[value].append(key)
tilepixdict = dict(tilepixdict)

end = time.time()
if args.verbose:
    print('Time taken to run loop = {}'.format(end - start))

filename = '_'.join(['preComputed', str(nside), "NN"]) + '.dat'
filename = os.path.join(args.outpath, filename)
with open(filename, "wb") as f:
    pickle.dump(tilepixdict, f)


