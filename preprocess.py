# Copyright (C) 2020 Jiahua He

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import mrcfile
import argparse
import numpy as np
from tqdm import tqdm

import interp3d

def get_args():
    parser = argparse.ArgumentParser(description="This script generates voxels in .npz format from .mrc file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mrc", "-m", type=str, required=True,
                        help="Input .mrc file")
    parser.add_argument("--npz", "-n", type=str, required=True,
                        help="Output .npz file")
    parser.add_argument("--apix", "-a", type=float, default=1.0,
                        help="Grid interval in angstroms")
    parser.add_argument("--threshold", "-t", type=float, default=0.0,
                        help="Density threshold level (default is 0 for simulated maps)")
    parser.add_argument("--ignorestart", "-i", type=bool, default=False,
                        help="Whether ignore start position, True for simulated map and False of experimental map")
    args = parser.parse_args()
    return args

def getmode(mapc, mapr, maps):
    if mapc == 1 and mapr == 2 and maps == 3:
        ordermode = 1
    elif mapc == 1 and mapr == 3 and maps == 2:
        ordermode = 2
    elif mapc == 2 and mapr == 1 and maps == 3:
        ordermode = 3
    elif mapc == 2 and mapr == 3 and maps == 1:
        ordermode = 4
    elif mapc == 3 and mapr == 1 and maps == 2:
        ordermode = 5
    elif mapc == 3 and mapr == 2 and maps == 1:
        ordermode = 6
    else:
        print("Input file gives malformed dimension ordering!")
        sys.exit()
    return ordermode

def main():
    vsize = 5
    args = get_args()
    mrc_file = args.mrc
    npz_file = args.npz
    apix = args.apix
    threshold = args.threshold
    ignorestart = args.ignorestart

    mrc = mrcfile.open(mrc_file, mode='r')

    data = mrc.data.copy()

    cella = mrc.header.cella.copy()
    [xlen, ylen, zlen] = cella.tolist()

    cellb = mrc.header.cellb.copy()
    [alpha, beta, gamma] = cellb.tolist()

    xpix = mrc.voxel_size.x
    ypix = mrc.voxel_size.y
    zpix = mrc.voxel_size.z

    nx = mrc.header.nx
    ny = mrc.header.ny
    nz = mrc.header.nz
    
    origin  = mrc.header.origin.copy()
    originx = origin.x
    originy = origin.y
    originz = origin.z

    nxstart = mrc.header.nxstart 
    nystart = mrc.header.nystart 
    nzstart = mrc.header.nzstart 

    mapc = mrc.header.mapc
    mapr = mrc.header.mapr
    maps = mrc.header.maps
    
    mrc.close()

    ordermode = getmode(mapc, mapr, maps)

    if ordermode == 1:
        xdim = nx
        ydim = ny
        zdim = nz
    elif ordermode == 2:
        xdim = nx
        ydim = nz
        zdim = ny
        data = data.transpose((1,0,2))
    elif ordermode == 3:
        xdim = ny
        ydim = nx
        zdim = nz
        data = data.transpose((0,2,1))
    elif ordermode == 4:
        xdim = nz
        ydim = nx
        zdim = ny
        data = data.transpose((1,2,0))
    elif ordermode == 5:
        xdim = ny
        ydim = nz
        zdim = nx
        data = data.transpose((2,0,1))
    elif ordermode == 6:
        xdim = nz
        ydim = ny
        zdim = nx
        data = data.transpose((2,1,0))

    nx = xdim
    ny = ydim
    nz = zdim
    
    if not ignorestart:
        if ordermode == 1: 
            originx += nxstart * xpix
            originy += nystart * ypix
            originz += nzstart * zpix
        elif ordermode == 2:
            originx += nxstart * xpix
            originy += nzstart * ypix
            originz += nystart * zpix
        elif ordermode == 3:
            originx += nystart * xpix
            originy += nxstart * ypix
            originz += nzstart * zpix
        elif ordermode == 4:
            originx += nzstart * xpix
            originy += nxstart * ypix
            originz += nystart * zpix
        elif ordermode == 5:
            originx += nystart * xpix
            originy += nzstart * ypix
            originz += nxstart * zpix
        elif ordermode == 6:
            originx += nzstart * xpix
            originy += nystart * ypix
            originz += nxstart * zpix
        origin.x = originx
        origin.y = originy
        origin.z = originz

    try:
        assert(alpha==beta==gamma==90.0)
    except AssertionError:
       print("[INFO] Input grid is not orthogonal.")
       sys.exit()

    try:
        assert(xpix==ypix==zpix==apix)
    except AssertionError:
        print("[INFO] Interpolating grid size from x:{:.2f} y:{:.2f} z:{:.2f} to a cubic grid of size {:.2f} using trilinear interpolation".format(xpix, ypix, zpix, apix))
        interp3d.interp3d.linear(data, zpix, ypix, xpix, apix, nz, ny, nx)
        data = interp3d.interp3d.mapout

        nz = interp3d.interp3d.pextz
        ny = interp3d.interp3d.pexty
        nx = interp3d.interp3d.pextx

    xlen = nx * apix
    ylen = ny * apix
    zlen = nz * apix

    cella.x = xlen
    cella.y = ylen
    cella.z = zlen

    nvoxels = 0
    for kx in range(0, nx):
        for ky in range(0, ny):
            for kz in range(0, nz):
                if data[kz, ky, kx] > threshold:
                    nvoxels += 1

    print("[INFO] {} voxels with density value greater than {:.5f} are retained.".format(nvoxels, threshold))

    voxels = np.zeros((nvoxels, 2 * vsize + 1, 2 * vsize + 1, 2 * vsize + 1), dtype = "float32") # zyx

    coords = np.zeros((nvoxels, 3), dtype = "int32") 
    nv = 0
    for kx in range(0, nx):
        for ky in range(0, ny):
            for kz in range(0, nz):
                if data[kz, ky, kx] > threshold:
                    x_lower = max(0, kx - vsize)
                    y_lower = max(0, ky - vsize)
                    z_lower = max(0, kz - vsize)

                    x_upper = min(nx - 1, kx + vsize)
                    y_upper = min(ny - 1, ky + vsize)
                    z_upper = min(nz - 1, kz + vsize)

                    dxl = x_lower - (kx - vsize)
                    dyl = y_lower - (ky - vsize)
                    dzl = z_lower - (kz - vsize)

                    dxu = (kx + vsize) - x_upper
                    dyu = (ky + vsize) - y_upper
                    dzu = (kz + vsize) - z_upper

                    coords[nv,:] = np.array([kz, ky, kx], dtype="int32")

                    voxels[                       nv,
                           dzl:(2 * vsize + 1 - dzu),
                           dyl:(2 * vsize + 1 - dyu),
                           dxl:(2 * vsize + 1 - dxu)] = data[z_lower:z_upper + 1, y_lower:y_upper + 1, x_lower:x_upper + 1]

                    vmax = np.max(voxels[nv])
                    vmin = np.min(voxels[nv])
                    voxels[nv] = (voxels[nv] - vmin) / (vmax - vmin) 

                    nv += 1
    del data

    print("[INFO] {} voxels are saved to {}".format(nv, npz_file))

    np.savez(npz_file, voxels = voxels[:nv, :],
                       coords = coords[:nv, :],
                                       nx = nx,
                                       ny = ny,
                                       nz = nz,
                                 cella = cella,
                               origin = origin,
                                   apix = apix,
            )
    del coords
    del voxels

if __name__ == "__main__":
    main()
