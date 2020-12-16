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

import interp3d

def get_args():
    parser = argparse.ArgumentParser(description="This script transfers a .mrc map file to a .situs map file with given pixel spacing", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mrc", "-m", type=str, required=True,
                        help="Input .mrc file")
    parser.add_argument("--situs", "-s", type=str, required=True,
                        help="Output .situs file")
    parser.add_argument("--apix", "-a", type=float, default=1.0,
                        help="Pixel spacing, default = 1.0)")
    parser.add_argument("--ignorestart", "-i", type=bool, default=False,
                        help="Whether ignore start position or not, 'True' for simulated maps and 'False' for experimental maps, default =  False")
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
        print("[EXIT] Input map file gives malformed dimension ordering.")
        sys.exit()
    return ordermode

def main():
    args = get_args()
    apix = args.apix
    mrc_file = args.mrc
    situs_file = open(args.situs, 'w')
    ignorestart = args.ignorestart

    mrc = mrcfile.open(mrc_file, mode='r')
    data = mrc.data.copy()

    cella = mrc.header.cella.copy()
    [xlen, ylen, zlen] = cella.tolist()

    xpix = mrc.voxel_size.x
    ypix = mrc.voxel_size.y
    zpix = mrc.voxel_size.z

    cellb = mrc.header.cellb.copy()
    [alpha, beta, gamma] = cellb.tolist()

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
        print("[EXIT] grid is not orthogonal.")
        sys.exit()

    try:
        assert(xpix==ypix==zpix==apix)
    except AssertionError:
        print("[INFO] interpolating grid size from x:{:.2f} y:{:.2f} z:{:.2f} to a cubic grid of size {:.2f} using trilinear interpolation".format(xpix, ypix, zpix, apix))
        interp3d.interp3d.linear(data, zpix, ypix, xpix, apix, nz, ny, nx)
        data = interp3d.interp3d.mapout

        nz = interp3d.interp3d.pextz
        ny = interp3d.interp3d.pexty
        nx = interp3d.interp3d.pextx

    xlen = nx * apix
    ylen = ny * apix
    zlen = nz * apix

    nn = 0
    print("{:8.6f} {:8.6f} {:8.6f} {:8.6f} {} {} {}".format(apix, originx, originy, originz, nx, ny, nz), file=situs_file)
    print("", file=situs_file)
    for kz in range(nz):
        for ky in range(ny):
            for kx in range(nx):
                nn += 1
                print("{:11.6f}".format(data[kz, ky, kx]), end=" ", file=situs_file)
                if nn == 10:
                    print("", file=situs_file)
                    nn = 0

if __name__ == "__main__":
    main()
