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

import torch
import mrcfile
import argparse
import numpy as np
from torch import nn
from tqdm import tqdm
from math import ceil
from torch import FloatTensor as FT
from torch.autograd import Variable as V

from DenseNetA import DenseNet

use_gpu = torch.cuda.is_available()
if use_gpu:
    torch.cuda.empty_cache()

def get_args():
    parser = argparse.ArgumentParser(description="This script predict main-chain probability and C-alpha probability for voxels in .npz format",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--npz", "-n", type=str, required=True,
                        help="Input .npz file")
    parser.add_argument("--mcmap", "-m", type=str, required=True,
                        help="Output predMC.mrc file")
    parser.add_argument("--camap", "-c", type=str, required=True,
                        help="Output predCA.mrc file")
    parser.add_argument("--bsize", "-b", type=int, default=128,
                        help="Number of points to predict in one batch")
    parser.add_argument("--model", "-M", type=str, default="../models/DenseNetA.pkl",
                        help="file name of trained model")
    args = parser.parse_args()
    return args

def main():
    vsize = 11
    args = get_args()
    npz_file = args.npz
    mcmap_file = args.mcmap
    camap_file = args.camap
    batch_size = args.bsize
    model_file = args.model

    print("[LOAD] loading data...")
    data = np.load(npz_file)
    coords = data['coords']

    X = data['voxels'].reshape(-1, vsize ** 3)

    nx = data['nx']
    ny = data['ny']
    nz = data['nz']

    cella = data['cella']
    origin = data['origin']
    apix = data['apix']

    origin_list = origin.tolist()
    [originx , originy , originz]= origin_list

    X = V(FT(X), requires_grad=False).view(-1, 1, vsize, vsize, vsize)

    data_num = len(X)

    n_batches = ceil(data_num/batch_size)

    model = DenseNet().cuda()
    model = nn.DataParallel(model)
    model.module.load_state_dict(torch.load(model_file))

    model.eval()

    CA_map = np.zeros((nz, ny, nx), dtype="float32")
    MC_map = np.zeros((nz, ny, nx), dtype="float32")

    print("[PREDICT] predicting...")
    for i in tqdm(range(n_batches)):
        coords_batch = coords[i * batch_size : (i + 1) * batch_size]
        coords_batch = tuple( map( tuple, coords_batch) )

        if use_gpu:
            X_batch = X[i * batch_size : (i + 1) * batch_size].cuda()
            y_pred = model(X_batch).cpu().detach().numpy()
        else:
            X_batch = X[i * batch_size : (i + 1) * batch_size]
            y_pred = model(X_batch).detach().numpy()
            

        for j in range( len( coords_batch ) ):
            MC_map[ coords_batch[j] ] = y_pred[j, 0]
            CA_map[ coords_batch[j] ] = y_pred[j, 1]

    new_mrc = mrcfile.new(mcmap_file, overwrite=True)
    new_mrc.voxel_size = apix
    new_mrc.header.origin = origin
    new_mrc.header.cella = cella
    new_mrc.set_data(MC_map)
    new_mrc.close()

    new_mrc = mrcfile.new(camap_file, overwrite=True)
    new_mrc.voxel_size = apix
    new_mrc.header.origin = origin
    new_mrc.header.cella = cella
    new_mrc.set_data(CA_map)
    new_mrc.close()

if __name__ == "__main__":
    main()
