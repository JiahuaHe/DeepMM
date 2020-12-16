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
import pandas as pd
from torch import nn
from tqdm import tqdm
from math import ceil
from torch import FloatTensor as FT
from torch.autograd import Variable as V

from DenseNetB import DenseNet

np.set_printoptions(threshold = np.inf) 
np.set_printoptions(suppress = True)

use_gpu = torch.cuda.is_available()
if use_gpu:
    torch.cuda.empty_cache()

def get_args():
    parser = argparse.ArgumentParser(description="This script predict amino acid and secondary structure type for LDPs",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mcv", "-m", type=str, required=True,
                        help="input voxel density file")
    parser.add_argument("--ldp", "-l", type=str, required=True,
                        help="input LDP file in .pdb format")
    parser.add_argument("--bsize", "-b", type=int, default=128,
                        help="number of data points to predict in one batch")
    parser.add_argument("--model", "-M", type=str, default="../models/DenseNetB.pkl",
                        help="file name of trained model")
    args = parser.parse_args()
    return args

def read_file(file_name):
	f = open(file_name,'r')
	source = []
	for line in f.readlines():
		source.append(line.strip())
	f.close()
	return source

def main():
    vsize = 10

    args = get_args()
    mcv_file = args.mcv
    ldp_file = args.ldp
    batch_size = args.bsize
    model_file = args.model

    ldp=read_file(ldp_file)    

    df = pd.read_csv(mcv_file, sep=' ', header=None)
    data = np.array(df.iloc[:,1:], dtype='float32')

    X = data.reshape(-1, vsize ** 3)

    X = V(FT(X), requires_grad=False).view(-1, 1, vsize, vsize, vsize)

    data_num = len(X)

    n_batches = ceil(data_num/batch_size)

    model = DenseNet().cuda()
    model = nn.DataParallel(model)
    model.module.load_state_dict(torch.load(model_file))

    model.eval()
    y = np.zeros((data_num,  7), dtype='float32')
    for i in tqdm(range(n_batches)):
        X_batch = X[i * batch_size : (i + 1) * batch_size].cuda()
        y_batch = model(X_batch).cpu().detach().numpy()
        y[i * batch_size : (i + 1) * batch_size] = y_batch
   
    i = 0
    for line in ldp:
        if line[0:4] == "ATOM":
            print("{}{:6d}{:6d}".format(line, np.argmax(y[i,:4]), np.argmax(y[i,4:])) )
            i += 1
        else:
            print(line)
            
if __name__ == "__main__":
    main()
