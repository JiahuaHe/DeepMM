'''
MIT License

Copyright (c) 2018 Geoff Pleiss
Copyright (c) 2020 Jiahua He

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.checkpoint as cp
from collections import OrderedDict
from torch.autograd import Variable as V

torch.manual_seed(2 ** 10)

def _bn_function_factory(norm, relu, conv):
    def bn_function(*inputs):
        concated_features = torch.cat(inputs, 1)
        bottleneck_output = conv(relu(norm(concated_features)))
        return bottleneck_output

    return bn_function


class _DenseLayer(nn.Module):
    def __init__(self, num_input_features, growth_rate, bn_size, drop_rate, efficient=False):
        super(_DenseLayer, self).__init__()
        self.add_module('norm1', nn.BatchNorm3d(num_input_features)),
        self.add_module('relu1', nn.ReLU(inplace=True)),
        self.add_module('conv1', nn.Conv3d(num_input_features, bn_size * growth_rate,
                        kernel_size=1, stride=1, bias=False)),
        self.add_module('norm2', nn.BatchNorm3d(bn_size * growth_rate)),
        self.add_module('relu2', nn.ReLU(inplace=True)),
        self.add_module('conv2', nn.Conv3d(bn_size * growth_rate, growth_rate,
                        kernel_size=3, stride=1, padding=1, bias=False)),
        self.drop_rate = drop_rate
        self.efficient = efficient

    def forward(self, *prev_features):
        bn_function = _bn_function_factory(self.norm1, self.relu1, self.conv1)
        if self.efficient and any(prev_feature.requires_grad for prev_feature in prev_features):
            bottleneck_output = cp.checkpoint(bn_function, *prev_features)
        else:
            bottleneck_output = bn_function(*prev_features)
        new_features = self.conv2(self.relu2(self.norm2(bottleneck_output)))
        if self.drop_rate > 0:
            new_features = F.dropout(new_features, p=self.drop_rate, training=self.training)
        return new_features

class _Transition(nn.Sequential):
    def __init__(self, num_input_features, num_output_features):
        super(_Transition, self).__init__()
        self.add_module('norm', nn.BatchNorm3d(num_input_features))
        self.add_module('relu', nn.ReLU(inplace=True))
        self.add_module('conv', nn.Conv3d(num_input_features, num_output_features,
                                          kernel_size=1, stride=1, bias=False))
        self.add_module('pool', nn.AvgPool3d(kernel_size=2, stride=2))

class _Transition_no_pool(nn.Sequential):
    def __init__(self, num_input_features, num_output_features):
        super(_Transition_no_pool, self).__init__()
        self.add_module('norm', nn.BatchNorm3d(num_input_features))
        self.add_module('relu', nn.ReLU(inplace=True))
        self.add_module('conv', nn.Conv3d(num_input_features, num_output_features,
                                          kernel_size=1, stride=1, bias=False))

class _DenseBlock(nn.Module):
    def __init__(self, num_layers, num_input_features, bn_size, growth_rate, drop_rate, efficient=False):
        super(_DenseBlock, self).__init__()
        for i in range(num_layers):
            layer = _DenseLayer(
                num_input_features + i * growth_rate,
                growth_rate=growth_rate,
                bn_size=bn_size,
                drop_rate=drop_rate,
                efficient=efficient,
            )
            self.add_module('denselayer_%d' % (i + 1), layer)

    def forward(self, init_features):
        features = [init_features]
        for name, layer in self.named_children():
            new_features = layer(*features)
            features.append(new_features)
        return torch.cat(features, 1)


class DenseNet(nn.Module):
    def __init__(self, growth_rate=12, shared_block_config=(8, 8), MC_block_config=(8, 8), CA_block_config=(8, 8), compression=0.5,
                 num_init_features=32, bn_size=4, drop_rate=0.2,
                 num_classes_MC=1, num_classes_CA=1, efficient=False):

        super(DenseNet, self).__init__()
        assert 0 < compression <= 1, 'compression of densenet should be between 0 and 1'

        # First convolution
        self.shared_blocks = nn.Sequential(OrderedDict([
            ('conv0', nn.Conv3d(1, num_init_features, kernel_size=5, stride=1, padding=2, bias=False)),
        ]))

        # Shared blocks
        num_features = num_init_features
        for i, num_layers in enumerate(shared_block_config):
            block = _DenseBlock(
                num_layers=num_layers,
                num_input_features=num_features,
                bn_size=bn_size,
                growth_rate=growth_rate,
                drop_rate=drop_rate,
                efficient=efficient,
            )
            self.shared_blocks.add_module('denseblock%d' % (i + 1), block)
            num_features = num_features + num_layers * growth_rate
#            if i != len(shared_block_config) - 1:
            if True:
                trans = _Transition_no_pool(num_input_features=num_features,
                                    num_output_features=int(num_features * compression))
                self.shared_blocks.add_module('transition%d' % (i + 1), trans)
                num_features = int(num_features * compression)

        num_shared_features = num_features

        self.MC_blocks=nn.Sequential()

        num_features=num_shared_features
        for i, num_layers in enumerate(MC_block_config):
            block = _DenseBlock(
                num_layers=num_layers,
                num_input_features=num_features,
                bn_size=bn_size,
                growth_rate=growth_rate,
                drop_rate=drop_rate,
                efficient=efficient,
            )
            self.MC_blocks.add_module('denseblock_%d' % (i + 1), block)
            num_features = num_features + num_layers * growth_rate
            if i != len(CA_block_config) - 1:
                trans = _Transition(num_input_features=num_features,
                                    num_output_features=int(num_features * compression))
                self.MC_blocks.add_module('transition_%d' % (i + 1), trans)
                num_features = int(num_features * compression)
        # Final batch norm
        self.MC_blocks.add_module('norm_final', nn.BatchNorm3d(num_features))
        # Final ReLU
        self.MC_blocks.add_module('relu_final', nn.ReLU(inplace=True))
        # Final pool
        self.MC_blocks.add_module('pool_final', nn.AdaptiveAvgPool3d(1))

        num_features_MC = num_features
         
        self.MC_classifier = nn.Linear(num_features_MC, num_classes_MC)

        self.CA_blocks=nn.Sequential()

        num_features=num_shared_features
        for i, num_layers in enumerate(CA_block_config):
            block = _DenseBlock(
                num_layers=num_layers,
                num_input_features=num_features,
                bn_size=bn_size,
                growth_rate=growth_rate,
                drop_rate=drop_rate,
                efficient=efficient,
            )
            self.CA_blocks.add_module('denseblock_%d' % (i + 1), block)
            num_features = num_features + num_layers * growth_rate
            if i != len(CA_block_config) - 1:
                trans = _Transition(num_input_features=num_features,
                                    num_output_features=int(num_features * compression))
                self.CA_blocks.add_module('transition_%d' % (i + 1), trans)
                num_features = int(num_features * compression)

        # Final batch norm
        self.CA_blocks.add_module('norm_final', nn.BatchNorm3d(num_features))

        # Final ReLU
        self.CA_blocks.add_module('relu_final', nn.ReLU(inplace=True))

        # Final pool
        self.CA_blocks.add_module('pool_final', nn.AdaptiveAvgPool3d(1))

        num_features_CA = num_features

        self.CA_classifier = nn.Linear(num_features_CA, num_classes_CA)

    def forward(self, x):
        out = self.shared_blocks(x)

        mc = self.MC_blocks(out)
        mc = self.MC_classifier(mc.view(mc.size(0), -1))

        ca = self.CA_blocks(out)
        ca = self.CA_classifier(ca.view(ca.size(0), -1))

        out = torch.cat((mc,ca), axis=1)
        return out
