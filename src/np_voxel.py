# -*- coding: cp1252 -*-
import numpy as np
import math
import os
from fastvoxel import ivec3
##
# Use multiple memmap to handle a huge matrix
def as_long_array(lst):
    return np.asarray(lst,dtype=np.long)
class np_voxel(object):
    def __init__(self,fastvoxel,shape=None,range_beg=None):
        self._fastvoxel=fastvoxel
        self.dtype=np.short
        if shape is None:
            shape=(fastvoxel.get_domain_size(),fastvoxel.get_domain_size(),fastvoxel.get_domain_size())
        self.shape=shape
        if range_beg is None:
            range_beg=(0,0,0)
        self._vrange_beg=ivec3(range_beg[0],range_beg[1],range_beg[2])
        self._filter=None
    def set_filter(self,filt_array):
        self._filter=filt_array
    def __getitem__(self, key):
        if len(key)!=3:
            raise(NotImplementedError("Only 3 dimensional matrix get value supported"))
        get_origin_pos=[0,0,0]
        get_destination_pos=[0,0,0]
        single_dim=[slice(None),slice(None),slice(None)]
        has_single=False
        for dim in range(3):
            if type(key[dim]) is slice:
                #Range
                range_get=key[dim].indices(self.shape[dim])
                get_origin_pos[dim]=range_get[0]
                get_destination_pos[dim]=range_get[1]
            else:
                #get only one item
                has_single=True 
                single_dim[dim]=0
                get_origin_pos[dim]=key[dim]
                get_destination_pos[dim]=key[dim]+1
        ret_shape=tuple(list(as_long_array(get_destination_pos)-as_long_array(get_origin_pos)))
        ret_value=np.empty(shape=ret_shape,dtype=self.dtype)
        vget_origin_pos=ivec3(get_origin_pos[0],get_origin_pos[1],get_origin_pos[2])
        if self._filter is None:
            self._fastvoxel.copy_matrix(ret_value,vget_origin_pos+self._vrange_beg)
        else:
            self._fastvoxel.copy_matrix_filtered(ret_value,vget_origin_pos+self._vrange_beg,self._filter)
        if not has_single:
            return ret_value
        else:
            return ret_value[tuple(single_dim)]
