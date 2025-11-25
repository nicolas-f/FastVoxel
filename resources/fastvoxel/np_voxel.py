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

if __name__ == '__main__':
    import fastvoxel as fv

    boxmin = fv.dvec3(0, 0, 0)
    boxmax = fv.dvec3(5, 5, 5)

    voxelizator = fv.TriangleScalarFieldCreator(0.5)

    voxelizator.first_step_params(boxmin, boxmax)
    # Définition des 8 sommets du cube           
    sommets = [fv.dvec3(5.0, 0.0, 0.0),
               fv.dvec3(0.0, 0.0, 0.0),
               fv.dvec3(0.0, 5.0, 0.0),
               fv.dvec3(5.0, 5.0, 0.0),
               fv.dvec3(0.0, 5.0, 5.0),
               fv.dvec3(5.0, 5.0, 5.0),
               fv.dvec3(0.0, 0.0, 5.0),
               fv.dvec3(5.0, 0.0, 5.0)
               ]
    # Définition des 12 faces triangulaire du cube
    # [ sommetA, sommetB, sommetC, idencombrement, idmateriau, idrecepteursurf ]
    faces = [[0, 1, 2, -1, 66, -1],
             [0, 2, 3, -1, 66, -1],
             [2, 4, 5, -1, 100, -1],
             [2, 5, 3, -1, 100, -1],
             [2, 6, 4, -1, 100, -1],
             [2, 1, 6, -1, 100, -1],
             [1, 0, 7, -1, 100, -1],
             [6, 1, 7, -1, 100, -1],
             [0, 3, 5, -1, 100, -1],
             [7, 0, 5, -1, 100, -1],
             [7, 5, 4, -1, 66, -1],
             [6, 7, 4, -1, 66, -1]
             ]
    for facedata in faces:
        voxelizator.second_step_pushtri(sommets[facedata[0]],
                                        sommets[facedata[1]],
                                        sommets[facedata[2]], facedata[4])
    voxelizator.third_step_volumescreator()
    cellid = voxelizator.get_cell_id_by_coord(fv.dvec3(2.5, 2.5, 2.5))
    cubevol = voxelizator.get_matrix_value(cellid)
    minv = fv.ivec3()
    maxv = fv.ivec3()
    voxelizator.get_cell_value_boundaries(minv, maxv, cubevol)
    minv -= fv.ivec3(1, 1, 1)
    maxv += fv.ivec3(1, 1, 1)
    extract_shape = maxv - minv
    vox_array = np_voxel(voxelizator, shape=(extract_shape[0], extract_shape[1], extract_shape[2]),
                         range_beg=(minv[0], minv[1], minv[2]))
    expected_res = np.asarray([[100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
                               [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
                               [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
                               [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
                               [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
                               [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
                               [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
                               [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
                               [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
                               [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
                               [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]], dtype=np.short)
    pieceof = vox_array[:, :, :]
    print(pieceof[:, :, extract_shape[2] / 2] == expected_res).all()
