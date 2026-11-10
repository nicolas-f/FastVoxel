# -*- coding: cp1252 -*-
import numpy as np
import fastvoxel as fv
import time

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
pieceof = np.zeros((extract_shape[0], extract_shape[1], extract_shape[2]), dtype=np.short)
voxelizator.copy_matrix(pieceof, minv)
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
print("Copy matrix ok ?", (pieceof[:, :, pieceof.shape[2] // 2] == expected_res).all())

# import matplotlib.pyplot as plt
# from pylab import show
# plt.matshow(pieceof[:,:,0],fignum=2)
# show()
