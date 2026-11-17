import fastvoxel as fv
import numpy as np
import time
from pathlib import Path


def run():
    resolution = 0.1
    doplot = True
    doexportvtk = False

    if doexportvtk:
        resolution = 0.5
    deb = time.time()
    # Init Voxelisation with precision in meter
    voxelizator = fv.TriangleScalarFieldCreator(resolution)
    # Polygonal file path
    ply_file = Path(__file__).resolve().with_name("elmia.ply")
    # Load PLY polygonal file
    if voxelizator.load_ply_model(str(ply_file)):
        print("Process done in %f sec" % (time.time() - deb))
        # Extract the cell i,j,k values corresponding to position in meter (a position in the volume to extract)
        cellid = voxelizator.get_cell_id_by_coord(fv.dvec3(7, -0.5, 1))
        # Get the corresponding volume id
        elmiavol = voxelizator.get_matrix_value(cellid)
        # Extract Boundary of this volume
        minv = fv.ivec3()
        maxv = fv.ivec3()
        voxelizator.get_cell_value_boundaries(minv, maxv, elmiavol)
        minv -= fv.ivec3(1, 1, 1)
        maxv += fv.ivec3(1, 1, 1)
        extract_shape = maxv - minv

        print("Matrix Size is %ix%ix%i" % (extract_shape[0], extract_shape[1], extract_shape[2]))

        # Plotting
        if doplot:
            print("voxelizator.get_first_volume_index():", voxelizator.get_first_volume_index())
            nbmarkers = voxelizator.get_first_volume_index() + voxelizator.get_volume_count()
            # Create a matrix to extract only an Y slice of the 3D voxelisation (save lot of memory)
            pieceof = np.zeros((extract_shape[0], 1, extract_shape[2]), dtype=np.short)
            # Transfer data from voxel to numpy matrix
            data_filt = np.arange(nbmarkers, dtype=np.short)
            data_filt[voxelizator.get_first_volume_index():] = -1  # Set -1 to all volumes
            data_filt[elmiavol] = 0  # Set 0 to air
            voxelizator.copy_matrix_filtered(pieceof, minv + fv.ivec3(0, cellid[1], 0), data_filt)
            import matplotlib.pyplot as plt
            from pylab import colorbar
            p = plt.matshow(np.rot90(pieceof[:, 0, :]), fignum=2, cmap=plt.get_cmap("jet", nbmarkers))
            colorbar(p, ticks=range(nbmarkers))
            output_file = Path(__file__).resolve().with_name("elmia_slice.png")
            plt.savefig(output_file, dpi=200, bbox_inches="tight")
            plt.close()
            print("Voxelized elmia slice saved in {}".format(output_file))
        elif doexportvtk:
            from evtk.hl import imageToVTK
            # Create a matrix to extract only an Y slice of the 3D voxelisation (save lot of memory)
            pieceof = np.zeros((extract_shape[0], extract_shape[1], extract_shape[2]), dtype=np.short)
            # Transfer data from voxel to numpy matrix
            voxelizator.copy_matrix(pieceof, minv)
            imageToVTK(r"C:\tmp\elmia", cellData={"materials": np.array(pieceof, dtype="float32")})
            time.sleep(1)


run()
