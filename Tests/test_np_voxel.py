import unittest
import numpy as np
import fastvoxel as fv
from fastvoxel.np_voxel import np_voxel


class TestNpVoxel(unittest.TestCase):
    """Unit tests for the np_voxel function of FastVoxel"""

    def setUp(self):
        """Initialize common data for all tests"""
        self.boxmin = fv.dvec3(0, 0, 0)
        self.boxmax = fv.dvec3(5, 5, 5)
        self.voxel_size = 0.5

        # Definition of the 8 vertices of the cube
        self.sommets = [
            fv.dvec3(5.0, 0.0, 0.0),
            fv.dvec3(0.0, 0.0, 0.0),
            fv.dvec3(0.0, 5.0, 0.0),
            fv.dvec3(5.0, 5.0, 0.0),
            fv.dvec3(0.0, 5.0, 5.0),
            fv.dvec3(5.0, 5.0, 5.0),
            fv.dvec3(0.0, 0.0, 5.0),
            fv.dvec3(5.0, 0.0, 5.0)
        ]

        # Definition of the 12 triangular faces of the cube
        # [vertexA, vertexB, vertexC, idencombrement, idmaterial, idrecepteursurf]
        self.faces = [
            [0, 1, 2, -1, 66, -1],
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

        # Expected result for a slice in the middle of the cube
        self.expected_res = np.asarray([
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
            [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
            [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
            [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
            [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
            [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
            [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
            [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
            [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
            [100, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100],
            [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
        ], dtype=np.short)

    def _create_voxelizator(self):
        """Creates and configures a voxelizator with the test cube"""
        voxelizator = fv.TriangleScalarFieldCreator(self.voxel_size)
        voxelizator.first_step_params(self.boxmin, self.boxmax)

        # Add the cube faces
        for facedata in self.faces:
            voxelizator.second_step_pushtri(
                self.sommets[facedata[0]],
                self.sommets[facedata[1]],
                self.sommets[facedata[2]],
                facedata[4]
            )

        voxelizator.third_step_volumescreator()
        return voxelizator

    def test_cube_middle_slice(self):
        """Test the slice in the middle of the voxelized cube"""
        voxelizator = self._create_voxelizator()

        # Get the cell at the center of the cube
        cellid = voxelizator.get_cell_id_by_coord(fv.dvec3(2.5, 2.5, 2.5))
        cubevol = voxelizator.get_matrix_value(cellid)

        # Calculate the boundaries of the region to extract
        minv = fv.ivec3()
        maxv = fv.ivec3()
        voxelizator.get_cell_value_boundaries(minv, maxv, cubevol)
        minv -= fv.ivec3(1, 1, 1)
        maxv += fv.ivec3(1, 1, 1)
        extract_shape = maxv - minv

        # Extract the voxel array
        vox_array = np_voxel(
            voxelizator,
            shape=(extract_shape[0], extract_shape[1], extract_shape[2]),
            range_beg=(minv[0], minv[1], minv[2])
        )

        # Check the middle slice
        middle_slice = vox_array[:, :, extract_shape[2] // 2]

        self.assertTrue(
            np.array_equal(middle_slice, self.expected_res),
            "The middle slice of the cube does not match the expected result"
        )

    def test_voxel_array_shape(self):
        """Test that the voxel array has the expected shape"""
        voxelizator = self._create_voxelizator()

        cellid = voxelizator.get_cell_id_by_coord(fv.dvec3(2.5, 2.5, 2.5))
        cubevol = voxelizator.get_matrix_value(cellid)

        minv = fv.ivec3()
        maxv = fv.ivec3()
        voxelizator.get_cell_value_boundaries(minv, maxv, cubevol)
        minv -= fv.ivec3(1, 1, 1)
        maxv += fv.ivec3(1, 1, 1)
        extract_shape = maxv - minv

        vox_array = np_voxel(
            voxelizator,
            shape=(extract_shape[0], extract_shape[1], extract_shape[2]),
            range_beg=(minv[0], minv[1], minv[2])
        )

        expected_shape = (extract_shape[0], extract_shape[1], extract_shape[2])
        self.assertEqual(
            vox_array.shape,
            expected_shape,
            f"La forme du tableau devrait être {expected_shape}"
        )

    def test_voxel_array_dtype(self):
        """Test that the voxel array has the correct data type"""
        voxelizator = self._create_voxelizator()

        cellid = voxelizator.get_cell_id_by_coord(fv.dvec3(2.5, 2.5, 2.5))
        cubevol = voxelizator.get_matrix_value(cellid)

        minv = fv.ivec3()
        maxv = fv.ivec3()
        voxelizator.get_cell_value_boundaries(minv, maxv, cubevol)
        minv -= fv.ivec3(1, 1, 1)
        maxv += fv.ivec3(1, 1, 1)
        extract_shape = maxv - minv

        vox_array = np_voxel(
            voxelizator,
            shape=(extract_shape[0], extract_shape[1], extract_shape[2]),
            range_beg=(minv[0], minv[1], minv[2])
        )

        self.assertEqual(
            vox_array.dtype,
            np.short,
            "The data type should be np.short"
        )

    def test_interior_exterior_values(self):
        """Test that the interior and exterior values are correct"""
        voxelizator = self._create_voxelizator()

        cellid = voxelizator.get_cell_id_by_coord(fv.dvec3(2.5, 2.5, 2.5))
        cubevol = voxelizator.get_matrix_value(cellid)

        minv = fv.ivec3()
        maxv = fv.ivec3()
        voxelizator.get_cell_value_boundaries(minv, maxv, cubevol)
        minv -= fv.ivec3(1, 1, 1)
        maxv += fv.ivec3(1, 1, 1)
        extract_shape = maxv - minv

        vox_array = np_voxel(
            voxelizator,
            shape=(extract_shape[0], extract_shape[1], extract_shape[2]),
            range_beg=(minv[0], minv[1], minv[2])
        )

        middle_slice = vox_array[:, :, extract_shape[2] // 2]

        # Check border values (100)
        self.assertTrue(
            np.all(middle_slice[0, :] == 100),
            "The top border values should be 100"
        )
        self.assertTrue(
            np.all(middle_slice[-1, :] == 100),
            "The bottom border values should be 100"
        )
        self.assertTrue(
            np.all(middle_slice[:, 0] == 100),
            "The left border values should be 100"
        )
        self.assertTrue(
            np.all(middle_slice[:, -1] == 100),
            "The right border values should be 100"
        )

        # Check interior values (102)
        interior = middle_slice[1:-1, 1:-1]
        self.assertTrue(
            np.all(interior == 102),
            "All interior values should be 102"
        )


if __name__ == '__main__':
    unittest.main()
