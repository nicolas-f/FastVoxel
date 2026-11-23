import pytest
import fastvoxel as fv


@pytest.mark.parametrize(
    "box_min, box_max, resolution",
    [
        # Even x Odd x Odd after +-2*res padding -> (8, 9, 5)
        ((0.0, 0.0, 0.0), (4.0, 5.0, 1.0), 1.0),
        # All-even grid after +-2*res padding -> (6, 8, 10)
        ((-1.0, -2.0, -3.0), (1.0, 2.0, 3.0), 1.0),
        # All-odd grid after +-2*res padding -> (7, 9, 11)
        ((0.0, 0.0, 0.0), (0.6, 1.0, 1.4), 0.2),
    ])

def test_grid_origin_alignment_and_coord_lookup(
        box_min, box_max, resolution):
    """
    The first voxel center must preserve the historical grid origin,
    and voxel-center coordinates must map back to their cell IDs.

    Parameters
    ----------
    box_min : tuple
        The minimum corner of the bounding box.
    box_max : tuple
        The maximum corner of the bounding box.
    resolution : float
        The voxel resolution.

    Returns
    -------
    None
    """

    voxelizer = fv.TriangleScalarFieldCreator(resolution)

    bmin = fv.dvec3(*box_min)
    bmax = fv.dvec3(*box_max)

    voxelizer.FirstStep_Params(bmin, bmax)

    shape = voxelizer.get_domain_shape()
    nx, ny, nz = shape.x, shape.y, shape.z

    # First voxel center corresponds to the lower padded bounding-box limit.
    expected_first_center = tuple(value - 2.0 * resolution for value in box_min)

    c_first = voxelizer.get_center_cell_coordinates(fv.ivec3(0, 0, 0))

    assert (c_first.x, c_first.y, c_first.z) == pytest.approx(
        expected = expected_first_center,
        abs = resolution * 1e-6)

    # Last voxel center follows from the grid dimensions.
    c_last = voxelizer.get_center_cell_coordinates(fv.ivec3(nx - 1, ny - 1, nz - 1))

    expected_last_center = tuple(expected_first_center[i] + (shape[i] - 1) * resolution
                                 for i in range(3))

    assert (c_last.x, c_last.y, c_last.z) == pytest.approx(
        expected = expected_last_center,
        abs = resolution * 1e-6)

    # Coordinates strictly inside a voxel must map to that voxel.
    sample_cells = [fv.ivec3(0, 0, 0),
                    fv.ivec3(nx // 2, ny // 2, nz // 2),
                    fv.ivec3(nx - 1, ny - 1, nz - 1)]

    for cell in sample_cells:
        center = voxelizer.get_center_cell_coordinates(cell_id = cell)

        for dx in (-0.25, 0.0, 0.25):
            for dy in (-0.25, 0.0, 0.25):
                for dz in (-0.25, 0.0, 0.25):
                    pt = fv.dvec3(center.x + dx * resolution,
                                  center.y + dy * resolution,
                                  center.z + dz * resolution)

                    mapped = voxelizer.get_cell_id_by_coord(position = pt)

                    assert (mapped.x, mapped.y, mapped.z) == (
                        cell.x, cell.y, cell.z)

def test_triangle_voxelization_matches_coord_lookup_on_even_axis():
    """Rasterizing a triangle inside voxel (0, 0, 0) must agree with get_cell_id_by_coord."""
    resolution = 1.0
    voxelizer = fv.TriangleScalarFieldCreator(resolution)
    # Padded size = (8, 9, 5): X is even (8), Y and Z are odd.
    voxelizer.FirstStep_Params(fv.dvec3(0.0, 0.0, 0.0), fv.dvec3(4.0, 5.0, 1.0))

    c0 = voxelizer.get_center_cell_coordinates(fv.ivec3(0, 0, 0))
    # Small non-degenerate triangle strictly inside voxel (0, 0, 0)
    a = fv.dvec3(c0.x - 0.2, c0.y - 0.2, c0.z)
    b = fv.dvec3(c0.x + 0.2, c0.y - 0.2, c0.z)
    c = fv.dvec3(c0.x, c0.y + 0.2, c0.z)
    voxelizer.SecondStep_PushTri(a, b, c, 1)

    # Point inside the triangle's voxel on the negative side of center
    query = fv.dvec3(c0.x - 0.25, c0.y, c0.z)
    cell = voxelizer.get_cell_id_by_coord(query)
    assert (cell.x, cell.y, cell.z) == (0, 0, 0)
    assert voxelizer.get_matrix_value(cell) == 1
