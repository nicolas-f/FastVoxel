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

def test_grid_center_symmetry_and_coord_lookup(box_min, box_max, resolution):
    """Grid center must be symmetrically centered around the bounding-box center
    and any point strictly inside a voxel (center +/- 0.25 * resolution) must be
    mapped to the same voxel."""
    voxelizer = fv.TriangleScalarFieldCreator(resolution)
    bmin = fv.dvec3(*box_min)
    bmax = fv.dvec3(*box_max)
    voxelizer.FirstStep_Params(bmin, bmax)

    shape = voxelizer.get_domain_shape()
    nx, ny, nz = shape.x, shape.y, shape.z

    # 1. Grid must be symmetrically centered around the bounding-box center
    expected_center = (0.5 * (box_min[0] + box_max[0]),
                       0.5 * (box_min[1] + box_max[1]),
                       0.5 * (box_min[2] + box_max[2]))
    c_first = voxelizer.get_center_cell_coordinates(fv.ivec3(0, 0, 0))
    c_last = voxelizer.get_center_cell_coordinates(fv.ivec3(nx - 1, ny - 1, nz - 1))
    actual_center = (0.5 * (c_first.x + c_last.x),
                     0.5 * (c_first.y + c_last.y),
                     0.5 * (c_first.z + c_last.z))
    assert actual_center == pytest.approx(expected_center, abs=1e-12)

    # 2. Any point strictly inside a voxel (center +/- 0.25 * resolution)
    #    must map back to that exact voxel id via get_cell_id_by_coord.
    sample_cells = [fv.ivec3(0, 0, 0),
                    fv.ivec3(nx // 2, ny // 2, nz // 2),
                    fv.ivec3(nx - 1, ny - 1, nz - 1)]
    for cell in sample_cells:
        center = voxelizer.get_center_cell_coordinates(cell)
        for dx in (-0.25, 0.0, 0.25):
            for dy in (-0.25, 0.0, 0.25):
                for dz in (-0.25, 0.0, 0.25):
                    pt = fv.dvec3(center.x + dx * resolution,
                                  center.y + dy * resolution,
                                  center.z + dz * resolution)
                    mapped = voxelizer.get_cell_id_by_coord(pt)
                    assert (mapped.x, mapped.y, mapped.z) == (cell.x, cell.y, cell.z)

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
