# FastVoxel

FastVoxel converts a 3D triangular mesh stored in a PLY file into a 3D NumPy-compatible voxel representation. Each voxel is associated with a material identifier, and enclosed volumes can be identified from the voxelized geometry.

FastVoxel does not require the input mesh to be manifold.

## Features

- Voxelization of triangular meshes from PLY files.
- Material identifiers associated with the voxelized geometry.
- Extraction of enclosed volumes.
- Support for **rectangular voxel grids**, allowing different dimensions along the X, Y, and Z axes.
- Double-precision coordinates and voxel resolution through the `dvec3` interface.
- NumPy-compatible access to the voxel matrix, including rectangular domains and filtered matrix extraction.
- Python exceptions for invalid arguments, out-of-range indices, and other runtime errors.

The voxel domain dimensions can be retrieved using:

```python
shape = voxelizer.get_domain_shape()
```

This method returns the domain shape as `(X, Y, Z)`. The historical `get_domain_size()` method is also available.

## Requirements

- Python 3.8 or later
- A C++ compiler with C++11 support
- CMake
- SWIG 3 or later
- Python development headers and library
- NumPy

Matplotlib is additionally required to run the PLY visualization demo.

## Building and installation

FastVoxel uses CMake and SWIG to compile its C++ extension and generate the Python wrapper.

First, create and activate a Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install NumPy:

```bash
python -m pip install numpy
```

Store the project directory:

```bash
PROJECT_DIR="$PWD"
```

Create the build directory :

```bash
mkdir -p "$PROJECT_DIR/build"
cd "$PROJECT_DIR/build"
```

Make sure that CMake, SWIG, and a compatible C++ compiler are available on your system.

Configure the project:

```bash
cmake "$PROJECT_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DPYTHON_EXECUTABLE="$(command -v python)" \
    -DCMAKE_INSTALL_PREFIX="$PROJECT_DIR/install"
```

Compile the extension:

```bash
cmake --build . -- -j4
```

Install FastVoxel Python package:

```bash
cmake --build . --target install
```

The installed package is located in the `install/fastvoxel` directory.

Return to the project directory

```bash
cd "$PROJECT_DIR"
```

To use this installation, add the installation prefix to `PYTHONPATH`:

```bash
export PYTHONPATH="$PWD/install:$PYTHONPATH"
```

You can then import FastVoxel from Python:

```python
import fastvoxel as fv
```

## Python interface

The main Python classes and vector types are:

- `fastvoxel.dvec3`: three-dimensional vector using double-precision coordinates,
- `fastvoxel.ivec3`: three-dimensional integer vector,
- `fastvoxel.ScalarFieldCreator`: scalar field creation interface,
- `fastvoxel.TriangleScalarFieldCreator`: scalar field creation from triangular meshes.

The historical `fastvoxel.vec3` name is retained as an alias for `fastvoxel.dvec3`.

The NumPy interface is available through:

```python
from fastvoxel.np_voxel import np_voxel
```

For example, a PLY mesh can be voxelized and converted to a NumPy-compatible representation as follows:

```python
import fastvoxel as fv
from fastvoxel.np_voxel import np_voxel

voxelizer = fv.TriangleScalarFieldCreator(0.1)

if voxelizer.load_ply_model("model.ply"):
    voxels = np_voxel(voxelizer)
    print("Voxel grid shape:", voxels.shape)
```

The voxel grid shape is derived from the actual domain dimensions, so rectangular grids are supported.

## Demos

The `demo` directory contains two Python scripts illustrating the FastVoxel API.

### Triangle feeding

`demo/test_triangle_feed.py` demonstrates how to construct a voxelized volume from triangular faces without loading a PLY file.

The script:
1. Defines a cubic domain.
2. Feeds the 12 triangular faces of a cube to the voxelizer.
3. Generates the enclosed volumes.
4. Extracts the voxel matrix.
5. Compares the resulting material identifiers against an expected matrix.

Run the example:

```bash
PYTHONPATH="$PWD/install:$PYTHONPATH" python demo/test_triangle_feed.py
```

Expected output:

```bash
Copy matrix ok ? True
```

### PLY mesh voxelization

`demo/import_ply.py` demonstrates how to load and voxelize a PLY mesh.

The script:
1. Loads `demo/elmia.ply`.
2. Voxelizes the mesh using a resolution of 0.1 m.
3. Extracts a slice of the voxel matrix.
4. Saves the resulting visualization as `demo/elmia_slice.png`.

Matplotlib is required for this example:

```bash
python -m pip install matplotlib
```

Run the example:

```bash
PYTHONPATH="$PWD/install:$PYTHONPATH" python demo/import_ply.py
```

The script prints the voxel grid dimensions and other information about the voxelization process.


## Visualization examples

The following figures illustrate FastVoxel voxelization results and the resulting voxel representation.

![FastVoxel voxel matrix example](https://github.com/user-attachments/assets/2f2e9a6c-82f9-4ef0-afb1-67f62e3232ac)
*3D visualization of the original PLY model used as input for voxelization.*

![FastVoxel voxelization example](https://github.com/user-attachments/assets/fb4db418-1386-4900-a688-5f20a304579c)
*2D vertical cross-section of the voxelized computational scene, showing the distribution of material identifiers within the voxel grid.*