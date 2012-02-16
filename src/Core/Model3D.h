/*
 *     This file is part of FastVoxel.
 *
 *     FastVoxel is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     FastVoxel is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *     along with FastVoxel.  If not, see <http://www.gnu.org/licenses/>.
 * FastVoxel is a voxelisation library of polygonal 3d model and do volumes identifications.
 * It is dedicated to finite element solvers
 * @author Nicolas Fortin
 * This project is the production of IFSTTAR (www.ifsttar.fr)
 * @copyright GNU Public License.V3
 * Official repository is https://github.com/nicolas-f/FastVoxel
 */
#ifndef _HMODEL
#define _HMODEL

/**
 * @file Model3D.h
 *
 * @brief Définition de types en rapport avec les modèles 3D
 */
#include <vector>

/**
 * @brief Indice d'une face
 *
 * Cette classe gère l'indexation d'une face
 */
class t_faceIndex {
public:
	t_faceIndex(void) : f(0), g(0) { }
	t_faceIndex(long _face,long _group) : f(_face), g(_group) { }
	t_faceIndex(const t_faceIndex *_faceIndex) : f(_faceIndex->f), g(_faceIndex->g) { }
	t_faceIndex(const t_faceIndex &_faceIndex) : f(_faceIndex.f), g(_faceIndex.g) { }

	int operator==(const t_faceIndex &_faceIndex) { return (this->f == _faceIndex.f && this->g == _faceIndex.g ); }
	int operator!=(const t_faceIndex &_faceIndex) { return !(*this == _faceIndex); }

	void Set(long _face,long _group) { f=_face;g=_group; }
	void Set(const t_faceIndex *_faceIndex) { f=_faceIndex->f;g=_faceIndex->g; }
	void Set(const t_faceIndex &_faceIndex) { f=_faceIndex.f;g=_faceIndex.g;}

	union {
		struct {long f,g;};
	};
};

/**
 * Structure d'un triangle
 */
struct triangleFace {
    float a[3];
	float b[3];
	float c[3];
};

/**
 * Structure complète d'une face en 3D
 */
struct SFace3D
{
	ivec3 Vertices;
	vec3 FaceNormals;
};

/**
 * Structure complète d'un groupe de surface
 */
struct SGroup3D
{
	std::vector<SFace3D> pFaces;
	long numOfVerts;
};

#endif
