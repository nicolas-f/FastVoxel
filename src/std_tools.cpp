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
#include "std_tools.hpp"
#include <math.h>
#include <float.h>
#ifdef _WIN32
	#include <direct.h>
#endif
#ifdef _UNIX
	#include <sys/stat.h>
    #include <sys/types.h>
#endif

int st_mkdir(const char* pathname, int perm)
{
	#ifdef _WIN32
		return mkdir(pathname);
	#endif
	#ifdef _UNIX
		return mkdir(pathname,perm);
	#endif
}


bool st_isfinite(const float& value)
{
	#ifdef _MSC_VER
		return _finite(value);
	#else
		return isfinite(value);
	#endif
}
