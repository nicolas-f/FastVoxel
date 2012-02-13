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
 * @author Nicolas Fortin , Judicaël Picaut judicael.picaut (home) ifsttar.fr
 * Official repository is https://github.com/nicolas-f/FastVoxel
 */
#include <string>


#ifndef __TEMPORARY_ENGLISH_NUMERIC__
#define __TEMPORARY_ENGLISH_NUMERIC__

	/**
	 * Cette classe permet d'utiliser les . comme séparateur de décimale le temps de l'execution du processus
	 */
	class EnglishTemporaryLocale
	{
	protected:
		std::string currentLocale;
	public:
		EnglishTemporaryLocale();
		~EnglishTemporaryLocale();
	};

#endif
