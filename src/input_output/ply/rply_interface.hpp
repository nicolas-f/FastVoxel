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


/*! \file ply.h
    \brief Implémentation de l'interpréteur de fichier modèle Ply (*.ply)
	@see
*/
#include <Core/mathlib.h>
#include <list>
#include <string>


#ifndef _HRPLY
#define _HRPLY

/*! \brief Implémentation de l'interpréteur de fichier modèle Poly (*.poly)
*/
namespace formatRPLY
{

	/**
	 * @brief Structure de données d'echange avec la classe
	 */
	struct t_face
	{
		t_face(const ivec3& _indicesSommets) : indicesSommets(_indicesSommets) {}
		ivec3 indicesSommets;
	};


	/**
	 * @brief Structure de données d'un calque
	 */
	struct t_layer
	{
		t_layer(){};
		t_layer(const std::string& _layerName):layerName(_layerName) {}
		std::string layerName;
	};

	/**
	 * @brief Structure de données du modèle
	 *
	 */
	struct t_model
	{
		std::list<t_face> modelFaces;
		std::list<vec3> modelVertices;
		std::list<t_layer> modelLayers;				 /*!< Liste des calques */
		std::list<std::size_t> modelFacesLayerIndex; /*!< Correspondance Indice de Face->Indice de calque*/
	};

/**
 *	\class Cply
 *	\brief Classe de sauvegarde et de chargement de fichier PLY (stanford)
 */
class CPly
{
public:
	/**
	 * Méthode d'importation d'un modèle 3D
	 */
	static bool ImportPly(t_model& sceneconst, std::string mfilename);
	/**
	 * Méthode d'exportation d'un modèle 3D
	 */
	static bool ExportPly(t_model& scene, std::string mfilename);
};





} //Fin namespace
#endif
