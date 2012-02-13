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

#include "spatial_discretization.hpp"
#include <vector>
#include <string>

#ifndef __SCALARFIELDBUILDERS__
#define __SCALARFIELDBUILDERS__

namespace ScalarFieldBuilders
{
	/**
	 * Retourne les coordonnées du centre du cube correspondant à l'indice en paramètre
	 */
	vec3 CellIdToCenterCoordinate( const ivec3& cell_id, const decimal& cellSize, const vec3& zeroCellCenter);
	/**
	 * Cette classe permet de générer un espace discrétisé en plusieurs volumes. C'est la première étape de la reconstruction du modèle.
	 */
	class ScalarFieldCreator
	{
	protected:
		PTR<SpatialDiscretization::weight_matrix> fieldData; //Données de la matrice X,Y,Z


		struct mainVolumeConstruction_t
		{
		    mainVolumeConstruction_t()
		    : cellSize(0.),cellCount(0),volumeCount(0),maximal_marker_index(0)
            {
            }
			decimal cellSize;
			unsigned int cellCount;
			vec3 mainVolumeCenter;
			vec3 cellHalfSize;
			vec3 zeroCellCenter;
			SpatialDiscretization::weight_t volumeCount;
			std::vector<decimal> volumeValue;
			vec3 boxMin;
			vec3 boxMax;
			SpatialDiscretization::weight_t maximal_marker_index;
		} volumeInfo;
		SpatialDiscretization::domainInformation_t domainInformation;
		decimal resolution;
		static void ComputeMatrixParams(const vec3& boxMin,const vec3& boxMax, const decimal& minResolution, mainVolumeConstruction_t& computedVolumeInfo);
		/**
		 * Initialise les données pour le volume extérieur
		 */
		void InitExteriorVolumeId();
		/**
		 * Propage les valeurs des identifiants de volume d'une cellule à l'autre
		 * Si la cellule cible est modifié cette méthode retourne vrai
		 */
		bool CellToCellVolumePropagation(const ivec2& destinationPropa,const ivec2& sourcePropa,const SpatialDiscretization::weight_t& volumeId);
		/**
		 * Propage un indice de volume dans toute la matrice
		 */
		void ExtandVolume(const SpatialDiscretization::weight_t& volumeId);
		/**
		 * Retourne la position de la première cellule avec la valeur en paramètre
		 */
		ivec3 GetFirstCellByWeight(const SpatialDiscretization::weight_t& weight,SpatialDiscretization::zcell** foundCell);

		/**
		 * Calcul pour chaque volume sa valeur en m^3
		 * @param[out] volumeValue Un tableau de dimension égale au nombre de volume dans le domaine. Dont la valeur est en m^3.
		 */
		void ComputeVolumesValue(std::vector<decimal>& volumeValue);


	public:
		/**
		 * Constructeur
		 * @param _resolution Dimension d'une cellule qui composera la matrice. Plus la résolution est élevée plus le model généré sera proche du modèle en entrée et plus de triangles seront générés.
		 */
		ScalarFieldCreator(const decimal& _resolution);
		/**
		 * Initialisation de la matrice selon la résolution et la boite englobante passé en paramètre.
		 * @param boxMin Coordonnées minimale des objets qui alimenteront la matrice
		 * @param boxMax Coordonnées maximale des objets qui alimenteront la matrice
		 */
		void FirstStep_Params(const vec3& boxMin,const vec3& boxMax);
		virtual ~ScalarFieldCreator();

		/**
		 * Une fois toutes les primitives renseignées. Cette méthode doit être appelée afin de détecter les volumes délimité par les limites.
		 */
		void ThirdStep_VolumesCreator();

		/**
		 * Retourne la valeur de la matrice selon les indices des cellules
		 * @param index Entier positif désignant le n° de cellule.
		 * @see GetDomainSize()
		 */
		SpatialDiscretization::weight_t GetMatrixValue(const ivec3& index);


		vec3 GetCenterCellCoordinates( const ivec3& cell_id) const;


        /**
         * Get the first cell value associated to matrix, this index start from the highest material index of materials
         */
        SpatialDiscretization::weight_t GetFirstVolumeIndex();

		/**
		 * Pour toutes les valeurs de Z pour un x,y donné. Retourne la valeur min,max.
		 */
		void GetMinMaxOnZ( const ivec2& xyCell, SpatialDiscretization::weight_t& minVolId,SpatialDiscretization::weight_t& maxVolId);
		bool IsContainsVol( const ivec2& xyCell, SpatialDiscretization::weight_t& volId);

        /**
         * Get the volume count
         */
		SpatialDiscretization::weight_t GetVolumeCount();
		/**
		 * Get the volume (m3) corresponding to volume ID [0-GetVolumeCount()[
		 */
		decimal GetVolumeValue(const SpatialDiscretization::weight_t& volId);

		void GetMinMax(vec3& minBox,vec3& maxBox);

		/**
		 * Exporte les indices et volumes des domaines
		 * @param fileName Nom et chemin du fichier de sortie
		 * @param volsLabelsFileName Nom et chemin du fichier d'entrée des noms utilisateur de volume (fichier texte avec "x y z NomDuVolume" à chaque ligne)
		 */
		void ExportVolsStats(const std::string& fileName, const std::string& volsLabelsFileName=std::string());
		///////////////////////////////
		// Debug Functions
		unsigned int count();
		void MakeXYZ(const std::string& filename,const SpatialDiscretization::weight_t& idVol);
		void ExportIJKData(const std::string& Infilename,const std::string& Outfilename);
		void MakeXYZ(const std::string& filename,const decimal& minVol);
		void ExportVTK(const std::string& filename,const SpatialDiscretization::weight_t& idVol=-1);
		std::size_t GetDomainSize();
		bool CheckDiscretisation();
		SpatialDiscretization::weight_t GetLargestVolumeId();

        /**
         * Copy data from spatial discretisation to array in argument
         * @param data Previously allocated data
         * @param ni x maximum size
         * @param nj y maximum size
         * @param nk z maximum size
         * @param extractPos zero position extraction
         */
        void CopyMatrix(SpatialDiscretization::weight_t* data,int ni,int nj,int nk,const ivec3& extractPos);
        void CopyMatrixFiltered(SpatialDiscretization::weight_t* data,int ni,int nj,int nk,const ivec3& extractPos,const SpatialDiscretization::weight_t* data_filter,int nindex );
        /**
		 * Retourne l'indice de la cellule contenant le point passé en paramètre
		 */
		ivec3 GetCellIdByCoord(const vec3& position);

		void GetCellValueBoundaries(ivec3& min,ivec3& max,const SpatialDiscretization::weight_t& volid);
	};

}

#endif
