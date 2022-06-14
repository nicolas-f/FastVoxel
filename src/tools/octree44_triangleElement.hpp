#include <math.h>

#ifndef __OCTREE44_TRIELEMENT_H__
#define __OCTREE44_TRIELEMENT_H__
/**
 * Indique si un triangle et un cube se superpose
 * @param[in] boxcenter Centre de la boite
 * @param[out] boxhalfsize Dimension de la moitiée de la boite
 * @param[out] triverts Position des sommets des triangles
 * @return 0 si aucun contact, 1 si il y a une superposition
 */
namespace boxtri_test
{
    int triBoxOverlap(double_t boxcenter[3],double_t boxhalfsize[3],double_t triverts[3][3]);
}

#endif
