#ifndef __OCTREE44_TRIELEMENT_H__
#define __OCTREE44_TRIELEMENT_H__

/**
 * Indique si un triangle et un cube se superposent.
 *
 * @param[in] boxcenter Centre de la boite
 * @param[in] boxhalfsize Dimension de la moitié de la boite
 * @param[in] triverts Position des sommets du triangle
 * @return 0 si aucun contact, 1 si il y a une superposition
 */
namespace boxtri_test
{
    int triBoxOverlap(
		double boxcenter[3],
		double boxhalfsize[3],
		double triverts[3][3]);
}
#endif