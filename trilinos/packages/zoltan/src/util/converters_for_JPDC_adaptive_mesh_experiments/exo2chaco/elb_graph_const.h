/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/

#ifndef _ELB_GRAPH_CONST_H_
#define _ELB_GRAPH_CONST_H_
#ifndef lint
static char *cvs_graphch_id = "$Id$";
#endif

#include "elb_const.h"

extern
int generate_graph(
  PROB_INFO_PTR prob,		/* Pointer to structure containing information
                                 * about the type of decomposition */
  MESH_INFO_PTR mesh,		/* Pointer to structure containing the mesh */
  GRAPH_INFO_PTR graph,		/* Pointer to structure in which to store
                                 * the graph. */
  WEIGHT_INFO_PTR weight,	/* Pointer to structure for graph weighting */
  SPHERE_INFO_PTR sphere	/* Pointer to sphere adjustment structure */
);

#endif /* _ELB_GRAPH_CONST_H_ */
