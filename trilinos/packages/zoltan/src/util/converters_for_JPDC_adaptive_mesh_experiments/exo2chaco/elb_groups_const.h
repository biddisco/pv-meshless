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

#ifndef _ELB_GROUPS_CONST_H_
#define _ELB_GROUPS_CONST_H_

#ifndef lint
static char *cvs_grpch_id = "$Id$";
#endif
#include "elb_const.h"

/* Function prototypes */
extern
int parse_groups(
  int *el_blk_ids,		/* array containing element block ids */
  int *el_blk_cnts,		/* array containing element block counts */
  MESH_INFO_PTR mesh,		/* Mesh information structure */
  PROB_INFO_PTR prob		/* Problem information */
);

extern
int get_group_info(
  MACHINE_PTR machine,
  PROB_INFO_PTR prob,
  MESH_INFO_PTR mesh,
  GRAPH_INFO_PTR graph,
  short elem2grp[],
  int nprocg[],
  int nelemg[],
  int *max_vtx,
  int *max_adj
);

#endif /* _ELB_GROUPS_CONST_H_ */
