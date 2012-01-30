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

#ifndef _ELB_EXO_CONST_H_
#define _ELB_EXO_CONST_H_

#ifndef lint
static char *cvs_exoch_id = "$Id$";
#endif
#include "elb_const.h"

/* Function prototypes */
extern
int read_exo_weights(
  PROB_INFO_PTR prob,		/* Pointer to problem info structure */
  WEIGHT_INFO_PTR weight	/* Pointer to weight info structure */
);

extern
int read_mesh_params(
  char exo_file[],		/* Name of ExodusII geometry file */
  PROB_INFO_PTR prob,		/* Pointer to problem info structure */
  MESH_INFO_PTR mesh,		/* Mesh information structure */
  SPHERE_INFO_PTR sphere	/* Sphere element info structure */
);

extern
int read_mesh(
  char exo_file[],		/* Name of ExodusII geometry file */
  PROB_INFO_PTR prob,		/* Problem information */
  MESH_INFO_PTR mesh,		/* Mesh information structure */
  WEIGHT_INFO_PTR weight	/* Weight specification structure */
);

extern
int init_weight_struct(PROB_INFO_PTR problem,	/* Problem information */
                       MESH_INFO_PTR mesh,      /* Mesh information structure */
                       WEIGHT_INFO_PTR weight); /* Weight specification structure */

#endif /* _ELB_EXO_CONST_H_ */
