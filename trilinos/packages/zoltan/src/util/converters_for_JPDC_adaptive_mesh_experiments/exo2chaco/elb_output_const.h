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

#ifndef _ELB_OUTPUT_CONST_H_
#define _ELB_OUTPUT_CONST_H_
#ifndef lint
static char *cvs_outputch_id = "$Id$";
#endif

#include "elb_const.h"

extern
int write_nemesis(char *filename,
                  MACHINE_PTR machine,
                  PROB_INFO_PTR problem,
                  MESH_INFO_PTR mesh,
                  LB_INFO_PTR lb,
                  SPHERE_INFO_PTR sphere
);

extern
int write_vis(char *vis_filename,
              char *exo_inp_filename,
              MACHINE_PTR machine,
              PROB_INFO_PTR problem,
              MESH_INFO_PTR mesh,
              LB_INFO_PTR lb
);

#endif /* _ELB_OUTPUT_CONST_H_ */
