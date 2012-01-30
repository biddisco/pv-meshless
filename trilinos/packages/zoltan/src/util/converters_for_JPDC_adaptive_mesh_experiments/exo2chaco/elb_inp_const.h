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

#ifndef _ELB_INP_CONST_H_
#define _ELB_INP_CONST_H_
#ifndef lint
static char *cvs_inp_consth_id = "$Id$";
#endif

#include "elb_const.h"

/* Prototype for command-line parsing function */
extern
int cmd_line_arg_parse(
  int argc, 		/* The command line argument count */
  char *argv[],		/* The command line arguments array */
  char ex_inp_file[],	/* The ExodusII input FEM file name */
  char as_inp_file[],	/* The ASCII input file name */
  char ne_out_file[],	/* The output NemesisI file name */
  MACHINE_PTR machine,	/* Pointer to structure in which to place machine
                         * information */
  LB_INFO_PTR lb,	/* Pointer to structure in which to place load
                         * balance parameters */
  PROB_INFO_PTR prob,	/* Pointer to structure in which to place general
                         * information about the run */
  SOLVE_INFO_PTR sol,	/* Pointer to structure in which to place parameters
                         * for the eigensolver */
  WEIGHT_INFO_PTR wgh	/* Pointer to structure in which to place parameters
                         * for the graph weighting scheme */
  );

/* Prototype for function which reads in the ASCII input file */
extern
int read_cmd_file(
  char as_inp_file[],	/* The ASCII input file name */
  char ex_inp_file[],	/* The ExodusII input FEM file name */
  char ne_out_file[],	/* The output NemesisI file name */
  MACHINE_PTR machine,	/* Pointer to structure in which to place machine
                         * information */
  LB_INFO_PTR lb,	/* Pointer to structure in which to place load
                         * balance parameters */
  PROB_INFO_PTR prob,	/* Pointer to structure in which to place general
                         * information about the run */
  SOLVE_INFO_PTR sol,	/* Pointer to structure in which to place parameters
                         * for the eigensolver */
  WEIGHT_INFO_PTR wght	/* Pointer to structure in which to place parameters
                         * for the eigensolver */
);

/* Prototype for function which checks the user specified input */
extern
int check_inp_specs(
  char ex_inp_file[],	/* The ExodusII input FEM file name */
  char ne_out_file[],	/* The output NemesisI file name */
  MACHINE_PTR machine,	/* Pointer to structure in which to place machine
                         * information */
  LB_INFO_PTR lb,	/* Pointer to structure in which to place load
                         * balance parameters */
  PROB_INFO_PTR prob,	/* Pointer to structure in which to place general
                         * information about the run */
  SOLVE_INFO_PTR sol,	/* Pointer to structure in which to place parameters
                         * for the eigensolver */
  WEIGHT_INFO_PTR wght	/* Pointer to structure in which to place parameters
                         * for the weighting scheme */
);

/* Various defines used by the input routines */
#define NONE -1

#endif /* _ELB_INP_CONST_H_ */
