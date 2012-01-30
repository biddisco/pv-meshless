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

#ifndef _ELB_ELM_CONST_H
#define _ELB_ELM_CONST_H
#ifndef lint
static char *cvs_elemch_id = "$Id$";
#endif

/* Define element types */
typedef enum {SPHERE, BAR2, BAR3, QUAD4, QUAD8, QUAD9,
              SHELL4, SHELL8, TRI3, TRI6, TSHELL3, TSHELL6, HEX8,
              HEX20, HEX27, HEXSHELL, TET4, TET10, TET8, WEDGE6, WEDGE15,
              WEDGE16, PYRAMID5, PYRAMID13, SHELL2, NULL_EL} E_Type;

extern char* elem_names[NULL_EL];

extern
E_Type get_elem_type(
  const char *elem_name,	/* ExodusII element name */
  const int   num_nodes,	/* Number of nodes in the element */
  const int   num_dim		/* Number of dimensions of the mesh */
);

extern
int get_elem_info(
  const int info,		/* The requested information */
  const E_Type elem_type	/* The element type */
);

extern
int numbermatch(
  int* sidenodes, 
  int i, 
  int j, 
  int k,
  int value );

extern
int get_side_id(
  const E_Type  etype,		/* The element type */
  const int *conn,		/* The element connectivity */
  const int  nsnodes,		/* The number of side nodes */
  int  side_nodes[],	        /* The list of side node IDs */
  const int  skip_check,        /* if set, then don't return error if
                                   some nodes aren't in connect table */
  const int partial_adj         /* if set, used partial method of 
                                   determining adjacency */
);

extern
int get_side_id_hex_tet(
  const E_Type  etype,          /* The element type */
  const int *conn,              /* The element connectivity */
  const int  nsnodes,           /* The number of side nodes */
  const int  side_nodes[]       /* The list of side node IDs */
);

extern
int ss_to_node_list(
  const E_Type  etype,		/* The element type */
  const int *connect,		/* The element connectivity */
  int  side_num,		/* The element side number */
  int  ss_node_list[]		/* The list of side node IDs */
);

extern
int get_ss_mirror(
  const E_Type etype,		/* The element type */
  const int *ss_node_list,	/* The list of side node IDs */
  int side_num,			/* The element side number */
  int mirror_node_list[]	/* The list of the mirror side node IDs */
);


/* Define element info requests */
#define NNODES		0
#define NQUAD		1
#define NDIM		2
#define NQUAD_SURF	3
#define NSIDE_NODES	4
#define NSIDES		5

/* Define for the maximum number of nodes on an element side/face */
#define MAX_SIDE_NODES	9
/*
 * Define for the maximum number of sides (and hence communications
 * entries) that an element can have
 */
#define MAX_ELEM_SIDES	6

#endif /* _ELB_ELM_CONST_H */
