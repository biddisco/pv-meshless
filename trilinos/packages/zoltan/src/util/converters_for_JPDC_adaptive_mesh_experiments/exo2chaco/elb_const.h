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

#ifndef _EXOIILB_CONST_H_
#define _EXOIILB_CONST_H_

#ifndef lint
static char *cvs_exoIIlbch_id = "$Id$";
#endif

/* This is for the maximum file name length */
#include <stdio.h>
#include <exodusII.h>
#include "elb_elem_const.h"

/* Maximum length of a filename */
#define MAX_FNL 1024

#define ELB_VERSION	"3.25"
#define UTIL_NAME	"nem_slice"
#define ELB_FALSE	0
#define ELB_TRUE	1

/* Macro for maximum value */
#ifndef MAX
#define MAX(x,y)  ((x > y) ? x : y)
#endif

/*
 * Constants for memory allocation of graph structures. The smaller these
 * values, the more efficient the code will be, memory wise. Larger values
 * will likely speed execution and prevent swap thrashing.
 */
#define SURND_ALLOC	8
#define ADJ_ALLOC	8
#define MEM_CHUNK_SIZE	16	/* Value MUST be >= 2 */
#define MEM_GROWTH      1.5

#define MAX_INP_LINE 10240

/* Prototype for timing function */
extern double get_time(void);

/* Structure used for the description of the machine for which the
 * load balance is to be constructed. */
struct Machine_Description
{
  int type;
  int num_dims;
  int dim[3];
  int num_boxes;     /* added for cluster type machines */
  int procs_per_box; /* added for cluster type machines, if there is only
                        one box, then this is the same as num_procs */
  int num_procs;
};

typedef struct  Machine_Description  MACHINE_TYPE;
typedef struct  Machine_Description *MACHINE_PTR;

/* Structure used for the description of what type of load balance is
 * to be performed. */
struct LB_Description
{
  int      type;
  int      refine;
  int      num_sects;
  int      cnctd_dom;
  int      outfile;
  char     file[MAX_FNL+1];

  /* Calculated quantities */
  short   *vertex2proc;

  /* Nodal */
  int    **int_nodes;
  int    **bor_nodes;
  int    **ext_nodes;
  int    **ext_procs;
  int     *num_int_nodes;
  int     *num_bor_nodes;
  int     *num_ext_nodes;

  /* Elemental */
  int  ***born_procs;
  int   **born_proc_cnts;
  int   **int_elems;
  int   **bor_elems;
  int   **e_cmap_elems;
  int   **e_cmap_sides;
  int   **e_cmap_procs;
  int   **e_cmap_neigh;
  int    *e_cmap_size;
  int    *num_int_elems;
  int    *num_bor_elems;

};

typedef struct LB_Description  LB_INFO;
typedef struct LB_Description *LB_INFO_PTR;

/* Structure for the problem description. */
struct Problem_Description
{
  int   type;
  int   read_coords;
  int   coarse_flag;
  int   alloc_graph;
  int   num_vertices;
  int   vis_out;
  int   skip_checks;      /* put in to skip some error checks for some meshes  */
  int   face_adj;         /* true if using face definition of adjacencies      */
  int   partial_adj;      /* true id allowing partial (3/4) of nodes to */
                          /* determine adjancencies */
  int   global_mech;      /* true if looking for mechanisms in original mesh   */
  int   local_mech;       /* true if looking for mechanisms in subdivided mesh */
  int   find_cnt_domains; /* true if finding number of connected domains in a graph */
  int   mech_add_procs;   /* adds processors in cases of mechanisms       */
  int   dsd_add_procs;    /* adds processors in cases of disconnected subdomains */
  int   no_sph;
  char *groups;
  int  *group_no;
  int   num_groups;
};

typedef struct Problem_Description  PROB_INFO;
typedef struct Problem_Description *PROB_INFO_PTR;

/* Structure for parameters needed for the Eigensolver in Chaco */
struct Solver_Description
{
  double tolerance;
  int    rqi_flag;
  int    vmax;
};

typedef struct Solver_Description  SOLVE_INFO;
typedef struct Solver_Description *SOLVE_INFO_PTR;

/* Structure used to store information about the weighting scheme, if
 * any, that is to be used. */
struct Weight_Description
{
  int   type;  /* See weight type below for possible types */
  int   ow_read; /* 1 if element block settings overwrite exodus file read */

  /* MAX_FNL, on most architectures, is defined in stdio.h. If not
   * the default value at the top of this file is used. */
  char  exo_filename[MAX_FNL+1];

  /* MAX_STR_LENGTH is defined in exodusII.h */
  char  exo_varname[MAX_STR_LENGTH];
  int   exo_tindx;
  int   exo_vindx;

  /* Variable parameters */
  int    nvals;

  /* vectors to hold element block weights */
  int   num_ebw;
  int   *elemblk;
  int   *elemblk_wgt;

  /* vector to indicate if weight value has already been overwritten */
  int   *ow;

  int   *vertices;

  float *edges;
};

typedef struct Weight_Description  WEIGHT_INFO;
typedef struct Weight_Description *WEIGHT_INFO_PTR;


/* Structure used to store information about the FEM mesh */
struct Mesh_Description
{
  int	  num_nodes;
  int     num_elems;
  int	  num_dims;
  int     num_el_blks;
  int    *eb_cnts;
  int     num_node_sets;
  int     num_side_sets;
  int     max_np_elem;
  int     ns_list_len;
  char    title[MAX_LINE_LENGTH+1];
  float  *coords;
  E_Type *elem_type;
  int   **connect;
};
typedef struct Mesh_Description  MESH_INFO;
typedef struct Mesh_Description *MESH_INFO_PTR;

/* Structure for handling meshes with spheres */
struct Sphere_Info
{
  int   num;
  int  *adjust;
  int  *begin;
  int  *end;
};
typedef struct Sphere_Info  SPHERE_INFO;
typedef struct Sphere_Info *SPHERE_INFO_PTR;

/* Structure used to store various information about the graph */
struct Graph_Description
{
  int     *adj;
  int     *start;
  long     nadj;
  long     adj_alloc;
  int    **sur_elem;
  int      max_nsur;
  int     *alloc_cnt;
  int     *nsur_elem;
};
typedef struct Graph_Description  GRAPH_INFO;
typedef struct Graph_Description *GRAPH_INFO_PTR;

/* Various constants */
#define NODAL 0
#define ELEMENTAL 1

#define UTIL_NAME "nem_slice"

/* Load balance types */
#define MULTIKL		0
#define SPECTRAL	1
#define INERTIAL	2
#define LINEAR 		3
#define RANDOM 		4
#define SCATTERED 	5
#define INFILE		6
#define KL_REFINE 	7
#define NO_REFINE 	8
#define NUM_SECTS	9
#define CNCT_DOM	10
#define OUTFILE		11
#define ZPINCH          12
#define BRICK           13
#define ZOLTAN_RCB      14
#define ZOLTAN_RIB      15
#define ZOLTAN_HSFC     16
#define GEN_CHACO       17

/* Machine types */
#define MESH 		0
#define HCUBE		1
#define HYPERCUBE	2
#define CLUSTER		3

/* Solver options */
#define TOLER 		0
#define USE_RQI 	1
#define VMAX		2

/* ISSUES options */

#define LOCAL_ISSUES 0
#define GLOBAL_ISSUES 1

/* Weighting options */
/*
 * NOTE: the position of NO_WEIGHT, READ_EXO, EL_BLK, and EWGT_ON
 * should not be changed. These are the possible values for the
 * "type" variable in the Weight struct. They need to b 0, 1, 2, & 4
 * to allow bit masking for the type. The other variables are not
 * currently used in the type, but are needed since they appear
 * on the command line.
 */
#define NO_WEIGHT	0
#define READ_EXO	1
#define EL_BLK		2
#define VAR_INDX	3
#define EDGE_WGT	4
#define TIME_INDX	5
#define VAR_NAME	6

extern int identify_mechanisms(
  MACHINE_PTR machine,
  PROB_INFO_PTR problem,
  MESH_INFO_PTR mesh,
  LB_INFO_PTR lb,
  GRAPH_INFO_PTR graph,
  int check_type
);

#endif /* _EXOIILB_CONST_H_ */
