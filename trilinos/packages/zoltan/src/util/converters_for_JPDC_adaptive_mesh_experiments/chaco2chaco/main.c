/*  Program to map a set a coordinates in one file to the closest coordinates
 *  in a second file.
 */

#include <stdio.h>
#include <stdlib.h>
#include "hsfc_hilbert_const.h"

#define MYDBLMAX  1000000000.
#define MYDBLMIN -1000000000.

struct KEY {
  int idx;
  double hsfcKey;
};

/***************************************************************************/
void readfile(
  char *fileroot,
  int filenum,
  int *n,
  struct KEY **key
)
{
char fname[65];
int j;
FILE *fp;
typedef double COORDINATES[3];
COORDINATES *coor;
COORDINATES max = {MYDBLMIN, MYDBLMIN, MYDBLMIN};
COORDINATES min = {MYDBLMAX, MYDBLMAX, MYDBLMAX}; 
COORDINATES delta;

  /* Get number of nodes in file from chaco graph file */
  sprintf(fname, "%s.%d.graph", fileroot, filenum);
  fp = fopen(fname, "r");
  fscanf(fp, "%d", n);
  fclose(fp);

  /* Get coordinate info from chaco coords file */
  printf("READING FILE %d  with %d nodes\n", filenum, *n);
  *key = (struct KEY *) malloc((1+*n) * sizeof(struct KEY));
  coor = (COORDINATES *) malloc(*n * sizeof(COORDINATES));
  sprintf(fname, "%s.%d.coords", fileroot, filenum);
  fp = fopen(fname, "r");
  for (j = 0; j < *n; j++) {
    fscanf(fp, "%lf %lf %lf", &(coor[j][0]), &(coor[j][1]), &(coor[j][2]));
    if (coor[j][0] < min[0]) min[0] = coor[j][0];
    if (coor[j][0] > max[0]) max[0] = coor[j][0];
    if (coor[j][1] < min[1]) min[1] = coor[j][1];
    if (coor[j][1] > max[1]) max[1] = coor[j][1];
    if (coor[j][2] < min[2]) min[2] = coor[j][2];
    if (coor[j][2] > max[2]) max[2] = coor[j][2];
  }
  delta[0] = max[0] - min[0];
  delta[1] = max[1] - min[1];
  delta[2] = max[2] - min[2];
  for (j = 0; j < *n; j++) {
    /* Scale to unit cube */
    coor[j][0] = (coor[j][0] - min[0]) / delta[0];
    coor[j][1] = (coor[j][1] - min[1]) / delta[1];
    coor[j][2] = (coor[j][2] - min[2]) / delta[2];
    (*key)[j].hsfcKey = Zoltan_HSFC_InvHilbert3d(coor[j]);
    (*key)[j].idx = j;
  }
  (*key)[*n].hsfcKey = MYDBLMAX;
  (*key)[*n].idx = *n;
  free(coor);
  fclose(fp);
}
/***************************************************************************/
/* Sorting values in increasing order. Criteria is double. */
static void quickpart_list_inc_double(
  struct KEY *list, int start, int end, int *equal, int *larger)
{
int i;
struct KEY key, change;

/*  key = list ? list[(end+start)/2] : 1; */
  key = list[(end+start)/2];

  *equal = *larger = start;
  for (i = start; i <= end; i++)
    if (list[i].hsfcKey < key.hsfcKey) {
      change            = list[i];
      list[i]           = list[*larger];
      list[(*larger)++] = list[*equal];
      list[(*equal)++]  = change;
    }
    else if (list[i].hsfcKey == key.hsfcKey) {
      list[i]           = list[*larger];
      list[(*larger)++] = key;
    }
}

void Zoltan_quicksort_list_inc_double(struct KEY *list, int start, int end)
{
int  equal, larger;

  if (start < end) {
    quickpart_list_inc_double (list, start, end, &equal, &larger);
    Zoltan_quicksort_list_inc_double (list, start, equal-1);
    Zoltan_quicksort_list_inc_double (list, larger, end);
  }
}

/***************************************************************************/


/***************************************************************************/
/* routine for bsearch locating the nearest keys to the given key */
int HSFC_compare (const void *key, const void *arg)
{
  if (((struct KEY*)key)->hsfcKey >= ((struct KEY*)(arg)+1)->hsfcKey) return 1;
  if (((struct KEY*)key)->hsfcKey < ((struct KEY*)arg)->hsfcKey) return -1;

  return 0;     /* key in interval [*arg, *(arg+1)] */
}


/***************************************************************************/
main(int argc, char *argv[])
{
int no, nn;
char fileroot[65];
char fname[65];
int numfiles;
struct KEY *ko, *kn;
FILE *fp;
struct KEY *closest;
int i, j, k;

  if (argc < 3) {
    printf("Usage:  a.out filenameroot maxfilenumber");
  }
  strcpy(fileroot, argv[1]);
  numfiles = atoi(argv[2]);

  /* Get info for file zero. */
  readfile(fileroot, 0, &no, &ko);
  Zoltan_quicksort_list_inc_double(ko, 0, no);

  for (i = 1; i <= numfiles; i++) {
    readfile(fileroot, i, &nn, &kn);

    /* For each node in file i, find nearest node in file i-1 */
    sprintf(fname, "%s.%d.nearest", fileroot, i);
    fp = fopen(fname, "w");

    for (j = 0; j < nn; j++) {
      /* Find old node closest to this node */
      closest = (struct KEY *) bsearch(&kn[j], ko, no, 
                                       sizeof(struct KEY), HSFC_compare);
      fprintf(fp, "%d\n", closest->idx + 1);
    }
    fclose(fp);
    
    /* Get ready for next file */
    free(ko);
    no = nn; ko = kn;
    if (i < numfiles) Zoltan_quicksort_list_inc_double(ko, 0, no);
  }
  free(ko);
}
