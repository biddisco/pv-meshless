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

#ifndef _ELB_UTIL_CONST_H_
#define _ELB_UTIL_CONST_H_
#ifndef lint
static char *cvs_util_ch_id = "$Id$";
#endif

/* Function prototypes */
extern
void print_usage(void);

extern
int token_compare(
  char *token,		/* The input character string */
  const char *key	/* The key to compare with token */
  );

extern
void strip_string(
  char inp_str[],	/* The string to strip */
  const char *tokens	/* The tokens to strip from the beginning and
                         * end of the input string */
  );

extern
void clean_string(
  char inp_str[],	/* The string to clean */
  const char *tokens	/* The tokens to strip multiple copies of */
  );

extern
void string_to_lower(
  char inp_str[],	/* The string to convert to lower case */
  const char cstop	/* Character where to stop */
  );

extern
void qsort4(int *v1, int *v2, int *v3, int *v4, int N);

extern
void qsort2(int *v1, int *v2, int N);

extern
void sort2_int_int(
  int  count,
  int *array1,
  int *array2
  );

extern
void sort3_int_int_int(
  int  count,
  int *array1,
  int *array2,
  int *array3
  );

extern
void sort4_iiii(
  int  count,
  int *array1,
  int *array2,
  int *array3,
  int *array4
  );

extern
void qsort4(
  int *v1,
  int *v2,
  int *v3,
  int *v4,
  int N
  );

extern
void find_first_last(
  int  value,
  int  vecsize,
  int *vector,
  int *first,
  int *last
  );

extern
int find_int(
  int  value1,
  int  value2,
  int  start,
  int  stop,
  int *vector1,
  int *vector2
  );

extern
int in_list(
  const int  search,		/* The value to search for */
  const int  count,		/* Number of elements in vector to search */
  int       *vector		/* The vector to search */
);

extern
int roundfloat(
  const float value		/* the value to be rounded */
);

extern int find_max (
  const int list_length,
  const int list[]
);

extern int find_min (
  const int list_length,
  const int list[]
);


extern
int find_inter (
  const int set1[],		/* the first set of integers */
  const int set2[],		/* the second set of integers */
  const int length1,		/* the length of the first set */
  const int length2,		/* the length of the second set */
  const int prob_type,		/* value indicating known info about lists */
  int inter_ptr[]		/* the values in the intersection */
);

#endif /* _ELB_UTIL_CONST_H_ */
