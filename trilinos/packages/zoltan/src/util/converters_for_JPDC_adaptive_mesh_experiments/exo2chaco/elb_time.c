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
#ifndef lint
static char *cvs_time_id = "$Id$";
#endif

#include <time.h>

#ifdef  HAS_GETRUSAGE
#include <sys/resource.h>
#endif

/*****************************************************************************/
/* This function returns the time
 *****************************************************************************/
double get_time(void)
{
  double timer;
#ifdef  HAS_GETRUSAGE
  struct rusage timeval;
  getrusage(RUSAGE_SELF, &timeval);
  timer = ((timeval.ru_utime.tv_sec + timeval.ru_stime.tv_sec) +
           1.0e-6 * (timeval.ru_utime.tv_usec + timeval.ru_stime.tv_usec));
#else
  /* ANSI timer, but lower resolution & wraps around after ~36 minutes. */
  timer = clock()/((double) CLOCKS_PER_SEC);
#endif
  return timer;
}
