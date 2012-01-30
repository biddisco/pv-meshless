#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/times.h>
#include <time.h>

#define __USE_XOPEN
#include <stdio.h>

void add_to_log(const char *my_name)
{
#define LEN 128
  char time_string[LEN];
  char log_string[LEN];

  double u_time, s_time;
  struct utsname sys_info;
  
  char *username = NULL;
  const char *codename = NULL;

  /* Don't log information if this environment variable is set */
  if (getenv("SEACAS_NO_LOGGING") != NULL) {
    fprintf(stderr, "SEACAS Audit logging disabled via SEACAS_NO_LOGGING setting.\n");
    return;
  }
  
  username = getlogin();
  if (username == NULL) {
    username = getenv("LOGNAME");
  }

  codename = strrchr (my_name, '/');
  if (codename == NULL)
    codename = my_name;
  else
    codename++;

  {
    time_t calendar_time = time(NULL);
    struct tm *local_time = localtime(&calendar_time);
    strftime(time_string, LEN, "%a %b %d %H:%M:%S %Z %Y", local_time);
  }

  {
    int ticks_per_second;
    struct tms time_buf;
    times(&time_buf);
    ticks_per_second = sysconf(_SC_CLK_TCK);
    u_time = (double)(time_buf.tms_utime + time_buf.tms_cutime) / ticks_per_second;
    s_time = (double)(time_buf.tms_stime + time_buf.tms_cstime) / ticks_per_second;
  }
  
  uname(&sys_info);

  sprintf(log_string, "%s %s %s %.3fu %.3fs 0:00.00 0.0%% 0+0k 0+0io 0pf+0w %s\n",
	  codename, username, time_string, u_time, s_time, sys_info.nodename);

  /* Now try to find the $ACCESS/etc/audit.log file */
  /* Don't need to try too hard since information is not critical; just useful */
  {
    char *access_dir = getenv("ACCESS");
    if (access_dir != NULL) {
      char filename[LEN];
      sprintf(filename, "%s/etc/audit.log", access_dir);
      if (0 == access(filename, W_OK)) {
	FILE *audit = fopen(filename, "a");
	if (audit != NULL) {
	  fprintf(audit, "%s", log_string);
	}
      }
    }
  }
}
