/******************************************************************************
**
** parse_cl.h
**
** Sun Oct 17 23:24:10 2010
** Linux 2.6.21.7-2.fc8xen (#1 SMP Fri Feb 15 12:39:36 EST 2008) i686
** root@domU-12-31-38-00-44-B2 (root)
**
** Header file for command line parser
**
** Automatically created by genparse v0.8.7
**
** See http://genparse.sourceforge.net for details and updates
**
******************************************************************************/

#include <stdio.h>

#ifndef bool
typedef enum bool_t
{
  false = 0, true
} bool;
#endif

/* customized structure for command line parameters */
struct arg_t
{
  bool h;
  bool V;
  bool p;
  bool z;
  bool b;
  float i;
  int t;
  int s;
  float f;
  float k;
  int r;
  bool e;
  int v;
  int optind;
};

/* function prototypes */
void Cmdline (struct arg_t *my_args, int argc, char *argv[]);
void usage (int status, char *program_name);
