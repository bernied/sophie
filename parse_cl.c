/******************************************************************************
**
** parse_cl.c
**
** Sun Oct 17 23:24:10 2010
** Linux 2.6.21.7-2.fc8xen (#1 SMP Fri Feb 15 12:39:36 EST 2008) i686
** root@domU-12-31-38-00-44-B2 (root)
**
** C file for command line parser
**
** Automatically created by genparse v0.8.7
**
** See http://genparse.sourceforge.net for details and updates
**
******************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include "parse_cl.h"

static struct option const long_options[] =
{
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'V'},
  {"print-path", no_argument, NULL, 'p'},
  {"initial-temp", required_argument, NULL, 'i'},
  {"steps-per-temp", required_argument, NULL, 't'},
  {"cooling-steps", required_argument, NULL, 's'},
  {"cooling-fraction", required_argument, NULL, 'f'},
  {"K", required_argument, NULL, 'k'},
  {"runs", required_argument, NULL, 'r'},
  {"max-verts", required_argument, NULL, 'v'},
  {"randomize", no_argument, NULL, 'z'},
  {"brute-force", no_argument, NULL, 'b'},
  {NULL, 0, NULL, 0}
};

/*----------------------------------------------------------------------------
**
** Cmdline ()
**
** Parse the argv array of command line parameters
**
**--------------------------------------------------------------------------*/

void Cmdline (struct arg_t *my_args, int argc, char *argv[])
{
  extern char *optarg;
  extern int optind;
  int c;
  int errflg = 0;

  my_args->h = false;
  my_args->V = false;
  my_args->p = false;
  my_args->z = false;
  my_args->b = false;

  optind = 0;
  while ((c = getopt_long (argc, argv, "hVpzbi:t:s:f:k:r:v:", long_options, &optind)) != - 1)
    {
      switch (c)
        {
        case 'h':
          my_args->h = true;
          usage (EXIT_SUCCESS, argv[0]);
          break;

        case 'V':
          my_args->V = true;
          break;

        case 'p':
          my_args->p = true;
          break;

        case 'i':
          my_args->i = atof (optarg);
          break;

        case 't':
          my_args->t = atoi (optarg);
          break;

        case 's':
          my_args->s = atoi (optarg);
          break;

        case 'f':
          my_args->f = atof (optarg);
          break;

        case 'k':
          my_args->k = atof (optarg);
          break;

        case 'r':
          my_args->r = atoi (optarg);
          break;

        case 'v':
          my_args->v = atoi (optarg);
          break;

		case 'z':
		  my_args->z = true;
		  break;

		case 'b':
		  my_args->b = true;
		  break;

        default:
          usage (EXIT_FAILURE, argv[0]);

        }
    } /* while */

  if (errflg)
    usage (EXIT_FAILURE, argv[0]);

  if (optind >= argc)
    my_args->optind = 0;
  else
    my_args->optind = optind;
}

/*----------------------------------------------------------------------------
**
** usage ()
**
** Print out usage information, then exit
**
**--------------------------------------------------------------------------*/

void usage (int status, char *program_name)
{
  if (status != EXIT_SUCCESS)
    fprintf (stderr, "Try `%s --help' for more information.\n",
            program_name);
  else
    {
      printf ("\
Usage: %s [OPTION]... [FILE]\n\
\n\
  -h, --help              display this help and exit\n\
  -V, --version           output version information and exit\n\
  -p, --print-path        Print tour for lowest time\n\
Options controlling simulated annealing heuristics:\n\
  -b, --brute-force       Use brute force to solve\n\
  -i, --initial-temp      Set starting temp for simulated annealing\n\
  -t, --steps-per-temp    Steps per temperature change\n\
  -s, --cooling-steps     Numer of cooling steps\n\
  -f, --cooling-fraction  Cooling fraction for exponent\n\
  -k, --K                 K factor\n\
  -r, --runs              Number of runs\n\
  -v, --max-verts         How many verticies before using heuristic\n\
  -z, --randomize         Randomize starting heuristic path\n\
\n", program_name);
    }
  exit (status);
}
