/* Copyright (C) 2010 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <dlfcn.h>
#include <sys/types.h>
#include <dirent.h>
#define SYSLOG
#include <exception.h>
#include <unistd.h>
#include <getopt.h>
#include <version.h>
#include <limits.h> /* PATH_MAX */

#include "bayesian-binning.h"

void
print_usage(char *pname, FILE *fp)
{

	(void)fprintf(fp,
                      "\nUsage: %s [OPTION]... SERVICE -- ...\n\n", pname);
	(void)fprintf(fp,
                      "Options:\n"
                      "   -L		list installed plugins (services) and exit\n"
                      "   -Y		use syslog\n"
                      "   --help	print help and exit\n"
                      "   --version	print version information and exit\n\n"
                      "Report bugs to <"EMAIL">.\n\n");

	return;

}

void
print_version(FILE *fp)
{

	(void)fprintf(fp,
                      "\n" PNAME " version " VERSION ", Copyright (C) 2010 Philipp Benner.\n\n"

                      "Compiled for " HOST ".\n\n"

                      "This is free software, and you are welcome to redistribute it\n"
                      "under certain conditions; see the source for copying conditions.\n"
                      "There is NO warranty; not even for MERCHANTABILITY or FITNESS\n"
                      "FOR A PARTICULAR PURPOSE.\n\n");

	return;

}

void
wrong_usage(const char *msg)
{

	if(msg != NULL) {
		(void)fprintf(stderr, "%s\n", msg);
	}
	(void)fprintf(stderr,
		      "Try `bayesian-binning --help' for more information.\n");

	exit(EXIT_FAILURE);

}

int
main(int argc, char *argv[])
{

	char *service;

	if(argc == 1) {
		wrong_usage("too few arguments");
		exit(EXIT_FAILURE);
	}

	for(;;) {
		int c, option_index = 0;
		static struct option long_options[] = {
			{ "help",	0, 0, 'h' },
			{ "version",	0, 0, 'v' }
		};

		c = getopt_long(argc, argv, "Y",
				long_options, &option_index);

		if(c == -1) {
			break;
		}

		switch(c) {
                case 'Y':
			use_syslog = 1;
			break;
                case 'h':
			print_usage(argv[0], stdout);
			exit(EXIT_FAILURE);
                case 'v':
			print_version(stdout);
			exit(EXIT_SUCCESS);
		default:
			wrong_usage(NULL);
			exit(EXIT_FAILURE);
		}
	}

	service = argv[optind];

	return 0;

}
