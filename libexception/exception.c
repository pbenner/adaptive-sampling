/* -- libexception_handle.c v.0.6 --
 *
 * Copyright (C) 2002, 2003, 2004, 2005 Philipp Benner
 *
 * This file is part of UpdateDD - http://updatedd.philipp-benner.de.
 *
 * UpdateDD is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * UpdateDD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with UpdateDD; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define _BSD_SOURCE 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#ifdef HAVE_SYSLOG_H
#include <syslog.h>
#endif /* HAVE_SYSLOG_H */
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif /* HAVE_NETDB_H */
#ifdef HAVE_ERRNO_H
#include <errno.h>
#endif /* HAVE_ERRNO_H */

#define GETMODE(mode, exit_code, err_type)		\
	exit_code = (mode>>2);	    /* ...xxxxx00 */	\
	err_type  = (mode&0x3);     /* ...00000xx */	\

#define STANDARD	0
#define PERROR		1
#define HERROR		2

int verbose = 0;

static char *
create_buffer(int err_type, const char *msg)
{
	const char *ptr;
	char *buffer;
	int buflen = strlen(msg)+1;

	switch(err_type) {
	case PERROR:
		ptr = strerror(errno);
		break;
	case HERROR:
#ifdef HAVE_HSTRERROR
		ptr = hstrerror(h_errno);
#else
                ptr = "";
#endif /* HAVE_HSTRERROR */
		break;
	default:
		exit(EXIT_FAILURE);
	}

	buflen += strlen(ptr)+2;
	if((buffer = (char *)malloc(buflen)) == NULL) {
		perror("malloc() failed");
		exit(EXIT_FAILURE);
	}
	(void)sprintf(buffer, "%s: %s", msg, ptr);
	*(buffer+buflen-1) = '\0';

	return buffer;
}

/*********************** STD WARNING ******************/

void
vstd_warn(int mode, const char *msg, va_list az)
{
	int err_type, exit_code;

	GETMODE(mode, exit_code, err_type);

	if(err_type) {
		char *buffer;
		buffer = create_buffer(err_type, msg);
		(void)vfprintf(stderr, buffer, az);
		free(buffer);
	} else {
		(void)vfprintf(stderr, msg, az);
	}
	(void)fprintf(stderr, "\n");
        (void)fflush(stderr);

	return;
}

void
std_warn(int mode, const char *msg, ...)
{
	va_list az;

	va_start(az, msg);
	vstd_warn(mode, msg, az);
	va_end(az);
	return;
}

/*********************** STD NOTICE *********************/

void
vstd_notice(int mode, const char *msg, va_list az)
{
	int err_type, exit_code;

	GETMODE(mode, exit_code, err_type);

	if(err_type) {
		char *buffer;
		buffer = create_buffer(err_type, msg);
		(void)vprintf(buffer, az);
		free(buffer);
	} else {
		(void)vprintf(msg, az);
	}
	(void)printf("\n");
        (void)fflush(stdout);

	return;
}


void
std_notice(int mode, const char *msg, ...)
{
	va_list az;

	va_start(az, msg);
	vstd_notice(mode, msg, az);
	va_end(az);
	return;
}

/*********************** LOG NOTICE ******************/

void
vlog_notice(int mode, const char *msg, va_list az)
{
	int err_type, exit_code;

	GETMODE(mode, exit_code, err_type);

	if(err_type) {
		char *buffer;
		buffer = create_buffer(err_type, msg);
#ifdef HAVE_VSYSLOG
		vsyslog(LOG_NOTICE, buffer, az);
#else
		(void)vprintf(buffer, az);
#endif /* HAVE_VSYSLOG */
		free(buffer);
	} else {
#ifdef HAVE_VSYSLOG
		vsyslog(LOG_NOTICE, msg, az);
#else
		(void)vprintf(msg, az);
#endif /* HAVE_VSYSLOG */
	}

	return;
}

void
log_notice(int mode, const char *msg, ...)
{
	va_list az;

	va_start(az, msg);
	vlog_notice(mode, msg, az);
	va_end(az);
	return;
}

/*********************** STD ERROR ******************/

void
vstd_err(int mode, const char *msg, va_list az)
{

	int err_type, exit_code;

	GETMODE(mode, exit_code, err_type);

	if(err_type) {
		char *buffer;
		buffer = create_buffer(err_type, msg);
		(void)vfprintf(stderr, buffer, az);
		free(buffer);
	} else {
		(void)vfprintf(stderr, msg, az);
	}
	(void)fprintf(stderr, "\n");
        (void)fflush(stderr);

	exit(exit_code);

	return;
}

void
std_err(int mode, const char *msg, ...)
{
	va_list az;

	va_start(az, msg);
	vstd_err(mode, msg, az);
	va_end(az);
	return;
}

/*********************** S WARN ******************/

void
vs_warn(char *buf, size_t bufsize, int mode, const char *msg, va_list az)
{

	int err_type, exit_code;

	GETMODE(mode, exit_code, err_type);

	if(err_type) {
		char *buffer;
		buffer = create_buffer(err_type, msg);
		(void)vsnprintf(buf, bufsize, buffer, az);
		free(buffer);
	} else {
		(void)vsnprintf(buf, bufsize, msg, az);
	}

	return;
}

void
s_warn(char *buf, size_t bufsize, int mode, const char *msg, ...)
{
	va_list az;

	va_start(az, msg);
	vs_warn(buf, bufsize, mode, msg, az);
	va_end(az);
	return;
}

/*********************** LOG ERROR ******************/

void
vlog_err(int mode, const char *msg, va_list az)
{
	int err_type, exit_code;

	GETMODE(mode, exit_code, err_type);

	if(err_type) {
		char *buffer;
		buffer = create_buffer(err_type, msg);
#ifdef HAVE_VSYSLOG
		vsyslog(LOG_ERR, buffer, az);
#else
		(void)vfprintf(stderr, buffer, az);
#endif /* HAVE_VSYSLOG */
		free(buffer);
	} else {
#ifdef HAVE_VSYSLOG
		vsyslog(LOG_ERR, msg, az);
#else
		(void)vfprintf(stderr, msg, az);
#endif /* HAVE_VSYSLOG */
	}

	exit(exit_code);

	return;
}

void
log_err(int mode, const char *msg, ...)
{
	va_list az;

	va_start(az, msg);
	vlog_err(mode, msg, az);
	va_end(az);
	return;
}

/*********************** LOG WARNING ******************/

void
vlog_warn(int mode, const char *msg, va_list az)
{
	int err_type, exit_code;

	GETMODE(mode, exit_code, err_type);

	if(err_type) {
		char *buffer;
		buffer = create_buffer(err_type, msg);
#ifdef HAVE_VSYSLOG
		vsyslog(LOG_WARNING, buffer, az);
#else
		(void)vfprintf(stderr, buffer, az);
#endif /* HAVE_VSYSLOG */
		free(buffer);
	} else {
#ifdef HAVE_VSYSLOG
		vsyslog(LOG_WARNING, msg, az);
#else
		(void)vfprintf(stderr, msg, az);
#endif /* HAVE_VSYSLOG */
	}

	return;
}

void
log_warn(int mode, const char *msg, ...)
{
	va_list az;

	va_start(az, msg);
	vlog_warn(mode, msg, az);
	va_end(az);
	return;
}
