/*
 *     $Id: safe.c 1398 2013-04-16 01:11:42Z apandini $
 *     Copyright (C) 2007-2013 Alessandro Pandini and Jens Kleinjung
 *
 *     This file is part of GSATools.
 *
 *     This is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     GSATools is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with GSATools.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "safe.h"

/*___________________________________________________________________________*/
/** safe file opening */
FILE *safe_open(const char *name, const char *mode)
{
    FILE *file = fopen(name, mode);
    if (file) {
		return file;
	} else {
		fprintf(stderr, "Error: Failed accessing file '%s'\n", name);
		exit(1);
	}
}

/*___________________________________________________________________________*/
/** safe memory allocation */
void *check_non_null(void *ptr)
{
    assert (ptr != 0);
    return ptr;
}

void *safe_malloc(size_t size)
{
    return check_non_null(malloc(size));
}

void *safe_realloc(void *ptr, size_t size)
{
    return check_non_null(realloc(ptr, size));
}

