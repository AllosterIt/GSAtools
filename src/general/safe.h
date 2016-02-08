/*
 *     $Id: safe.h 1398 2013-04-16 01:11:42Z apandini $
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

#ifndef SAFE_H
#define SAFE_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/*___________________________________________________________________________*/
/* file */
FILE *safe_open(const char *name, const char *mode);

/*___________________________________________________________________________*/
/* allocation */
void *check_non_null(void *ptr);
void *safe_malloc(size_t size);
void *safe_realloc(void *ptr, size_t size);

#endif

