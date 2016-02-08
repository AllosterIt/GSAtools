/*
 *     $Id: mergesort.h 1398 2013-04-16 01:11:42Z apandini $
 *     Copyright (C) 2008-2013 Alessandro Pandini
 *
 *     This file is part of GSATools.
 *
 *     GSATools is free software: you can redistribute it and/or modify
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

#ifndef MERGESORT_H
#define MERGESORT_H

#include <stdio.h>
#include <string.h>

#include "safe.h"

/*____________________________________________________________________________*/
/* prototypes */
void MergeSort(void *array, size_t size, size_t esize, int (*compare) (const void *key1, const void *key2));

#endif

