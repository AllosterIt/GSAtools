/*
 *     $Id: fragment_colour.h 1402 2013-04-16 01:19:10Z apandini $
 *     Copyright (C) 2011-2013 Alessandro Pandini
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

#ifndef FRAGMENT_COLOUR_H
#define FRAGMENT_COLOUR_H

#include "../general/safe.h"

typedef struct {
  float r;
  float g;
  float b;
} colour_rgb;

typedef struct {
	char *alphabetName;
	int nFragment;
	colour_rgb *colour;
} SA_colour;

/*____________________________________________________________________________*/
/* prototypes */
SA_colour *load_colour_set(int *ptr_nSet);
void free_colour_set(SA_colour *SA_colour_set);

#endif
