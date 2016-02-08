/*
 *     $Id: fragment_colour.c 1402 2013-04-16 01:19:10Z apandini $
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

#include "fragment_colour.h"
#include <string.h>

int nStructuralAlphabets = 1;
int nFragmentVector[] = {25};

SA_colour *load_colour_set(int *ptr_nSet){

	int i;
	SA_colour *SA_colour_set;

    *ptr_nSet = nStructuralAlphabets;

	SA_colour_set = safe_malloc(sizeof(SA_colour) * nStructuralAlphabets);

	/* M32K25 colours */
	SA_colour_set[0].alphabetName = safe_malloc(sizeof(char) * (strlen("M32K25") + 1));
	SA_colour_set[0].colour = safe_malloc(sizeof(colour_rgb) * nFragmentVector[0]);

	strcpy(SA_colour_set[0].alphabetName, "M32K25");
	SA_colour_set[0].alphabetName[strlen("M32K25")] = '\0';
	SA_colour_set[0].nFragment =  nFragmentVector[0];

	/* A [0.00, 0.00, 1.00 ] */
	SA_colour_set[0].colour[0].r = 0.00;
	SA_colour_set[0].colour[0].g = 0.00;
	SA_colour_set[0].colour[0].b = 1.00;

	/* B [0.00, 0.21, 1.00 ] */
	SA_colour_set[0].colour[1].r = 0.00;
	SA_colour_set[0].colour[1].g = 0.21;
	SA_colour_set[0].colour[1].b = 1.00;

	/* C [0.00, 0.42, 1.00 ] */
	SA_colour_set[0].colour[2].r = 0.00;
	SA_colour_set[0].colour[2].g = 0.42;
	SA_colour_set[0].colour[2].b = 1.00;

	/* D [0.00, 0.62, 1.00 ] */
	SA_colour_set[0].colour[3].r = 0.00;
	SA_colour_set[0].colour[3].g = 0.62;
	SA_colour_set[0].colour[3].b = 1.00;

	/* E [0.00, 0.83, 1.00 ] */
	SA_colour_set[0].colour[4].r = 0.00;
	SA_colour_set[0].colour[4].g = 0.83;
	SA_colour_set[0].colour[4].b = 1.00;

	/* F [0.02, 1.00, 0.96 ] */
	SA_colour_set[0].colour[5].r = 0.02;
	SA_colour_set[0].colour[5].g = 1.00;
	SA_colour_set[0].colour[5].b = 0.96;

	/* G [0.12, 1.00, 0.75 ] */
	SA_colour_set[0].colour[6].r = 0.12;
	SA_colour_set[0].colour[6].g = 1.00;
	SA_colour_set[0].colour[6].b = 0.75;

	/* H [0.23, 1.00, 0.54 ] */
	SA_colour_set[0].colour[7].r = 0.23;
	SA_colour_set[0].colour[7].g = 1.00;
	SA_colour_set[0].colour[7].b = 0.54;

	/* I [0.33, 1.00, 0.33 ] */
	SA_colour_set[0].colour[8].r = 0.33;
	SA_colour_set[0].colour[8].g = 1.00;
	SA_colour_set[0].colour[8].b = 0.33;

	/* J [0.44, 1.00, 0.12 ] */
	SA_colour_set[0].colour[9].r = 0.44;
	SA_colour_set[0].colour[9].g = 1.00;
	SA_colour_set[0].colour[9].b = 0.12;

	/* K [0.54, 1.00, 0.00 ] */
	SA_colour_set[0].colour[10].r = 0.54;
	SA_colour_set[0].colour[10].g = 1.00;
	SA_colour_set[0].colour[10].b = 0.00;

	/* L [0.64, 1.00, 0.00 ] */
	SA_colour_set[0].colour[11].r = 0.64;
	SA_colour_set[0].colour[11].g = 1.00;
	SA_colour_set[0].colour[11].b = 0.00;

	/* M [0.75, 1.00, 0.00 ] */
	SA_colour_set[0].colour[12].r = 0.75;
	SA_colour_set[0].colour[12].g = 1.00;
	SA_colour_set[0].colour[12].b = 0.00;

	/* N [0.85, 1.00, 0.00 ] */
	SA_colour_set[0].colour[13].r = 0.85;
	SA_colour_set[0].colour[13].g = 1.00;
	SA_colour_set[0].colour[13].b = 0.00;

	/* O [0.96, 1.00, 0.00 ] */
	SA_colour_set[0].colour[14].r = 0.96;
	SA_colour_set[0].colour[14].g = 1.00;
	SA_colour_set[0].colour[14].b = 0.00;

	/* P [1.00, 0.95, 0.00 ] */
	SA_colour_set[0].colour[15].r = 1.00;
	SA_colour_set[0].colour[15].g = 0.95;
	SA_colour_set[0].colour[15].b = 0.00;

	/* Q [1.00, 0.88, 0.00 ] */
	SA_colour_set[0].colour[16].r = 1.00;
	SA_colour_set[0].colour[16].g = 0.88;
	SA_colour_set[0].colour[16].b = 0.00;

	/* R [1.00, 0.81, 0.00 ] */
	SA_colour_set[0].colour[17].r = 1.00;
	SA_colour_set[0].colour[17].g = 0.81;
	SA_colour_set[0].colour[17].b = 0.00;

	/* S [1.00, 0.73, 0.00 ] */
	SA_colour_set[0].colour[18].r = 1.00;
	SA_colour_set[0].colour[18].g = 0.73;
	SA_colour_set[0].colour[18].b = 0.00;

	/* T [1.00, 0.66, 0.00 ] */
	SA_colour_set[0].colour[19].r = 1.00;
	SA_colour_set[0].colour[19].g = 0.66;
	SA_colour_set[0].colour[19].b = 0.00;

	/* U [1.00, 0.54, 0.00 ] */
	SA_colour_set[0].colour[20].r = 1.00;
	SA_colour_set[0].colour[20].g = 0.54;
	SA_colour_set[0].colour[20].b = 0.00;

	/* V [1.00, 0.40, 0.00 ] */
	SA_colour_set[0].colour[21].r = 1.00;
	SA_colour_set[0].colour[21].g = 0.40;
	SA_colour_set[0].colour[21].b = 0.00;

	/* W [1.00, 0.27, 0.00 ] */
	SA_colour_set[0].colour[22].r = 1.00;
	SA_colour_set[0].colour[22].g = 0.27;
	SA_colour_set[0].colour[22].b = 0.00;

	/* X [1.00, 0.13, 0.00 ] */
	SA_colour_set[0].colour[23].r = 1.00;
	SA_colour_set[0].colour[23].g = 0.13;
	SA_colour_set[0].colour[23].b = 0.00;

	/* Y [1.00, 0.00, 0.00 ] */
	SA_colour_set[0].colour[24].r = 1.00;
	SA_colour_set[0].colour[24].g = 0.00;
	SA_colour_set[0].colour[24].b = 0.00;

	return SA_colour_set;

}

void free_colour_set(SA_colour *SA_colour_set){

	int i;

	for(i = 0; i < nStructuralAlphabets; ++i){
		free(SA_colour_set[i].alphabetName);
		free(SA_colour_set[i].colour);
	}
	free(SA_colour_set);

}
