/*
 *     $Id: fragments.h 1402 2013-04-16 01:19:10Z apandini $
 *     Copyright (C) 2007-2013 Alessandro Pandini and Jens Kleinjung
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

#ifndef FRAGMENTS_H
#define FRAGMENTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>

#include "pdb_structure.h"
#include "quatRMSD.h"
#include "vec2mat.h"
#include "../general/matrix.h"
#include "../general/safe.h"

/*____________________________________________________________________________*/
/* structures */

/* Fragment set */

typedef struct {
    char *setname;
    int lFragment;
    int nFragment;
    char *codeOrder;
    float ***coord_values;
} FragmentSet;

/*____________________________________________________________________________*/
/* prototypes */
int initialise_fragment_set(FragmentSet *selected_set, char *setname);
Str *load_fragment_data(char *pdbCodeFileName, char *userCodeString, char *stralphabetName, FragmentSet *fragment_set, int pdbCode);
float frag_overlap_rmsd(Str *fragment1, Str *fragment2, int lFragment, int overlap);
float local_rmsd(Vec *subStructure, Str *fragment, int lSubstr, int overlap);
int get_best_fragment(Vec *subStructure, FragmentSet *fragment_set, Str *fragment_str, float *ptr_rmsd);
void free_fragment_sets(FragmentSet *fragment_sets);

#endif
