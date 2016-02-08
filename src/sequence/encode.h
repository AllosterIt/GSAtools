/*
 *     $Id: encode.h 1400 2013-04-16 01:15:21Z apandini $
 *     Copyright (C) 2008-2013 Alessandro Pandini and Jens Kleinjung
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

#ifndef ENCODE_H
#define ENCODE_H

#include "../general/safe.h"
#include "../general/mergesort.h"
#include "../structure/fragments.h"
#include "../structure/pdb_structure.h"
#include "../structure/transform_segment.h"

/*____________________________________________________________________________*/
/* structures */

typedef struct {
    int natom;
    char *encodedString;
    Vec *coord;
    float rmsd;
} ReconStructure;

/*____________________________________________________________________________*/
/* prototypes */
void free_ReconStructure(ReconStructure *reconStr, int nRecStr);
ReconStructure *initialise_ReconStructure(int lFragment, int natom,
        int nRecStr);
void globalfit_encode(ReconStructure *recStr, ReconStructure *tempRecStr,
                      Vec *refcoord, FragmentSet *fragment_set,
                      Str *fragment_str, int pos, int heapsize);
void localfit_encode(ReconStructure *recStr, Vec *refcoord,
                     FragmentSet *fragment_set, Str *fragment_str,
                     float *localfit_rmsd_array);

#endif
