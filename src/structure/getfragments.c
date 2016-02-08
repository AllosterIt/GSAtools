/*
 *     $Id: getfragments.c 1402 2013-04-16 01:19:10Z apandini $
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

#include "getfragments.h"

void read_n_fragments(FILE *fragmentfile, FragmentSet *fragment_set, char *extended_codeOrder)
{
    unsigned int i,j;
    char line[80] = "";
    unsigned int nFragment = 0;
    unsigned int lFragment = 0;
    unsigned int nCA = 0;
    char atomName[8];
    char *pfgets = 0;

    /*____________________________________________________________________________*/
    /** determine number of fragment and fragment length */
    while(! feof(fragmentfile))
    {
         pfgets = fgets(line, 80, fragmentfile);

        /* increase frament count if record is MODEL */
         if(strncmp(line, "MODEL ", 6) == 0)
         {
             ++ nFragment;
             lFragment = 0;
         }

        /* skip if record is not ATOM */
        if(strncmp(line, "ATOM  ", 6) != 0)
            continue;

         /* atom name */
         for (i = 12, j = 0; i < 16; )
             atomName[j++] = line[i++];
         atomName[j] = '\0';

         /* recording only Calpha atoms */
         if (strncmp(atomName, " CA ", 4) != 0)
             continue;

         ++ lFragment;
         ++ nCA;
    }

    /* check consistency of numerb of CA, nFragment and lFragment */
    if (nCA / nFragment != lFragment)
    {
        fprintf(stderr, "Inconsistency in the fragment file provided by the user\n");
        exit(1);
    }

    assert(strlen(extended_codeOrder) >= nFragment);

    fragment_set->lFragment = lFragment;
    fragment_set->nFragment = nFragment;
    fragment_set->codeOrder = safe_malloc((nFragment + 1) * sizeof(char));
    strncpy(fragment_set->codeOrder, extended_codeOrder, nFragment);
    fragment_set->codeOrder[nFragment] = '\0';
    fragment_set->coord_values = alloc_float_matrix3D(fragment_set->coord_values, fragment_set->nFragment, fragment_set->lFragment, 3);
}

void read_fragments(FILE *fragmentfile,  FragmentSet *fragment_set)
{
    unsigned int i,j;
    char line[80] = "";
    unsigned int frag_idx = -1;
    unsigned int atom_idx = 0;
    char atomName[8];
    char *pfgets = 0;

    rewind(fragmentfile); /* rewind file to beginning */

    /* read the user defined fragments from file */
    while(! feof(fragmentfile))
    {
         pfgets = fgets(line, 80, fragmentfile);

        /* increase frament index if record is MODEL */
         if(strncmp(line, "MODEL ", 6) == 0)
         {
             ++ frag_idx;
            atom_idx = 0;
         }

        /* skip if record is not ATOM */
        if(strncmp(line, "ATOM  ", 6) != 0)
            continue;

        /* atom name */
        for (i = 12, j = 0; i < 16; )
            atomName[j++] = line[i++];
        atomName[j] = '\0';

        /* recording only Calpha atoms */
        if (strncmp(atomName, " CA ", 4) != 0)
            continue;

        /* coordinates */
        fragment_set->coord_values[frag_idx][atom_idx][0] = atof(&line[30]);
        fragment_set->coord_values[frag_idx][atom_idx][1] = atof(&line[38]);
        fragment_set->coord_values[frag_idx][atom_idx][2] = atof(&line[46]);

        ++ atom_idx;
    }
}
