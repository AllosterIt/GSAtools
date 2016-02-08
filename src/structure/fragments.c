/*
 *     $Id: fragments.c 1402 2013-04-16 01:19:10Z apandini $
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

#include "fragments.h"
#include "fragments_coords.h"
#include "getfragments.h"

int initialise_fragment_set(FragmentSet *selected_set, char *setname){

    int i, j, k;

    /*____________________________________________________________________________*/
    /* initialise fragment sets with constant values for the available alphabets */

    ConstantFragmentSet *constant_set;

    for (constant_set = constant_fragment_sets; constant_set < constant_fragment_sets + sizeof constant_fragment_sets / sizeof *constant_fragment_sets; constant_set ++) 
        if (strcmp(constant_set->setname, setname) == 0){
            selected_set->setname = constant_set->setname;
            selected_set->lFragment = constant_set->lFragment;
            selected_set->nFragment = constant_set->nFragment;
            selected_set->codeOrder = constant_set->codeOrder;
            selected_set->coord_values = alloc_float_matrix3D(selected_set->coord_values, selected_set->nFragment, selected_set->lFragment, 3);
            for (i = 0; i < selected_set->nFragment; ++ i)
                for (j = 0; j < selected_set->lFragment; ++ j)
                    for (k = 0; k < 3; ++ k)
                        selected_set->coord_values[i][j][k] = constant_set->coord_values[i][j][k];
            return 0;
        }

    fprintf(stderr, "\nExiting: Fragment set '%s' not implemented!\nAvailable sets:\n", setname);
    for (constant_set = constant_fragment_sets; constant_set < constant_fragment_sets + sizeof constant_fragment_sets / sizeof *constant_fragment_sets; constant_set ++)
        fprintf(stderr, "\t%s\n", constant_set->setname);

	return 1;
}

Str *load_fragment_data(char *pdbCodeFileName, char *userCodeString, char *stralphabetName, FragmentSet *fragment_set, int pdbcode){

    unsigned int i, j; /* indices */
    Str *fragment_str; /** structure array for alphabet fragments */
    FILE *pdbCodeFile;

    /*____________________________________________________________________________*/
    /** initialise fragment sets */
    if (pdbcode) {
        pdbCodeFile = safe_open(pdbCodeFileName, "r");
        read_n_fragments(pdbCodeFile, fragment_set, userCodeString);
        read_fragments(pdbCodeFile, fragment_set);
        fclose(pdbCodeFile);
    }else{
        if (initialise_fragment_set(fragment_set, stralphabetName) != 0)
			exit(1);
    }

    /*____________________________________________________________________________*/
    /** get structure alphabet (stralphabet) */
    fragment_str = safe_malloc(fragment_set->nFragment * sizeof(Str));
    for (i = 0; i < fragment_set->nFragment; ++ i)
    {
        fragment_str[i].atom = safe_malloc(fragment_set->lFragment * sizeof(Atom));
        for (j =0; j < fragment_set->lFragment; ++ j)
        {
            fragment_str[i].atom[j].pos.x = fragment_set->coord_values[i][j][0]; 
            fragment_str[i].atom[j].pos.y = fragment_set->coord_values[i][j][1]; 
            fragment_str[i].atom[j].pos.z = fragment_set->coord_values[i][j][2]; 
        }
    }

    return fragment_str;
}

float frag_overlap_rmsd(Str *fragment1, Str *fragment2, int lFragment, int overlap){

    float rmsd = -1; /* rmsd value of fragment superpositioning */
    double **matFrag1, **matFrag2; /* matrix required by quadRMSD */
    double coeff[3]; /* vector required by quadRMSD */

    matFrag1 = MatInit(3, overlap);
    matFrag2 = MatInit(3, overlap);

    /** calculate RMSD with Quaternions */
    Str2Mat(fragment1, matFrag1, overlap, (lFragment - overlap));
    Str2Mat(fragment2, matFrag2, overlap, 0);
    rmsd = QCP_rmsd(matFrag1, matFrag2, overlap, coeff);

    MatDestroy(matFrag1);
    MatDestroy(matFrag2);

    return(rmsd);
}

float local_rmsd(Vec *subStructure, Str *fragment, int lSubStr, int overlap){

    float rmsd = -1; /* rmsd value of fragment superpositioning */
    double **matSub, **matFrag; /* matrix required by quadRMSD */
    double coeff[3]; /* vector required by quadRMSD */

    matSub = MatInit(3, lSubStr);
    matFrag = MatInit(3, lSubStr);

    /** calculate RMSD with Quaternions */
    xyzVec2Mat(subStructure, matSub, lSubStr, (lSubStr - overlap));
    Str2Mat(fragment, matFrag, lSubStr, 0);
    rmsd = QCP_rmsd(matSub, matFrag, lSubStr, coeff);

    MatDestroy(matSub);
    MatDestroy(matFrag);

    return(rmsd);
}

int get_best_fragment(Vec *subStructure, FragmentSet *fragment_set, Str *fragment_str, float *ptr_rmsd){

    unsigned int j; /* indices */

    float rmsd = -1; /* rmsd value of fragment superpositioning */
    float minRmsd = 100.; /* lowest rmsd value */
    int bestFragment_index = -1; /* stralphabet fragment with lowest rmsd value */
    double **matSub, **matFrag; /* matrix required by quadRMSD */
    double coeff[3]; /* vector required by quadRMSD */

    matSub = MatInit(3, fragment_set->lFragment);
    matFrag = MatInit(3, fragment_set->lFragment);

    /** against alphabet fragments of length 'lFragment' */
    for (j = 0; j < fragment_set->nFragment; ++ j )
    {
        /** calculate RMSD with Quaternions */
        xyzVec2Mat(subStructure, matSub, fragment_set->lFragment, 0);
        Str2Mat(&(fragment_str[j]), matFrag, fragment_set->lFragment, 0);
        rmsd = QCP_rmsd(matSub, matFrag, fragment_set->lFragment, coeff);
        /** select stralphabet fragment with lowest rmsd to PDB fragment */
        if ((minRmsd = GSL_MIN(rmsd, minRmsd)) == rmsd)
            bestFragment_index = j;
    }

    *ptr_rmsd = minRmsd;

    MatDestroy(matSub);
    MatDestroy(matFrag);

    return bestFragment_index;
}

void free_fragment_sets(FragmentSet *fragment_set)
{
    free_float_matrix3D(fragment_set->coord_values, fragment_set->nFragment, fragment_set->lFragment);
    free(fragment_set);
}
