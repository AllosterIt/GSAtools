/*
 *     $Id: encode.c 1400 2013-04-16 01:15:21Z apandini $
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

#include "encode.h"

extern float CADistCutoff;

/*____________________________________________________________________________*/
/** compare reconstructed structures */
int RecStrCmp(ReconStructure *recStrA, ReconStructure *recStrB) {

    float rmsd_A = (float) recStrA->rmsd;
    float rmsd_B = (float) recStrB->rmsd;

    return ((rmsd_A > rmsd_B) ? 1 : (rmsd_A < rmsd_B) ? -1 : 0);
}

/*____________________________________________________________________________*/
/** copy reconstructed structures */
void RecStrCpy(ReconStructure *recStrDest, ReconStructure *recStrSource) {

    int i; /* index */

    /* copy encoded string */
    strcpy(recStrDest->encodedString, recStrSource->encodedString);
    /* copy RMSD value */
    recStrDest->rmsd = recStrSource->rmsd;
    /* copy coordinates */
    for(i = 0; i < recStrSource->natom; ++ i) {
        recStrDest->coord[i].x = recStrSource->coord[i].x;
        recStrDest->coord[i].y = recStrSource->coord[i].y;
        recStrDest->coord[i].z = recStrSource->coord[i].z;
    }
}

/*____________________________________________________________________________*/
/** free memory from reconstructed structures */
void free_ReconStructure(ReconStructure *reconStr, int nRecStr) {

    int i; /* counter */

    for(i = 0; i < nRecStr; ++ i) {
        free(reconStr[i].encodedString);
        free(reconStr[i].coord);
    }

    free(reconStr);
}

/*____________________________________________________________________________*/
/** allocate memory for reconstructed structures */
ReconStructure *initialise_ReconStructure(int lFragment, int natom,
        int nRecStr) {

    int i; /* counter */
    int nChar; /* n character in encoded string */
    ReconStructure *recStr; /* array of reconstructed structures */

    nChar = natom - lFragment + 1;

    recStr = safe_malloc(nRecStr * sizeof(ReconStructure));

    for(i = 0; i < nRecStr; ++ i) {
        recStr[i].natom = natom;
        recStr[i].encodedString = safe_malloc((nChar + 1) * sizeof(char));
        strcpy(recStr[i].encodedString, "");
        recStr[i].coord = safe_malloc(natom * sizeof(Vec));
        recStr[i].rmsd = FLT_MAX;
    }

    return(recStr);
}

/*____________________________________________________________________________*/
/** encode string by global fit */
void globalfit_encode(ReconStructure *recStr, ReconStructure *tempRecStr,
                      Vec *refcoord, FragmentSet *fragment_set,
                      Str *fragment_str, int pos, int heapsize) {

    int i, j, k; /* counter */
    int overlap = 3; /* size of overlap */
    Vec *fragCoord; /* fragment coordinates */
    int heapCount = 0; /* counter for structure in the heap */
    int (*fcmp)() = &RecStrCmp;
    int idx;

    /** symmetry operators for structural superpositioning */
    gsl_matrix *U = 0; /* rotation matrix */
    gsl_vector *t = 0; /* translation vector */
    gsl_matrix *X = 0; /* structure fragment matrix query */
    gsl_matrix *Y = 0; /* structure fragment matrix template */

    /* allocate memory for fragment coordinates */
    fragCoord = safe_malloc(fragment_set->lFragment * sizeof(Vec));

    /* at the beginning of the reconstruction the head of the heap is
     * filled with the fragments */
    if (pos == 0) {
        /* allocate memory for kabsch matrices */
        kabsch_alloc(&U, &t, &X, &Y, fragment_set->lFragment);
        for(i = 0; i < fragment_set->nFragment; ++ i) {
            for (j = 0; j < fragment_set->lFragment; j ++) {
                fragCoord[j].x = fragment_str[i].atom[j].pos.x;
                fragCoord[j].y = fragment_str[i].atom[j].pos.y;
                fragCoord[j].z = fragment_str[i].atom[j].pos.z;
            }
            recStr[i].rmsd = superimpose_segment(X, Y, U, t, &refcoord[pos],
                                                 fragCoord,
                                                 &recStr[i].coord[pos],
                                                 fragment_set->lFragment);
            recStr[i].encodedString[pos] = fragment_set->codeOrder[i];
            recStr[i].encodedString[pos + 1] = '\0';
        }
    } else {
        /* allocate memory for kabsch matrices */
        kabsch_alloc(&U, &t, &X, &Y, overlap);
        for(k = 0; k < heapsize; ++ k) {
            if (strcmp(recStr[k].encodedString, "") != 0) {
                for(i = 0; i < fragment_set->nFragment; ++ i) {
                    idx = k * fragment_set->nFragment + i;
                    for (j = 0; j < fragment_set->lFragment; j ++) {
                        fragCoord[j].x = fragment_str[i].atom[j].pos.x;
                        fragCoord[j].y = fragment_str[i].atom[j].pos.y;
                        fragCoord[j].z = fragment_str[i].atom[j].pos.z;
                    }
                    tempRecStr[idx].rmsd = extend_segment(X, Y, U, t,
                                                          recStr[k].coord,
                                                          fragCoord,
                                                          tempRecStr[idx].coord,
                                                          refcoord,
                                                          pos + overlap,
                                                          fragment_set->lFragment,
                                                          overlap);
                    strcpy(tempRecStr[idx].encodedString,
                           recStr[k].encodedString);
                    tempRecStr[idx].encodedString[pos] =
                        fragment_set->codeOrder[i];
                    tempRecStr[idx].encodedString[pos + 1] = '\0';
                    heapCount += 1;
                }
            }
        }
        /* check if tentative structures are more than what the heap can hold */
        if (heapCount < heapsize) {
            /* sort temporary reconstructed structures array by RMSD */
            MergeSort(tempRecStr, heapCount, sizeof(ReconStructure), fcmp);
            /* copy from temporary reconstructed structures array */
            for(k = 0; k < heapsize; ++ k)
                if (strcmp(tempRecStr[k].encodedString, "") != 0)
                    RecStrCpy(&recStr[k], &tempRecStr[k]);
        } else {
            /* sort temporary reconstructed structures array by RMSD */
            MergeSort(tempRecStr, heapCount, sizeof(ReconStructure), fcmp);
            /* copy best RMSD structures to reconstructed structures array */
            for(k = 0; k < heapsize; ++ k)
                RecStrCpy(&recStr[k], &tempRecStr[k]);
        }
    }

    /* free memory from kabsch matrices */
    kabsch_free(U, t, X, Y);

    /* free memory from fragment coordinates */
    free(fragCoord);
}

/*____________________________________________________________________________*/
/** encode string by local fit */
void localfit_encode(ReconStructure *recStr, Vec *refcoord,
                     FragmentSet *fragment_set, Str *fragment_str,
                     float *localfit_rmsd_array) {

    unsigned int i; /* indices */
    float rmsd; /* local RMSD */
    float rmsdSum = 0; /* sum of local RMSDs */

    /** structure superpositioning parameters*/
    unsigned int shift = 1; /* peptide shift */
    int bestFragment = -1; /* stralphabet fragment with lowest rmsd value */
	float CADist = 0; /* C-alpha - C-alpha distance */

    /*________________________________________________________________________*/
    /** symmetry operators for structural superpositioning */
    gsl_matrix *U = 0; /* rotation matrix */
    gsl_vector *t = 0; /* translation vector */
    gsl_matrix *X = 0; /* structure fragment matrix query */
    gsl_matrix *Y = 0; /* structure fragment matrix template */

    /*________________________________________________________________________*/
    /* allocate memory for kabsch superpositioning */
    kabsch_alloc(&U, &t, &X, &Y, fragment_set->lFragment);

    /*________________________________________________________________________*/
    /** superimpose PDB structure with alphabet fragments */
    /** scan fragments of length 'lFragment' of PDB structure */
    for (i = 0; i <= recStr->natom - fragment_set->lFragment; i += shift)
    {
        CADist = coord_rmsd(&(refcoord[i+2]), &(refcoord[i+3]));
        if (CADist <= CADistCutoff){
			bestFragment = get_best_fragment(&refcoord[i], fragment_set,
                                         fragment_str, &rmsd);

			/** translate best fragment to stralphabet character */
			recStr->encodedString[i] = fragment_set->codeOrder[bestFragment];

			/** update local RMSD array */
			localfit_rmsd_array[i] = rmsd;

			/** add local RMSD to sum */
			rmsdSum += rmsd;
		} else {
			/*fprintf(stderr, "Chain break: Not encoding positions %d, %d, %d\n", i+1, i+2, i+3);*/ 
			recStr->encodedString[i] = ' ';
			recStr->encodedString[i+1] = ' ';
			recStr->encodedString[i+2] = ' ';
			i += 2;
		}
    }

    /*________________________________________________________________________*/
    /* free memory from kabsch superpositioning matrices */
    kabsch_free(U, t, X, Y);

    /*________________________________________________________________________*/
    /* properly terminate encoded string and print it if not silent */
    recStr->encodedString[i] = '\0';

    /*________________________________________________________________________*/
    /* update average rmsd value */
    recStr->rmsd = (rmsdSum / (recStr->natom - fragment_set->lFragment + 1) );
}
