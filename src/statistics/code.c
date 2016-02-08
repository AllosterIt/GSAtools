/*
 *     $Id: code.c 1401 2013-04-16 01:17:37Z apandini $
 *     Copyright (C) 2010-2013 Alessandro Pandini
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

#include "code.h"

/*____________________________________________________________________________*/
/** reset probabilities  */
void reset_probabilities(Set *codeSet) {

    int i; /* index */

    for(i = 0; i < codeSet->nElements; ++i )
        codeSet->element[i].prob = 0.0;

}

/*____________________________________________________________________________*/
/** record probabilities code from string */
void record_probabilities_from_string(char *string, Set *codeSet) {

    int i,j; /* index */
    int nchar = 0;

    nchar = strlen(string);

    for(i = 0; i < nchar; ++i)
        for(j = 0; j < codeSet->nElements; ++j )
            if (string[i] == codeSet->element[j].code)
                codeSet->element[j].prob += (1.0 / nchar);

}

/*____________________________________________________________________________*/
/** generate code from string */
void extract_code_from_string(char *string, Set *codeSet) {

    int i,j; /* index */
    int allocated = 5;
    int added = 0;
    int nchar = 0;
    int (*fcmp)() = &compare_elements;

    nchar = strlen(string);

    codeSet->element = safe_malloc(allocated * sizeof(Element));
    codeSet->nElements = 0;

    for(i = 0; i < nchar; ++i) {

        for(j = 0; j < codeSet->nElements; ++j ) {
            if (string[i] == codeSet->element[j].code) {
                codeSet->element[j].prob += (1.0 / nchar);
                added = 1;
            }
        }

        if (added == 0) {
            codeSet->element[codeSet->nElements].code = string[i];
            codeSet->element[codeSet->nElements].prob = (1.0 / nchar);
            codeSet->nElements ++;

            if (codeSet->nElements == allocated) {
                allocated += 5;
                codeSet->element = safe_realloc(codeSet->element,
                                                allocated * sizeof(Element));
            }
        }

        added = 0;
    }

    MergeSort(codeSet->element, codeSet->nElements, sizeof(Element), fcmp);
}

/*____________________________________________________________________________*/
/** initialize code from string */
void initialize_code_from_string(char *string, Set *codeSet) {

    int i,j; /* index */
    int allocated = 1;
    int added = 0;
    int nchar = 0;
    int (*fcmp)() = &compare_elements;

    nchar = strlen(string);

    for(i = 0; i < nchar; ++i) {
        for(j = 0; j < codeSet->nElements; ++j ) {
            if (string[i] == codeSet->element[j].code) {
                fprintf(stderr,
                        "Warning: duplicated character %c in input string.\n",
                        string[i]);
                added = 1;
            }
        }

        if (added == 0) {
            codeSet->element[codeSet->nElements].code = string[i];
            codeSet->element[codeSet->nElements].prob = 0.0;
            codeSet->nElements ++;

            if (codeSet->nElements == allocated) {
                allocated += 5;
                codeSet->element = safe_realloc(codeSet->element,
                                                allocated * sizeof(Element));
            }
        }

        added = 0;
    }

    MergeSort(codeSet->element, codeSet->nElements, sizeof(Element), fcmp);
}
