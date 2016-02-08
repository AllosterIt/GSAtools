/*
 *     $Id: sequence.h 1400 2013-04-16 01:15:21Z apandini $
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

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../general/safe.h"
#include "../structure/pdb_structure.h"

/*____________________________________________________________________________*/
/* structures */

typedef struct {
    Seq *sequence; /* array of sequences */
    int nSequences; /* number of elements */
} SeqSet;

/*____________________________________________________________________________*/
/** prototypes */
void read_inputFileEnsemble(FILE *inputFile, SeqSet *inputSequenceSet);
char *get_string_from_column(SeqSet *fastaSequenceSet, int column);
void add_sequence_to_set(SeqSet *localfastaSequenceSet, char *frameDesc,
                         char *encodedString, int seqLength, int fastaIndex);

#endif /* SEQUENCE_H_ */
