/*
 *     $Id: sequence.c 1400 2013-04-16 01:15:21Z apandini $
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

#include "sequence.h"

/*____________________________________________________________________________*/
/** read input sequences from file */
void read_inputFileEnsemble(FILE *inputFile, SeqSet *inputSequenceSet) {

    int i;
	const int lLine = 8192;
    char line[8192] = "";
    int allocated_seq = 64;
    int allocated_res = 64;

    /*________________________________________________________________________*/
    /* allocate memory for sequences */
    inputSequenceSet->sequence = safe_malloc(allocated_seq * sizeof(Seq));
    inputSequenceSet->nSequences = 0;

    while(fgets(line, lLine, inputFile) != NULL ) {
		/* new sequence */
		inputSequenceSet->sequence[inputSequenceSet->nSequences].res = safe_malloc(allocated_res * sizeof(char));
        inputSequenceSet->sequence[inputSequenceSet->nSequences].length = 0;
        /* sequence name */
		inputSequenceSet->sequence[inputSequenceSet->nSequences].name = safe_malloc(sizeof(char) * 3);
        strcpy(inputSequenceSet->sequence[inputSequenceSet->nSequences].name, "NA");
		/* sequence consistency */
		assert((strlen(line) < lLine) && "increase input string length limit in code");
		if (inputSequenceSet->nSequences > 0){
			assert(((strlen(line) - 1) == inputSequenceSet->sequence[inputSequenceSet->nSequences - 1].length) && "sequences must be of equal length");
                }

        /*________________________________________________________________*/
		/* sequence residues */
        for(i = 0; i < strlen(line); ++ i) {
            /* if character is not end of line include it */
            if (line[i] != '\n') {
				/* accept only characters within alphabet range */
                if (line[i] < 65 || line[i] > 122 || (line[i] >= 91 && line[i] <= 96)) { 
                    fprintf(stderr, "%s:%d: Exiting: unusable character '%c' (ASCII %d) in seq. %d, col. %d\n",
						__FILE__, __LINE__, line[i], line[i], inputSequenceSet->nSequences, i);
					exit(1);
				}

				/* assign residue */
                inputSequenceSet->sequence[inputSequenceSet->nSequences].res[inputSequenceSet->sequence[inputSequenceSet->nSequences].length] = line[i];
                ++ inputSequenceSet->sequence[inputSequenceSet->nSequences].length;

				/* allocate more memory for sequence residues if needed */
				if (inputSequenceSet->sequence[inputSequenceSet->nSequences].length == allocated_res) {
					allocated_res += 64;
					inputSequenceSet->sequence[inputSequenceSet->nSequences].res =
						safe_realloc(inputSequenceSet->sequence[inputSequenceSet->nSequences].res, allocated_res * sizeof(char));
				}
            }
        }

		/*____________________________________________________________________*/
		/* complete new sequence */
		inputSequenceSet->sequence[inputSequenceSet->nSequences].res[inputSequenceSet->sequence[inputSequenceSet->nSequences].length] = '\0';
        ++ inputSequenceSet->nSequences;

        /* allocate more memory for sequences if needed */
		if (inputSequenceSet->nSequences == allocated_seq) {
			allocated_seq += 64;
			inputSequenceSet->sequence = safe_realloc(inputSequenceSet->sequence, allocated_seq * sizeof(Seq));
		}
    }
#ifdef DEBUG
	fprintf(stderr, "number of input sequences: %d\n", inputSequenceSet->nSequences);
#endif
}

/*____________________________________________________________________________*/
/** get string from column of sequences */
char *get_string_from_column(SeqSet *fastaSequenceSet, int column) {

    int i;
    char *column_string;

    column_string = safe_malloc(sizeof(char) * (fastaSequenceSet->nSequences + 1));

    /* read characters from column */
    for(i = 0; i < fastaSequenceSet->nSequences; ++ i) {
        column_string[i] = fastaSequenceSet->sequence[i].res[column];
                assert((column_string[i] >= 65 && column_string[i] <= 122) && (column_string[i] < 91 || column_string[i] > 96) && "unusable character in alignment column");
	}
    column_string[i] = '\0';

    return(column_string);

}

/*____________________________________________________________________________*/
/** add string to sequence set */
void add_sequence_to_set(SeqSet *localfastaSequenceSet, char *frameDesc,
                         char *encodedString, int seqLength, int fastaIndex) {

    /*________________________________________________________________________*/
    /* allocate memory for new sequence */
    localfastaSequenceSet->nSequences++;
    localfastaSequenceSet->sequence =
        safe_realloc(localfastaSequenceSet->sequence,
                     sizeof(Seq) *
                     localfastaSequenceSet->nSequences);

    /*________________________________________________________________*/
    /* allocate memory and initialize res string */
    localfastaSequenceSet->sequence[fastaIndex].res =
        safe_malloc(sizeof(char) * seqLength);
    localfastaSequenceSet->sequence[fastaIndex].name =
        safe_malloc(sizeof(char) * (strlen(frameDesc) + 1));
    localfastaSequenceSet->sequence[fastaIndex].length = seqLength;

    strncpy(localfastaSequenceSet->sequence[fastaIndex].name, frameDesc,
            strlen(frameDesc));
    localfastaSequenceSet->sequence[fastaIndex].name[strlen(frameDesc)] = '\0';
    strncpy(localfastaSequenceSet->sequence[fastaIndex].res, encodedString,
            seqLength);

}
