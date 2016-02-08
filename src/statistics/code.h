/*
 *     $Id: code.h 1401 2013-04-16 01:17:37Z apandini $
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

#ifndef CODE_H_
#define CODE_H_

#include <string.h>

#include "../general/mergesort.h"
#include "probability.h"

/*____________________________________________________________________________*/
/* data structures */

/*____________________________________________________________________________*/
/* prototypes */

void reset_probabilities(Set *codeSet);
void record_probabilities_from_string(char *string, Set *codeSet);
void extract_code_from_string(char *string, Set *codeSet);
void initialize_code_from_string(char *string, Set *codeSet);

#endif /* CODE_H_ */
