/*
 *     $Id: vec2mat.h 1402 2013-04-16 01:19:10Z apandini $
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

#if !defined VEC2MAT_H
#define VEC2MAT_H

#include "pdb_structure.h"

/*____________________________________________________________________________*/
/* structures */

/*____________________________________________________________________________*/
/* prototypes */
void xyzVec2Mat(Vec *vec, double **matrix, int cols, int shift);
void Str2Mat(Str *str, double **matrix, int cols, int shift);

#endif
