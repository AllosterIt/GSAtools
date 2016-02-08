/*
 *     $Id: transform_segment.h 1402 2013-04-16 01:19:10Z apandini $
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

#ifndef TRANSFORM_STRUCTURE_H
#define TRANSFORM_STRUCTURE_H

#include <gsl/gsl_math.h>

#include "kabsch.h"
#include "pdb_structure.h"
#include "../general/safe.h"

/*____________________________________________________________________________*/
/* prototypes */
float coord_rmsd(Vec *v1, Vec *v2);
void kabsch_free(gsl_matrix *U, gsl_vector *t, gsl_matrix *X, gsl_matrix *Y);
void kabsch_alloc(gsl_matrix **ptr_U, gsl_vector **ptr_t, gsl_matrix **ptr_X, gsl_matrix **ptr_Y, int npoints);
void define_points(gsl_matrix *X, gsl_matrix *Y, Vec *coords_A, Vec *coords_B, int npoints); 
float sigma_square_average(gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *coords_C, int npoints, int pos, int lStruct);
float transform_coordinates(gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *coords_C, int npoints);
float superimpose_segment(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *coords_C, int npoints);
float grow_chain(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *refcoords, int npoints, int mpoints, int overlap);
float fitgrow_chain(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *refcoords, int npoints, int mpoints, int overlap);
float extend_segment(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *coords_D, Vec *refcoords, int npoints, int mpoints, int overlap);

#endif
