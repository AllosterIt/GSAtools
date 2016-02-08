/*
 *     $Id: quatRMSD.h 1402 2013-04-16 01:19:10Z apandini $
 *     Copyright (c) 2005-2008 Douglas L. Theobald
 *     Copyright (c) 2008-2013 Alessandro Pandini
 *
 *     This file is part of GSATools. 
 *
 *     This code is a modification of QuatCharPoly.c from Douglas L. Theobald
 *     The original code can be found at:
 *     http://monkshood.colorado.edu/QCP
 *
 *     If you use this code in a publication, please reference: 
 *
 *        Douglas L. Theobald (2005)
 *        "Rapid calculation of RMSD using a quaternion-based characteristic
 *         polynomial."
 *        Acta Crystallographica A 61(4):478-480.
 *
 *        Horn, B. K. P. (1987).
 *        "Closed-form solution of absolute orientation using unit quaternions."
 *        J Opt Soc Am A 4(4):629Ð642.
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

#if !defined QUATRMSD_H
#define QUATRMSD_H 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*____________________________________________________________________________*/
/* structures */

/*____________________________________________________________________________*/
/* prototypes */
void PrintCoords(const double **coords, const int len);
double **MatInit(const int rows, const int cols);
void MatDestroy(double **matrix);
double QuatCharPoly(const double **coords1, const double **coords2, const int len, double *coeff);
double QCP_rmsd(double **coords1, double **coords2, const int len, double *coeff);
void CenterCoords(double **coords, const int len);

#endif
