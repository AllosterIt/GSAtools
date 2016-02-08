/*
 *     $Id: pdb_structure.h 1402 2013-04-16 01:19:10Z apandini $
 *     Copyright (C) 2004-2008 Jens Kleinjung
 *     Copyright (C) 2008-2013 Alessandro Pandini
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

#ifndef PDB_STRUCTURE_H
#define PDB_STRUCTURE_H

/*____________________________________________________________________________*/
/* structures */

typedef struct
{
   float x, y, z;
   /*float tx, ty, tz;*/ /* disabled since 06.01.09 */
} Vec;

/* sequence */
typedef struct  
{
    char *name; /* sequence name */
    char *res; /* array of residues = sequence */
    int length; /* length of sequence */
} Seq;

/* atom : definition of PDB atom format, numbers indicate columns */
typedef struct
{
	/* PDB data */
	char recordName[8]; /* Record type; 1 -  6*/
	int atomNumber; /* Atom serial number;  7 - 11 */
	char atomName[8]; /* Atom name; 13 - 16 */
	char alternativeLocation[2]; /* Alternate location indicator; 17 */
	char residueName[4]; /* Residue name; 18 - 20 */
	char chainIdentifier[2]; /* Chain identifier; 22 */
	int residueNumber; /* Residue sequence number; 23 - 26 */
	Vec pos; /* position vector */
	char description[32]; /* everything before coordinates */
	Vec tpos; /* transformed position vector */
} Atom;

/* molecular structure */
typedef struct
{
	Atom *atom; /* array of selected (CA) atoms constituting structure */
	int natom; /* number of selected (CA) atoms */
	Seq sequence; /* sequence of structure */
} Str;

#endif
