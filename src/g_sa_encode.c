/*
 *     $Id: g_sa_encode.c 1418 2013-04-23 16:33:35Z apandini $
 *     Copyright (C) 2011-2013 Alessandro Pandini
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

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "pdbio.h"
#include "tpxio.h"
#include "xvgr.h"
#include "matio.h"

#include "general/safe.h"
#include "sequence/encode.h"
#include "sequence/sequence.h"
#include "statistics/code.h"
#include "statistics/probability.h"
#include "structure/fragments.h"
#include "structure/fragment_colour.h"
#include "structure/pdb_structure.h"
#include "structure/transform_segment.h"

/*____________________________________________________________________________*/
/* Warning cutoff for C-alpha - C-alpha distance */
/* float CADistCutoff = 4.1; */
float CADistCutoff = 6.0;

/*____________________________________________________________________________*/
/** encode structure */
void encode_structure(FragmentSet *fragment_set, Str *fragment_str,
		t_topology *top, t_trxframe *fr, char *title,
		FILE *rmsdlocalFile, FILE *localFile, FILE *fastalocalFile,
		FILE *rmsdglobalFile, FILE *globalFile, FILE *fastaglobalFile,
		FILE *pdbglobalFile, FILE *logFile, int heapsize, gmx_bool fasta,
		gmx_bool globalfit, gmx_bool verbose, gmx_bool entropy, gmx_bool xpmoutput,
		int *ptr_fastaIndex,
		SeqSet *localfastaSequenceSet, SeqSet *globalfastaSequenceSet,
		int *chainBreak){

    int i = 0; 
	int j = 0; /* indexes */
	int k = 0;
    Vec *refcoord; /* array of ref coords */
    ReconStructure *recStr; /* array of reconstructed structures */
    ReconStructure *tempRecStr = 0; /* temporary array for attempted recon */
    float *verboseRMSDvec = 0; /* RMSD vector for partial reconstruction */
    float *localfit_rmsd_array; /* local fit fragment RMSD array */
    float *globalfit_rmsd_array = 0; /* global fit CA RMSD array */
    char ***gnames = 0; /* index group names */
    t_blocka *grps = 0; /* index groups */
    int isize = 0; /* index group size */
    atom_id *index = NULL; /* atom indeces */
    int nCalpha; /* number of C-alpha */
    int NM2ANG = 10.00; /* scaling factor */
    float CADist = 0; /* C-alpha - C-alpha distance */
    char frameDesc[64]; /* frame description */
    char strtitle[256]; /* pdb title */
    rvec *xmem = NULL; /* coord vector */
    t_trxframe frout; /* output frame */
    t_atoms useatoms; /* atom info */
	int allocated = 64; /* allocated memory */

    /*________________________________________________________________________*/
    /* get C-alpha indices */
	if (globalfit){
	    snew(gnames, 1);
	    snew(grps, 1);
	    snew(grps->index, 1);
	    /* generate standard groups */
	    analyse(&(top->atoms), grps, gnames, FALSE, FALSE);
	    /* find C-alpha group */
	    while(strcmp(gnames[0][i], "C-alpha"))
	    	i++;
	    /* get C-alpha indices */
	    isize = grps->index[i+1]-grps->index[i];
	    snew(index, isize);
	    for(j = 0; j < isize; j++)
	    	index[j] = grps->a[grps->index[i]+j];
	}

    /*________________________________________________________________________*/
    /* count C-alpha */
    nCalpha = 0;
    for(i = 0; i < top->atoms.nr; ++i)
    	if (strcmp(*(top->atoms.atomname[i]), "CA") == 0)
    		nCalpha ++;
    if (globalfit && isize != nCalpha){
    	fprintf(stderr, "Error on number of C-alpha atoms\n");
    	exit(1);
    }

    /*________________________________________________________________________*/
    /* allocate memory for reference, reconstructed coords and RMSD arrays */
    refcoord = safe_malloc(nCalpha * sizeof(Vec));
    localfit_rmsd_array = safe_malloc((nCalpha - fragment_set->lFragment + 1) * sizeof(float));
    if (globalfit)
	    globalfit_rmsd_array = safe_malloc(nCalpha * sizeof(float));

    /*________________________________________________________________________*/
    /* initialise RMSD arrays */
    for (i = 0; i < nCalpha - fragment_set->lFragment; ++ i){
        localfit_rmsd_array[i] = 0.0;
    }
    if (globalfit){
	    for (i = 0; i < nCalpha; ++ i){
	        globalfit_rmsd_array[i] = 0.0;
	    }
    }

    /*________________________________________________________________________*/
    /* set refcoords vector to CA frame coords */
    i = 0;
    for (j = 0; j < top->atoms.nr; ++ j){
    	if (strcmp(*(top->atoms.atomname[j]), "CA") == 0){
	        refcoord[i].x = fr->x[j][XX] * NM2ANG;
	        refcoord[i].y = fr->x[j][YY] * NM2ANG;
	        refcoord[i].z = fr->x[j][ZZ] * NM2ANG;

	        i++;
    	}
    }
    /*________________________________________________________________________*/
    /* check for chain breaks */
	/* print this only once */
	if (*chainBreak == 0) {
		for (j = 0; j < (nCalpha - 1); ++ j){
			CADist = coord_rmsd(&(refcoord[j]), &(refcoord[j+1]));
			if (CADist > CADistCutoff){
					fprintf(stderr, "\nWarning: distance CA%d - CA%d is %6.3f A (expected 3.8 A, cutoff is %f A).\n",
																		j+1, j+2, CADist, CADistCutoff);

					fprintf(stderr, "         This indicates (1) a chain break or (2) transgression of a periodic boundary.\n");
					fprintf(stderr, "         (1) The last three residues of each chain are by definition not encoded.\n");
					fprintf(stderr, "               Remember that when mapping correlations to the initial structure,\n");
					fprintf(stderr, "               because that structure contains 3 residues more per additional chain\n");
					fprintf(stderr, "               than the correlation matrix.\n");
					fprintf(stderr, "         (2) Use the program trjconv to centre the molecule in the periodic box.\n\n");

					fprintf(stderr, "Assuming chain break: Not encoding positions %d, %d, %d.\n\n", j+1, j+2, j+3);

					++ (*chainBreak);
			}
		}
	}

    /*________________________________________________________________________*/
    /** allocate memory for reconstructed structures */
    recStr = initialise_ReconStructure(fragment_set->lFragment, nCalpha, heapsize);
    if (globalfit)
	    tempRecStr = initialise_ReconStructure(fragment_set->lFragment, nCalpha, heapsize * fragment_set->nFragment);
    /** if verbose mode reconstructed structures RMSD are logged */
    if (verbose && globalfit)
        verboseRMSDvec = safe_malloc(heapsize * (nCalpha - fragment_set->lFragment + 1) * sizeof(float));

    /*________________________________________________________________________*/
    /* local fit */
    localfit_encode(recStr, refcoord, fragment_set, fragment_str, localfit_rmsd_array);
    sprintf(frameDesc, "%8.3f", fr->time);

    /*________________________________________________________________________*/
    /* add sequence to sequence set */
    if (entropy || xpmoutput)
		add_sequence_to_set(localfastaSequenceSet, frameDesc, recStr[0].encodedString, (nCalpha - fragment_set->lFragment + 1), *ptr_fastaIndex);

    /* record local fit average RMSD */
    fprintf(logFile, "Time: %10s\tRMSD(LF):%8.3f\tencoding: %s\n", frameDesc, recStr[0].rmsd, recStr[0].encodedString);
    fprintf(localFile, "%s\n", recStr[0].encodedString);
    if (fasta)
	    fprintf(fastalocalFile, ">Time %8.3f\n%s\n", fr->time, recStr[0].encodedString);

    /*________________________________________________________________________*/
    /** output local RMSD array */
    for (i = 0; i < (nCalpha - fragment_set->lFragment + 1); ++ i){
        fprintf(rmsdlocalFile, "%8.3f\n", localfit_rmsd_array[i]);
    }

    fflush(rmsdlocalFile);
    fflush(localFile);
    fflush(logFile);
    if (fasta)
	    fflush(fastalocalFile);

    if (globalfit){
    	/*____________________________________________________________________*/
	    /* brute force global fit */
    	for(i = 0; i < (nCalpha - fragment_set->lFragment + 1); ++i){
    		globalfit_encode(recStr, tempRecStr, refcoord, fragment_set, fragment_str, i, heapsize);
    		/* if verbose mode record reconstructed structures RMSD */
		    if (verbose){
		    	for(j = 0; j < heapsize; ++j)
		    		verboseRMSDvec[i * heapsize + j] = tempRecStr[0].rmsd;
		    }
    	}

	    /** update global RMSD array */
	    for (i = 0; i < nCalpha; ++ i){
	        globalfit_rmsd_array[i] = coord_rmsd(&(recStr[0].coord[i]), &(refcoord[i]));
		        fprintf(rmsdglobalFile, "%8.3f\n", globalfit_rmsd_array[i]);
	    }

	    /*________________________________________________________________________*/
	    /* get atom info from topology */
	    init_t_atoms(&useatoms, nCalpha, FALSE);
	    for(i = 0; i < nCalpha; i++) {
		useatoms.resinfo[i] = top->atoms.resinfo[i];
	    	useatoms.atomname[i] = top->atoms.atomname[index[i]];
	    	useatoms.atom[i] = top->atoms.atom[index[i]];
	    }
	    useatoms.nres = nCalpha;
	    useatoms.nr = nCalpha;

	    /*________________________________________________________________________*/
	    /* get coord from reconstructed structure */
	    snew(xmem, nCalpha);
	    frout.natoms = nCalpha;
	    frout.x = xmem;

	    for (i = 0; i < nCalpha; ++ i){
		        frout.x[i][XX] = recStr[0].coord[i].x / NM2ANG;
		        frout.x[i][YY] = recStr[0].coord[i].y / NM2ANG;
		        frout.x[i][ZZ] = recStr[0].coord[i].z / NM2ANG;
	    }

	    /*________________________________________________________________________*/
	    /* add sequence to sequence set */
		add_sequence_to_set(globalfastaSequenceSet, frameDesc, recStr[0].encodedString, (nCalpha - fragment_set->lFragment + 1), *ptr_fastaIndex);

	    /*________________________________________________________________________*/
	    /* save best fitting result */
		fprintf(pdbglobalFile, "REMARK    GENERATED BY G_SA_ENCODING\n");
		fprintf(pdbglobalFile, "REMARK    GLOBAL FIT RMSD: %8.3f\n", recStr[0].rmsd);
 		write_pdbfile(pdbglobalFile, "RECONSTRUCTED", &useatoms, frout.x, fr->ePBC, fr->box, ' ', fr->step, NULL, TRUE);

	    /* record best fitting result */
	    fprintf(logFile, "Time: %10s\tRMSD(GF):%8.3f\tencoding: %s\n", frameDesc, recStr[0].rmsd, recStr[0].encodedString);
	    fprintf(globalFile, "%s\n", recStr[0].encodedString);
	    if (fasta)
		    fprintf(fastaglobalFile, ">Time %8.3f\n%s\n", fr->time, recStr[0].encodedString);

	    fflush(rmsdglobalFile);
	    fflush(globalFile);
	    fflush(logFile);
	    if (fasta)
		    fflush(fastaglobalFile);
    }

    /*________________________________________________________________________*/
    /** update number of sequences */
    *ptr_fastaIndex = *(ptr_fastaIndex) + 1;

    /*________________________________________________________________________*/
    /** free the memory */

    /** free array of reconstructed structures */
    free_ReconStructure(recStr, heapsize);
    if (globalfit){
	    free_ReconStructure(tempRecStr, heapsize * fragment_set->nFragment);
	    free_t_atoms(&useatoms, FALSE);
    }

    /** free ref and reconstructed coords structure */
    free(refcoord);

    /** reconstructed structure RMSD vector */
    if (verbose && globalfit)
        free(verboseRMSDvec);

    /** RMDS arrays **/
    free(localfit_rmsd_array);
    if (globalfit)
	    free(globalfit_rmsd_array);

    /* memory for global fit reconstruction */
    if (globalfit){
        sfree(gnames);
        sfree(grps);
        sfree(index);
        sfree(xmem);
    }

}

/*____________________________________________________________________________*/
/** initialize xpm mapping */
void initialize_xpm_SA_mapping(t_mapping *SA_mapping, FragmentSet *fragment_set, SA_colour *colours){

	int i;
	char SAstring[2];

	for(i = 0; i < fragment_set->nFragment; ++i){
                SAstring[0] = fragment_set->codeOrder[i];
                SAstring[1] = '\0';
		SA_mapping[i].code.c1 = fragment_set->codeOrder[i];
		SA_mapping[i].code.c2 = 0;
		SA_mapping[i].desc = strdup(SAstring);
		if (colours != NULL){
			SA_mapping[i].rgb.r = colours->colour[i].r;
			SA_mapping[i].rgb.g = colours->colour[i].g;
			SA_mapping[i].rgb.b = colours->colour[i].b;
		}else{
			SA_mapping[i].rgb.r = ((float) (fragment_set->nFragment - (i + 1))) / fragment_set->nFragment;
			SA_mapping[i].rgb.g = 0.0;
			SA_mapping[i].rgb.b = 1.0;
		}
	}

}

/*____________________________________________________________________________*/
/** initialize xpm matrix */
void initialize_xpm_matrix(t_matrix *xpm_mat, int nx, int ny){

	int i;

	xpm_mat->flags = 0;
	xpm_mat->nx = nx;
	xpm_mat->ny = ny;
	xpm_mat->bDiscrete = TRUE;
	xpm_mat->axis_x = 0;
	xpm_mat->axis_y = 0;
	srenew(xpm_mat->axis_x, nx);
	for(i = 0; i < nx; ++i)
		xpm_mat->axis_x[i] = (float) i;
	srenew(xpm_mat->axis_y, ny);
	for(i = 0; i < ny; ++i)
		xpm_mat->axis_y[i] = (float) i;

}

/*____________________________________________________________________________*/
/** free xpm matrix */
void free_xpm_matrix(t_matrix *xpm_mat, int nmatcol){

    int i;

    sfree(xpm_mat->axis_x);
    sfree(xpm_mat->axis_y);

    sfree(xpm_mat->map);

	for(i = 0; i < nmatcol; ++i)
		sfree(xpm_mat->matrix[i]);
	sfree(xpm_mat->matrix);

}

/*____________________________________________________________________________*/
/** create and save xpm matrix */
void save_SA_xpm_matrix(FILE *xpmout, SeqSet *fastaSequenceSet, int sequenceLength, FragmentSet *fragment_set, SA_colour *SA_colour_set, int nColourSets){

	int i,j; /* index */
	t_matrix xpm_mat; /* matrix for xpm printout */
	t_xpmelmt c; /* placeholder for mapping instance */
	int colourIndex = -1;

	initialize_xpm_matrix(&xpm_mat, fastaSequenceSet->nSequences, sequenceLength);
	xpm_mat.map = safe_malloc(sizeof(t_mapping) * fragment_set->nFragment);

	for(i = 0; i < nColourSets; ++i)
		if (strcmp(SA_colour_set[i].alphabetName, fragment_set->setname) == 0)
			colourIndex = i;
	if (colourIndex > -1)
		initialize_xpm_SA_mapping(xpm_mat.map, fragment_set, &(SA_colour_set[colourIndex]));
	else
		initialize_xpm_SA_mapping(xpm_mat.map, fragment_set, NULL);

	xpm_mat.nmap = fragment_set->nFragment;

	sprintf(xpm_mat.title, "SA xpm Matrix");
	sprintf(xpm_mat.label_x, "frame");
	sprintf(xpm_mat.label_y, "Fragment");
	sprintf(xpm_mat.legend, "%s", fragment_set->setname);

	snew(xpm_mat.matrix, fastaSequenceSet->nSequences);
	for(i = 0; i < fastaSequenceSet->nSequences; ++i)
		snew(xpm_mat.matrix[i], sequenceLength);

	c.c2 = 0;
	for(i = 0; i < fastaSequenceSet->nSequences; ++i)
		for(j = 0; j < sequenceLength; ++j){
			c.c1 = fastaSequenceSet->sequence[i].res[j];
			xpm_mat.matrix[i][j] = max(0, searchcmap(xpm_mat.nmap, xpm_mat.map, c));
		}

	write_xpm_m(xpmout,xpm_mat);


    free_xpm_matrix(&xpm_mat, fastaSequenceSet->nSequences);
}

/*____________________________________________________________________________*/
/** calculate profile */
void calculate_profile(const char *outFilename, int sequenceLength, SeqSet *fastaSequenceSet, FragmentSet *fragment_set, gmx_bool xpmoutput){

	int i,j = 0; /* index */
	char *iString; /* string placeholders */
	Set iCodeSet; /* alphabet code set for position i */
	FILE *outprofile; /* file handle for entropy data */

    /* initialize codeSets */
    iCodeSet.nElements = 0;
    iCodeSet.element = safe_malloc(sizeof(Element));

    initialize_code_from_string(fragment_set->codeOrder, &iCodeSet);

	outprofile = safe_open(outFilename, "w");
	fprintf(outprofile, "         ");
	for(j = 0; j < fragment_set->nFragment; ++j){
		fprintf(outprofile, " %8c", fragment_set->codeOrder[j]);
	}
	fprintf(outprofile, "\n");

	for (i = 0; i < sequenceLength; ++ i){
		fprintf(outprofile, " %8d", (i + 1));
		iString = get_string_from_column(fastaSequenceSet, i);
		record_probabilities_from_string(iString, &iCodeSet);

		for(j = 0; j < iCodeSet.nElements; ++j){
			fprintf(outprofile, " %8.3f", iCodeSet.element[j].prob);
		}
		fprintf(outprofile, "\n");

	    free(iString);
		reset_probabilities(&iCodeSet);

	}

	free(iCodeSet.element);
	fclose(outprofile);

}

/*____________________________________________________________________________*/
/** calculate position entropy */
void entropy_per_position(const char *outFilename, int sequenceLength, SeqSet *fastaSequenceSet, const output_env_t oenv){

	int i = 0; /* index */
	char *iString; /* string placeholders */
	Set iCodeSet; /* alphabet code set for position i */
	float H = 0.0; /* Shannon's entropy */
	FILE *outEntropy; /* file handle for entropy data */

	outEntropy = xvgropen(outFilename, "Entropy", "Fragment", "H / bit", oenv);
	for (i = 0; i < sequenceLength; ++ i){
		iString = get_string_from_column(fastaSequenceSet, i);
		extract_code_from_string(iString, &iCodeSet);
		H = Shannon(&iCodeSet);
		fprintf(outEntropy, "%5d %8.3f\n", (i + 1), H);
	    free(iCodeSet.element);
	    free(iString);
	}
	fclose(outEntropy);
}


/*____________________________________________________________________________*/
/** main */
int main(int argc,char *argv[])
{
    /*________________________________________________________________________*/
    /* static variables */
    const char *desc[] = {
        "[TT]g_sa_encode[tt] reads a trajectory file ([TT]-f[tt]) and translate it into an alignment",
        "of structural strings using a code made by a set of representative 4-residue",
        "fragments (Structural Alphabet or SA). The M32K25 SA ([TT]-alphabet[tt]) is currently",
        "supported (Pandini et al., BMC Bioinformatics [BB]11[bb], 97 (2010)).[PAR]",

        "The letter assignment for each position and time is written to an ASCII file",
        "[TT].out[tt] ([TT]-strlf[tt] for local fit and [TT]-strgf[tt] for global fit mode). The assignment is",
        "done by matching the best-fitting alphabet fragment to stretches of 4",
        "consecutive C-alpha atoms in the protein (local fit, default). This allows for",
        "non-exact fragment overlaps. A seamless reconstruction of the overall protein",
        "structure with exact overlaps can be enforced by selecting the global fit",
        "option ([TT]-global[tt]). For global fit encoding, a heap size can be specified",
        "([TT]-heapsize[tt]), with larger values producing more accurate reconstructions (Park",
        "and Levitt, J Mol Biol [BB]249[bb], 493-507 (1995)). Please be aware that, depending",
        "on the heap size value, the global fit mode can be significantly more",
        "time-consuming than the local fit one.[PAR]",

        "In both modes, the RMSD between each fragment in the protein structure and the",
        "assigned letter is written for each frame to an [TT].xvg[tt] file ([TT]-rmsdlf[tt] for local",
        "fit and [TT]-rmsdgf[tt] for global fit mode).[PAR]",
        
        "The encoding is performed on the C-alpha atom group defined from the input",
        "reference structure ([TT]-s[tt]). A protein of [TT]n-[tt]residues is encoded in strings of",
        "length [TT]n-3[tt], with the [TT]i[tt]-th letter encoding for the C-alpha[TT][i][tt] to",
        "C-alpha[TT][i + 3][tt] fragment. Each line in the alignment corresponds to a different frame",
        "in the trajectory. Multiple chains or selections with non-consecutive C-alpha",
        "atoms are not currently supported. Each chain should be encoded separately.[PAR]",
        
        "Option [TT]-fasta[tt] writes the alignment in fasta format.[PAR]",
        
        "Option [TT]-xpm[tt] writes the alignment in [TT].xpm[tt] format: this can be converted to",
        "postscript with [TT]xpm2ps[tt] for an immediate visualization of the time evolution",
        "of the encoding for each position. The name of the [TT].xpm[tt] file can be specified",
        "with the [TT]-xpmlf[tt] (local fit) or [TT]-xpmgf[tt] (global fit) flag.[PAR]", 
        
        "Basic statistical analyses of the alignment can also be requested:[PAR]",

        "Option [TT]-entropy[tt] writes the Shannon entropy of each position to a [TT].xvg[tt] file", 
        "([TT]-Hlf[tt] for local fit or [TT]-Hgf[tt] for global fit).[PAR]",

        "Option [TT]-profile[tt] writes the relative frequency of each letter at each position",
        "to an ASCII [TT].dat[tt] file ([TT]-prolf[tt] for local fit or [TT]-progf[tt] for global fit)."
    };

    static char* userCodeString =
    		"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    static int heapsize = 1000;
    static gmx_bool fasta = FALSE;
    static gmx_bool globalfit = FALSE;
    static gmx_bool entropy = FALSE;
    static gmx_bool profile= FALSE;
    static gmx_bool xpmoutput= FALSE;
    static gmx_bool verbose = FALSE;
    int alphabet=0;
    static char *alphabetname[] = {NULL, "M32K25", NULL};
    enum {a_null, a_M32K25, a_nr};
    Str *fragment_str; /** structure array for alphabet fragments */
    FragmentSet *fragment_set; /** data structure for fragment set */
    FILE *rmsdlocalFile = 0; /* rmsd local filehandle */
    FILE *localFile = 0; /* local encoding  filehandle */
    FILE *fastalocalFile = 0; /* fasta local encoding  filehandle */
    FILE *rmsdglobalFile = 0; /* global encoding filehandle */
    FILE *globalFile = 0; /* fasta global encoding filehandle */
    FILE *fastaglobalFile = 0; /* rmsd global filehandle */
    FILE *pdbglobalFile = 0; /* pdb global filehandle */
    FILE *logFile = 0; /* log filehandle */
    FILE *xpmout = 0; /* xpm filehandle */
    char basenamelf[64]; /* basename */
    char basenamegf[64]; /* basename */
    char fastalocalFilename[64]; /* fasta local filename */
    char fastaglobalFilename[64]; /* fasta local filename */
	SeqSet localfastaSequenceSet; /* local fit fasta sequence set */
	SeqSet globalfastaSequenceSet; /* global fit fasta sequence set */
	int fastaIndex = 0; /* index of fasta sequences */
	int sequenceLength = 0; /* length of encoded sequence */
	SA_colour *SA_colour_set; /* placeholder for colour set  */
	int nColourSets; /* n colour set available */
	int chainBreak = 0; /* flag for chain breaks */

    /*________________________________________________________________________*/
    /* Extra arguments */
  
    t_pargs pa[] = {
    		{ "-alphabet", FALSE, etENUM, {&alphabetname},
    				"Structural Alphabet"
    		},
    		{ "-heapsize", FALSE, etINT, {&heapsize},
    				"Number of paths to grow in global fit mode"
    		},
    		{ "-fasta", FALSE, etBOOL, {&fasta},
    				"Generate also fasta output"
    		},
    		{ "-global", FALSE, etBOOL, {&globalfit},
    				"Perform also global fit encoding"
    		},
    		{ "-entropy", FALSE, etBOOL, {&entropy},
    				"Calculate entropy per position"
    		},
    		{ "-profile", FALSE, etBOOL, {&profile},
    				"Calculate frequency profile"
    		},
    		{ "-xpm", FALSE, etBOOL, {&xpmoutput},
    				"Output xpm matrix"
    		},
    		{ "-verbose", FALSE, etBOOL, {&verbose},
    				"Verbose output"
    		}
    };
  
    /*________________________________________________________________________*/
    /* gromacs variables */

    t_topology top;
    t_atoms    *atoms=NULL;
    int        ePBC;
    char       title[STRLEN];
    t_trxframe fr;
    rvec       *xtop;
    matrix     box;
    t_trxstatus *status;
    int        flags = TRX_READ_X;
    int        i,j;
    output_env_t oenv;

    /*________________________________________________________________________*/
    /* gromacs file variables */

    t_filenm fnm[] = {
    		{ efTPS,  NULL,  NULL, ffREAD },             /* topology */
    		{ efTRX, "-f", NULL, ffREAD },               /* trajectory */
    		{ efXVG, "-rmsdlf", "lf_rmsd", ffOPTWR},     /* local fit rmsd */
    		{ efOUT, "-strlf", "lf_str", ffWRITE},       /* local fit */
    		{ efXPM, "-xpmlf", "lf", ffOPTWR},           /* local xpm matrix */
    		{ efXVG, "-Hlf", "lf_entropy", ffOPTWR},     /* local fit entropy */
    		{ efDAT, "-prolf", "lf_prof", ffOPTWR},      /* local fit profile */
    		{ efXVG, "-rmsdgf", "gf_rmsd", ffOPTWR},     /* global fit rmsd */
    		{ efOUT, "-strgf", "gf_str", ffOPTWR},       /* global fit */
    		{ efXPM, "-xpmgf", "gf", ffOPTWR},           /* local xpm matrix */
    		{ efXVG, "-Hgf", "gf_entropy", ffOPTWR},     /* global fit entropy */
    		{ efDAT, "-progf", "gf_prof", ffOPTWR},      /* global fit profile */
    		{ efPDB, "-pdbgf", "gf", ffOPTWR},           /* global fit pdb */
    		{ efLOG, "-log", "log.log", ffOPTWR}         /* log */
    };

#define NFILE asize(fnm)

    /*________________________________________________________________________*/
    /* Print (C) */
    //CopyRight(stderr,argv[0]);

    /*________________________________________________________________________*/
    /* Parse args */
    parse_common_args(&argc, argv, PCA_CAN_TIME, NFILE, fnm,
    		asize(pa), pa, asize(desc), desc, 0, NULL,
                &oenv);

    /*________________________________________________________________________*/
    /* check alphabet */
    alphabet=1;
    while (alphabet < a_nr &&
    		strcasecmp(alphabetname[0], alphabetname[alphabet])!= 0)
    	alphabet++;
    if (alphabet == a_nr)
    	gmx_fatal(FARGS,"Invalid alphabet");

    fprintf(stderr,"Using %s alphabet for encoding\n", alphabetname[0]);

    /*________________________________________________________________________*/
    /* Read topology */
  
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
    sfree(xtop);

    /*________________________________________________________________________*/
    /** Load the fragment set data into fragment_set and fragment_str */
    fragment_set = safe_malloc(sizeof(FragmentSet));
    fragment_str = load_fragment_data(0, userCodeString, alphabetname[0],
    		fragment_set, 0);

    /* The first time we read data is a little special */
    read_first_frame(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
    fprintf(stderr, "\n");

    /* open log, output and rmsd output files */
    logFile = opt2FILE("-log", asize(fnm), fnm, "w");
    rmsdlocalFile = xvgropen(opt2fn("-rmsdlf", asize(fnm), fnm), "local fit RMSD", "Fragment", "RMSD (\\cE\\C)", oenv);
    localFile = opt2FILE("-strlf", asize(fnm), fnm, "w");
    if (fasta){
	    strncpy(basenamelf, opt2fn("-strlf", asize(fnm), fnm), strlen(opt2fn("-strlf", asize(fnm), fnm)) - 4);
    	basenamelf[strlen(opt2fn("-strlf", asize(fnm), fnm)) - 4] = '\0';
    	strcpy(fastalocalFilename, basenamelf);
	    fastalocalFile = safe_open(strcat(fastalocalFilename, ".fasta"), "w");
    }
    if (globalfit){
	    rmsdglobalFile = xvgropen(opt2fn("-rmsdgf", asize(fnm), fnm), "global fit RMSD", "Residue", "RMSD (\\cE\\C)", oenv);
	    globalFile = opt2FILE("-strgf", asize(fnm), fnm, "w");
	    if (fasta){
		    strncpy(basenamegf, opt2fn("-strgf", asize(fnm), fnm), strlen(opt2fn("-strgf", asize(fnm), fnm)) - 4);
	    	basenamegf[strlen(opt2fn("-strgf", asize(fnm), fnm)) - 4] = '\0';
	    	strcpy(fastaglobalFilename, basenamegf);
		    fastaglobalFile = safe_open(strcat(fastaglobalFilename, ".fasta"), "w");
	    }
	    pdbglobalFile = opt2FILE("-pdbgf", asize(fnm), fnm, "w");
    }

    /*________________________________________________________________________*/
	/* allocate memory for new sequence */
	localfastaSequenceSet.sequence = safe_malloc(sizeof(Seq));
	localfastaSequenceSet.nSequences = 0;

	if (globalfit){
		globalfastaSequenceSet.sequence = safe_malloc(sizeof(Seq));
		globalfastaSequenceSet.nSequences = 0;
	}

    /* This is the main loop over frames */
    do {
	    if (globalfit)
	        fprintf(rmsdglobalFile, "&\n");
	    encode_structure(fragment_set, fragment_str, &top, &fr, title,
	    		rmsdlocalFile, localFile, fastalocalFile,
	    		rmsdglobalFile, globalFile, fastaglobalFile, pdbglobalFile,
	    		logFile, heapsize, fasta, globalfit, verbose, entropy, xpmoutput,
	    		&fastaIndex, &localfastaSequenceSet, &globalfastaSequenceSet, &chainBreak);
        fprintf(rmsdlocalFile, "&\n");
    }while(read_next_frame(oenv,status,&fr));

    thanx(stderr);

    /* sequence length */
   	sequenceLength = localfastaSequenceSet.sequence[0].length;

    /*________________________________________________________________________*/
    /* calculate per position entropy */
    if (entropy){
		entropy_per_position(opt2fn("-Hlf", asize(fnm), fnm), sequenceLength, &localfastaSequenceSet, oenv);
    	if (globalfit){
			entropy_per_position(opt2fn("-Hgf", asize(fnm), fnm), sequenceLength, &globalfastaSequenceSet, oenv);
    	}
    }

    /*________________________________________________________________________*/
    /* calculate profile */
    if (profile){
		calculate_profile(opt2fn("-prolf", asize(fnm), fnm), sequenceLength, &localfastaSequenceSet, fragment_set, xpmoutput);
    	if (globalfit){
			calculate_profile(opt2fn("-progf", asize(fnm), fnm), sequenceLength, &globalfastaSequenceSet, fragment_set, xpmoutput);
    	}
    }

    /*________________________________________________________________________*/
    /* print xpm map of encoding */

	SA_colour_set = load_colour_set(&nColourSets);

	if (xpmoutput){
		xpmout = opt2FILE("-xpmlf", asize(fnm), fnm, "w");
		save_SA_xpm_matrix(xpmout, &localfastaSequenceSet, sequenceLength, fragment_set, SA_colour_set, nColourSets);
		fclose(xpmout);
		if (globalfit){
			xpmout = opt2FILE("-xpmgf", asize(fnm), fnm, "w");
			save_SA_xpm_matrix(xpmout, &globalfastaSequenceSet, sequenceLength, fragment_set, SA_colour_set, nColourSets);
			fclose(xpmout);
		}
	}


    /*________________________________________________________________________*/
	/** close file handles */
    fclose(rmsdlocalFile);
    fclose(localFile);
    fclose(logFile);
    if (globalfit){
	    fclose(rmsdglobalFile);
	    fclose(globalFile);
	    fclose(pdbglobalFile);
    }
    if (fasta){
    	fclose(fastalocalFile);
    	if (globalfit)
	    	fclose(fastaglobalFile);
    }

    /*________________________________________________________________________*/
    /* free memory */

    /* structure */
    for (i = 0; i < fragment_set->nFragment; ++ i)
        free(fragment_str[i].atom);
    free(fragment_str);

    /* fragment set */
    free_float_matrix3D(fragment_set->coord_values, fragment_set->nFragment, fragment_set->lFragment);
    free(fragment_set);

    /* colours */
	free_colour_set(SA_colour_set);

    /* fasta sequence set */
    for(i = 0; i < localfastaSequenceSet.nSequences; ++i){
    	free(localfastaSequenceSet.sequence[i].res);
    	free(localfastaSequenceSet.sequence[i].name);
    }
    free(localfastaSequenceSet.sequence);
    if (globalfit){
	    for(i = 0; i < globalfastaSequenceSet.nSequences; ++i){
	    	free(globalfastaSequenceSet.sequence[i].res);
	    	free(globalfastaSequenceSet.sequence[i].name);
	    }
	    free(globalfastaSequenceSet.sequence);
    }

    /* filenames */
    for(i = 0; i<NFILE; i++){
        for(j = 0; j<fnm[i].nfiles; j++){
            free(fnm[i].fns[j]);
        }   
    }   

    /* topology */
    done_top(&top);

    return 0;
}

