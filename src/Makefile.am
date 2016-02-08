AUTOMAKE_OPTIONS=subdir-objects
LIBS = -lgromacs -lm -lgsl -lgslcblas
AM_CFLAGS = -Wall -Wno-unused -funroll-all-loops
AM_CPPFLAGS = -I$(GMXDIR)/include -I$(GMXDIR)/include/gromacs
AM_LDFLAGS = -L$(GMXDIR)/lib -L/opt/local/lib $(LIBS)
#bin_PROGRAMS = g_sa_encode g_sa_analyze
bin_PROGRAMS = g_sa_analyze
#g_sa_encode_SOURCES = g_sa_encode.c\
#	general/matrix.c general/matrix.h\
#	general/safe.c general/safe.h\
#	general/mergesort.c general/mergesort.h\
#	sequence/encode.c sequence/encode.h\
#	sequence/sequence.c sequence/sequence.h\
#	statistics/code.c statistics/code.h\
#	statistics/probability.c statistics/probability.h\
#	structure/fragment_colour.c structure/fragment_colour.h\
#	structure/fragments.c structure/fragments.h\
#	structure/fragments_coords.h\
#	structure/getfragments.c structure/getfragments.h\
#	structure/kabsch.c structure/kabsch.h\
#	structure/pdb_structure.h\
#	structure/quatRMSD.c structure/quatRMSD.h\
#	structure/transform_segment.c structure/transform_segment.h\
#	structure/vec2mat.c structure/vec2mat.h\
#	legacyheaders/macros.h\
#	includeheaders/types/simple.h\
#	legacyheaders/typedefs.h
g_sa_analyze_SOURCES= g_sa_analyze.c\
	general/matrix.c general/matrix.h\
	general/safe.c general/safe.h\
	general/mergesort.c general/mergesort.h\
	sequence/sequence.c sequence/sequence.h\
	statistics/code.c statistics/code.h\
	statistics/probability.c statistics/probability.h\
	structure/fragment_colour.c structure/fragment_colour.h\
	structure/fragments.c structure/fragments.h\
	structure/fragments_coords.h\
	structure/getfragments.c structure/getfragments.h\
	structure/kabsch.c structure/kabsch.h\
	structure/pdb_structure.h\
	structure/quatRMSD.c structure/quatRMSD.h\
	structure/transform_segment.c structure/transform_segment.h\
	structure/vec2mat.c structure/vec2mat.h\
	partition/label.c partition/label.h\
	partition/object.c partition/object.h\
	partition/partition.c partition/partition.h\
	partition/score.c partition/score.h\
	partition/value.c partition/value.h\
	legacyheaders/sysstuff.h\
	legacyheaders/typedefs.h\
	includeheaders/molfile_plugin.h\
	includeheaders/vmdio.h\
	includeheaders/vmdplugin.h\
	legacyheaders/types/simple.h\
	includeheaders/types/enums.h\
	includeheaders/types/block.h\
	includeheaders/types/symtab.h\
	includeheaders/types/idef.h\
	includeheaders/types/atoms.h\
	includeheaders/types/topology.h\
	includeheaders/types/energy.h\
	includeheaders/types/inputrec.h\
	includeheaders/types/ishift.h\
	includeheaders/types/graph.h\
	includeheaders/types/nrnb.h\
	includeheaders/types/nblist.h\
	includeheaders/types/nsgrid.h\
	includeheaders/types/forcerec.h\
	includeheaders/types/fcdata.h\
	includeheaders/types/mdatom.h\
	includeheaders/types/pbc.h\
	includeheaders/types/ifunc.h\
	includeheaders/types/group.h\
	includeheaders/types/state.h\
	includeheaders/types/shellfc.h\
	includeheaders/types/constr.h\
	includeheaders/types/matrix.h\
	includeheaders/types/oenv.h\
	includeheaders/types/trx.h\
	includeheaders/vec.h\
	includeheaders/maths.h\
	legacyheaders/macros.h\
	includeheaders/mpelogging.h\
	includeheaders/physics.h\
	includeheaders/statutil.h
	includeheaders/types/filenm.h\
	includeheaders/readinp.h\
	includeheaders/wman.h\
	includeheaders/pdbio.h\
	includeheaders/oenv.h\
	includeheaders/gmxfio.h\
	includeheaders/warninp.h\
	includeheaders/symtab.h\
	includeheaders/atomprop.h\
	includeheaders/xdrf.h\
	includeheaders/futil.h\
	legacyheaders/types/commrec.h\
	legacyheaders/gmxmpi.h\
	legacyheaders/types/idef.h
