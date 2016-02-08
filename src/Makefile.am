LIBS = -lmd -lgmx -lm -lgsl -lgslcblas
AM_CFLAGS = -Wall -Wno-unused -funroll-all-loops
AM_CPPFLAGS = -I$(GMXDIR)/include -I$(GMXDIR)/include/gromacs
AM_LDFLAGS = -L$(GMXDIR)/lib -L/opt/local/lib $(LIBS)
bin_PROGRAMS = g_sa_encode g_sa_analyze
g_sa_encode_SOURCES = g_sa_encode.c\
	general/matrix.c general/matrix.h\
	general/safe.c general/safe.h\
	general/mergesort.c general/mergesort.h\
	sequence/encode.c sequence/encode.h\
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
	structure/vec2mat.c structure/vec2mat.h
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
	partition/value.c partition/value.h
