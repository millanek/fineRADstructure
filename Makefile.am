ACLOCAL_AMFLAGS = -I m4
AM_LDFLAGS = -lz

bin_PROGRAMS = finestructure RADpainter

#OPTIMIZATION = -O1 -g -Wall
OPTIMIZATION = -O3 -Wall -mfpmath=sse -msse -msse2 -funroll-loops -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations -lgsl -lgslcblas
finestructure_CXXFLAGS = $(GSL_CFLAGS) $(OPTIMIZATION) -Wall 
finestructure_LDFLAGS = $(OPTIMIZATION) 
finestructure_LDADD = $(GSL_LIBS)
finestructure_SOURCES = safegetline.h data.h fsxml.h state.h node.h inf1.h infadmixture.h infextract.h infextract2.h infextract3.h infextract4.h infextract5.h infmcmc.h fines.h finesfunctions.h rng.h safegetline.cpp data.cpp fsxml.cpp state.cpp node.cpp inf1.cpp infadmixture.cpp infextract.cpp infextract2.cpp infextract3.cpp infextract4.cpp infextract5.cpp infmcmc.cpp finesfunctions.cpp fines.cpp rng.cpp

RADpainter_CXXFLAGS = $(GSL_CFLAGS) -std=c++0x
RADpainter_LDADD = $(GSL_LIBS)
RADpainter_SOURCES = tagPainter.h paintSql.h utils.h tagPainter.cpp paintSql.cpp utils.cpp hapsFromVCF.h hapsFromVCF.cpp gzstream.h gzstream.cpp

check_SCRIPTS = profilerun.sh
TESTS = $(check_SCRIPTS)


profilerun.sh:
	echo './finestructure -s 532772982 -x 100 -y 50 testout' > profilerun.sh
	chmod +x profilerun.sh
