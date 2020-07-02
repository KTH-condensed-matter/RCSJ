################# -*-makefile-*-
# Program parameters
#################

PROJECT  = JJarray
TARGETS	 = rcsj analyze rcsj_voltages find-phase-slips plot-v
PROGRAMS = $(TARGETS) show xiv xxiv slurm

DEST	= run

HDRS	= mt19937.h util.h graph.h Estimate.h TriSolve.h Table.h Average.h
SRCS	= rcsj.cc graph.cc icon.cc analyze.cc find-phase-slips.cc plot-v.cc
OTHER	= Makefile xrcsj show xrcsj-v xrcsj.tmpl go xrcsj-v-photons xiv

#################
# Targets
#################

all:	$(TARGETS)

rcsj:	rcsj.cc graph.o
	$(LINK.cc) -DPBC rcsj.cc graph.o -o $@ $(LDLIBS)

rcsj_voltages:	rcsj.cc graph.o
	$(LINK.cc) -DPRINT_VOLTAGE rcsj.cc graph.o -o $@ $(LDLIBS)

# analyze needs to be compiled without optimization flags...
analyze: COMPFLAGS=

.PHONY:	all

#################
# Compiler parameters
#################

INSTALL=install

### Gnu GCC compiler collection

IPATH	= -I/usr/X11/include
LDLIBS	= -L/usr/X11R6/lib -lX11
COMPFLAGS = -O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -funit-at-a-time -pipe -mtune=native
CCFLAGS	= $(COMPFLAGS) $(IPATH)
FFLAGS	= $(COMPFLAGS)
CFLAGS	= $(COMPFLAGS)
FC	= f77
LD	= c++
CCC	= c++
CC	= gcc

### Compiler flags for MacOS:

# IPATH   = -I/opt/X11/include 
# LDLIBS  = -L/opt/X11/lib -lX11
# COMPFLAGS = -Ofast -ffast-math -funroll-loops -funit-at-a-time -pipe -march=native -m64
# CCFLAGS = $(COMPFLAGS) $(IPATH)
# FFLAGS  = $(COMPFLAGS)
# CFLAGS  = $(COMPFLAGS)
# FC      = f77
# LD      = c++
# CCC     = c++
# CC      = cc

#################
# Target Rules
#################

clean:;		@rm -f *.o core

clobber:;	@rm -f *.o core tags TAGS $(TARGETS)

depend:;	makedepend -- $(CCFLAGS) $(CPPFLAGS) -- $(SRCS)

echo:;		@echo $(HDRS) $(SRCS) $(OTHER)

tags TAGS:	$(HDRS) $(SRCS); @etags $(HDRS) $(SRCS)

tar $(PROJECT).tar: $(HDRS) $(SRCS) $(OTHER)
		tar cvhlf $(PROJECT).tar $(MAKEFILE) $(HDRS) $(SRCS) $(OTHER)

install::	$(PROGRAMS:%=$(DEST)/%)

$(DEST)/%::	%
		@echo Installing $< in $(DEST).
		@mkdir -p $(DEST);
		@cp -fp $< $(DEST)

### Make rules

.SUFFIXES: .cc .o

COMPILE.cc=$(CCC) $(CCFLAGS) $(CPPFLAGS)  -c
LINK.cc=$(CCC) $(CCFLAGS) $(CPPFLAGS) $(LDFLAGS) 

.cc:
	$(LINK.cc) -o $@ $< $(LDLIBS)

.cc.o:
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

.cc.a:
	$(COMPILE.cc) -o $% $<
	$(AR) $(ARFLAGS) $@ $%
	$(RM) $%

###
