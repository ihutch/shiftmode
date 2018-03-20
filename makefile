#########################################################################
# This makefile should build on linux provided it has access to the X11 library
# which typically requires the development package, a fortran compiler,
# and git. It won't work on MSWindows.
SHELL=/bin/bash
ACCISPARENT= $(HOME)/src/
ACCISHOME=${ACCISPARENT}accis/
ACCISX=$(ACCISHOME)libaccisX.a
MODULES=shiftmode.o
LIBRARIES = -L$(ACCISHOME) -laccisX -lX11
COMPILE-SWITCHES = -Wall -O2
#########################################################################
ifeq ("$(FORTRAN)","")
# Configure compiler. Mostly one long continued bash script.
# Preference order mpif90, ftn, gfortran, f77
 FORTRAN:=\
$(shell \
 if which mpif90 >/dev/null 2>&1; then echo -n "mpif90";else\
  if which ftn >/dev/null 2>&1 ; then echo -n "ftn";else\
    if which gfortran >/dev/null 2>&1; then echo -n "gfortran";else\
     if which f77 >/dev/null 2>&1 ; then echo -n "f77";else\
	echo "Unable to decide compiler. Specify via FORTRAN=..." >&2; exit 1;\
     fi;\
    fi;\
  fi;\
 fi;\
)
endif
export FORTRAN
#########################################################################
 ACCISCHECK:=\
$(shell echo >&2 "Checking accis library...";\
 if [ -f "${ACCISX}" ] ; then echo "Library ${ACCISX} exists."; else\
   if [ -d "${ACCISPARENT}" ] ; then echo -n "src directory exists. ";\
     else mkdir ${ACCISPARENT} ; fi;\
   if [ -d "${ACCISHOME}" ] ; then echo -n "accis directory exists. ";\
     else cd ${ACCISPARENT};\
	git clone git@github.com:ihutch/accis.git; cd - >/dev/null; fi;\
   cd ${ACCISHOME}; make; cd - >/dev/null;\
   if [ -f "${ACCISX}" ] ; then echo -n "Made ${ACCISX}";\
     else echo "Error making ${ACCISX}"; fi;\
 fi;\
)
#########################################################################
default : $(ACCISX)
	@echo "$(ACCISCHECK)"

$(ACCISX) : $(ACCISHOME)Makefile
	cd $(ACCISHOME); make; cd -

%.o : %.f makefile ;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f

%.o : %.f90 makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f90

%.o : %.F makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.F

% : %.f $(ACCISX) ;
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f  $(LIBRARIES)

% : %.f90  makefile $(ACCISX) $(MODULES);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90 $(MODULES) $(LIBRARIES)

% : %.F$ (ACCISX) ;
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.F  $(LIBRARIES)


