#
# Makefile template
#
# Special targets:
#    make install bindir=../../run exe=lotos-euros.x
#


## ------------------------------------------------
## compiler and flags (template)
## ------------------------------------------------
#
## compiler and linker:
#FC = f90
#F77 = f77
#LINKER = $(FC)
#
## compile flags:
#FFLAGS = 
#
## linker flags:
#LDFLAGS =  
#
## library flags:
#LIBS = 


# ------------------------------------------------
# compiler and flags (included)
# ------------------------------------------------

# flags are written to include file by pycasso scripts:
include Makefile_flags


# ------------------------------------------------
# implicit rules
# ------------------------------------------------

# how to form object files from F90 source:
%.o: %.F90
	$(FC) -c -o $@ $(FFLAGS) $<
	@echo ' '

# how to form object files from f90 source:
%.o: %.f90
	$(FC) -c -o $@ $(FFLAGS) $<
	@echo ' '

# how to form object files from F source:
%.o: %.F
	$(F77) -c -o $@ $(FFLAGS) $<
	@echo ' '

# how to form object files from f source:
%.o: %.f
	$(F77) -c -o $@ $(FFLAGS) $<
	@echo ' '


## ------------------------------------------------
## explicit rules (template)
## ------------------------------------------------
#
#le_sia_isorropia2_isocom.o:
#	$(FC) -c -o $@ $(FFLAGS) -warn nodeclarations $<
#	@echo " "
#
#le_sia_isorropia2_isofwd.o:
#	$(FC) -c -o $@ $(FFLAGS) -warn nodeclarations $<
#	@echo " "
#
#le_sia_isorropia2_isorev.o:
#	$(FC) -c -o $@ $(FFLAGS) -warn nodeclarations $<
#	@echo " "


# ------------------------------------------------
# explicit rules (included)
# ------------------------------------------------

# explicit rules are written to include file by pycasso scripts:
include Makefile_rules


# ------------------------------------------------
# dependencies
# ------------------------------------------------

# dependencies are written to include file by pycasso scripts:
include Makefile_deps


# ------------------------------------------------
# targets
# ------------------------------------------------

install:
	@if [ ! -d $(bindir) ]; then /bin/mkdir -p $(bindir) ; fi
	/bin/cp $(exe) $(bindir)

clean:
	/bin/rm -f *.o *.mod
	/bin/rm -f *.x

# ------------------------------------------------
# end
# ------------------------------------------------

