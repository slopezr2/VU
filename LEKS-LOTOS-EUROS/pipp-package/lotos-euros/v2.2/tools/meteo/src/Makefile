#
# Makefile
#

# --- settings --------------------------------------

# host depended settings,
# probably created by script calling make:
include Makefile_settings

# collected include flags:
INCS = $(NETCDF_INCS)

# target executables:
c3po_select_x            = c3po_select.x
c3po_spectral_gradient_x = c3po_spectral_gradient.x
c3po_spg2rgg_x           = c3po_spg2rgg.x
c3po_coarselevs_x        = c3po_coarselevs.x


# --- rules -------------------------------------------

%.o: %.F90
	$(FC) -c -o $@ $(FFLAGS) $(INCS) $<
	@echo ' '


# --- help --------------------------------------------

help:
	@echo ''
	@echo 'Usage:'
	@echo '  make $(c3po_select_x)'
	@echo '  make $(c3po_spectral_gradient_x)'
	@echo '  make $(c3po_coarselevs_x)'
	@echo '  make clean'
	@echo ''

	
# --- targets ------------------------------------------

binas.o : binas.F90

BINAS_OBJS = binas.o

# *

go_fu.o : go_fu.F90 go.inc
go_print.o : go_print.F90 go.inc go_fu.o
go_system.o : go_system.F90 go.inc go_print.o
go_string.o : go_string.F90 go.inc go_print.o
go_date.o : go_date.F90 go.inc go_string.o go_print.o
go.o : go.F90 go.inc go_date.o go_string.o go_system.o go_print.o go_fu.o

GO_OBJS = go.o go_date.o go_string.o go_fu.o go_print.o go_system.o

# *

num_tools.o : binas.o
num.o : num_tools.o

NUM_OBJS = num.o num_tools.o

# *

c3po_grid_tools.o: binas.o
c3po_grid_rgg.o: go.o c3po_grid_tools.o num.o
c3po_grid_singleton.o:
c3po_grid_spg.o: go.o binas.o c3po_grid_singleton.o
c3po_grid.o: c3po_grid_rgg.o c3po_grid_spg.o

C3PO_GRID_OBJS = c3po_grid.o c3po_grid_rgg.o c3po_grid_singleton.o c3po_grid_spg.o c3po_grid_tools.o

# *

c3po_levs_hyb.o :
c3po_levs.o : c3po_levs_hyb.o

C3PO_LEVS_OBJS = c3po_levs.o c3po_levs_hyb.o

# *

c3po_file_nc.o: go.o c3po_levs.o
c3po_file_rgg.o: go.o c3po_file_nc.o c3po_grid.o
c3po_file_spg.o: go.o c3po_file_nc.o c3po_grid.o
c3po_file_convert.o: go.o c3po_file_spg.o c3po_file_rgg.o c3po_grid.o
c3po_file.o: c3po_file_nc.o c3po_file_rgg.o c3po_file_spg.o c3po_file_convert.o

C3PO_FILE_OBJS = c3po_file.o c3po_file_convert.o c3po_file_nc.o c3po_file_rgg.o c3po_file_spg.o

# *

c3po_select_main.o: go.o c3po_grid.o c3po_file.o
c3po_select.o: go.o c3po_select_main.o

C3PO_SELECT_OBJS = c3po_select.o c3po_select_main.o $(C3PO_FILE_OBJS) $(C3PO_LEVS_OBJS) $(C3PO_GRID_OBJS) $(NUM_OBJS) $(GO_OBJS)

$(c3po_select_x): c3po_select.o
	$(LD) -o $@ $(LDFLAGS) $(C3PO_SELECT_OBJS) $(NETCDF_LIBS)

# *

c3po_spectral_gradient_main.o: go.o c3po_grid.o c3po_file.o
c3po_spectral_gradient.o: go.o c3po_spectral_gradient_main.o

C3PO_SPECTRAL_GRADIENT_OBJS = c3po_spectral_gradient.o c3po_spectral_gradient_main.o $(C3PO_FILE_OBJS) $(C3PO_GRID_OBJS) $(GO_OBJS)

$(c3po_spectral_gradient_x): c3po_spectral_gradient.o
	$(LD) -o $@ $(LDFLAGS) $(C3PO_SPECTRAL_GRADIENT_OBJS) $(NETCDF_LIBS)

# *

c3po_spg2rgg_main.o: go.o c3po_grid.o c3po_file.o
c3po_spg2rgg.o: go.o c3po_spg2rgg_main.o

C3PO_SPG2RGG_OBJS = c3po_spg2rgg.o c3po_spg2rgg_main.o $(C3PO_FILE_OBJS) $(C3PO_GRID_OBJS) $(NUM_OBJS) $(GO_OBJS) $(BINAS_OBJS)

$(c3po_spg2rgg_x): c3po_spg2rgg.o
	$(LD) -o $@ $(LDFLAGS) $(C3PO_SPG2RGG_OBJS) $(NETCDF_LIBS)

# *

c3po_coarselevs_main.o: go.o c3po_grid.o c3po_levs.o c3po_file.o
c3po_coarselevs.o: go.o c3po_coarselevs_main.o

C3PO_COARSELEVS_OBJS = c3po_coarselevs.o c3po_coarselevs_main.o $(C3PO_FILE_OBJS) $(C3PO_LEVS_OBJS) $(C3PO_GRID_OBJS) $(NUM_OBJS) $(GO_OBJS)

$(c3po_coarselevs_x): c3po_coarselevs.o
	$(LD) -o $@ $(LDFLAGS) $(C3PO_COARSELEVS_OBJS) $(NETCDF_LIBS)

# *

clean:
	/bin/rm -f *.o *.mod *.x


# --- end ----------------------------------------------
