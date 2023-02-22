
# ------------------------------------------------
# help
# ------------------------------------------------

"""
Generate chemistry code for LOTOS-EUROS.
"""

# ------------------------------------------------
# common
# ------------------------------------------------

# get common stuff:
#   - logger instance
from genes_common import *


# ------------------------------------------------
# function
# ------------------------------------------------


def WriteSource( rcf, basekey, specs, reacs, labeling, labelled_specs, with_vbs, srcdir ) :

    """
    Write source files.
    """
    
    # modules:
    import subprocess
    
    # list of files to be written:
    srclist = rcf.get( basekey+'.files' )
    srcfiles = srclist.split()
    
    # loop over all files to be written:
    for srcfile in srcfiles :
        # info ...
        logger.info( '        write %s ...' % srcfile )
        # switch:
        if srcfile == 'le_indices.inc' :
            WriteSource_Indices( srcfile, rcf, specs )
        elif srcfile == 'le_chem_work.F90' :
            WriteSource_Work( srcfile, rcf, specs, reacs, labeling, labelled_specs, with_vbs )
        else :
            logger.error( 'do not know how to generate source file "%s"' % srcfile )
            raise Exception
        #endif
    #endfor
    
    # show generated files ?
    command = rcf.get( basekey+'.show.command' )
    # show ?
    if len(command) > 0 :
        # run:
        subprocess.check_call( command, shell=True )
    #endif
    
    ## testing ..
    #logger.warning( '' )
    #logger.warning( '*********************************************' )
    #logger.warning( 'break after writing code' )
    #logger.warning( '*********************************************' )
    #logger.warning( '' )
    #raise Exception
    
#enddef   # WriteSource


# ***


def WriteSource_Indices( indices_inc, rcf, specs ) :
    
    # modules:
    import os
    
    # tools:
    import utopya
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # settings
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # advect radicals ?
    advect_radicals = rcf.get( 'advection.radicals', 'bool' )
    
    # properties for which indices groups should be created:
    indgroup_props = rcf.get( 'genes.indgroup.props' ).split()
    
    # properties for which flag arrays should be created:
    flaglist_props = rcf.get( 'genes.flaglist.props' ).split()
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # indices include files
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    ## count radicals:
    #nrad = 0
    #for ispec in range(specs.nspec) :
    #    if 'radical' in specs.props[ispec] : nrad = nrad + 1
    ##endfor
    
    # count number of enabled tracers:
    nspec_all = 0
    for itr in range(specs.nspec) :
        if specs.enabled[itr] :
            nspec_all = nspec_all + 1
        #endif
    #endfor

    # count number of enabled tracers, excluding the special:
    nspec = 0
    for itr in range(specs.nspec) :
        # skip ?
        if not specs.enabled[itr] : continue
         # accumulated tracers are computed on output:
        if 'accum'    in specs.props[itr] : continue
        # biascorrected tracers used to be skipped,
        # but now included as 'normal' tracers:
        #if 'biascorr' in specs.props[itr] : continue
        # increase counter:
        nspec = nspec + 1
    #endfor
  
    # init empty:
    lines = []
    
    lines.append( '!#######################################################################' )
    lines.append( '!' )
    lines.append( '! LE indices' )
    lines.append( '!' )
    lines.append( '! Do not edit, automatically generated ...' )
    lines.append( '!' )
    lines.append( '!###############################################################################' )
    lines.append( '' )
    # counters:
    lines.append( '  ! number of tracers, including special:' )
    lines.append( '  integer, parameter  ::  nspec_all = %i' % nspec_all )
    lines.append( '' )
    lines.append( '  ! number of tracers, excluding special:' )
    lines.append( '  integer, parameter  ::  nspec = %i' % nspec )
    lines.append( '' )
    
    # indices:
    lines.append( '  ! tracer indices;' )
    lines.append( '  ! value -1 for disabled tracers:' )
    ispec = 0
    idum = 0
    for itr in range(specs.nspec) :
        # enabled ?
        if specs.enabled[itr] :
            # increase counter:
            ispec = ispec + 1
            # index for this tracer:
            indx = ispec
        else :
            # decrease index:
            idum = idum - 1
            # index for this tracer:
            indx = idum
        #endif
        # add line:
        lines.append( '  integer, parameter  ::  ispec_%-8s = %3i  ! %s' % (specs.name[itr],indx,specs.desc[itr]) )
    #endfor
    lines.append( '' )
    
    # indices:
    lines.append( '  ! tracer indices as variables, for backwards compatibility:' )
    ispec = 0
    idum = 0
    for itr in range(specs.nspec) :
        # enabled ?
        if specs.enabled[itr] :
            # increase counter:
            ispec = ispec + 1
            # index for this tracer:
            indx = ispec
        else :
            # decrease index:
            idum = idum - 1
            # index for this tracer:
            indx = idum
        #endif
        # add line:
        lines.append( '  integer             ::  i_%-8s = %3i  ! %s' % (specs.name[itr],indx,specs.desc[itr]) )
    #endfor
    lines.append( '' )
    
    # names:
    mlen = max(map(len,specs.name))
    lines.append( '  ! tracer names:' )
    lines.append( '  character(len=*), parameter  ::  specname(nspec_all) = (/ &' )
    # loop over species:
    sep = ', &'
    ispec = 0
    for itr in range(specs.nspec) :
        # enabled ?
        if not specs.enabled[itr] : continue
        # increase counter:
        ispec = ispec + 1
        # last one ?
        if ispec == nspec_all : sep = ' /)'
        # add line:
        fmt = "          '%-"+str(mlen)+"s'"+sep+"    ! %s"
        lines.append( fmt % (specs.name[itr].lower(),specs.name[itr]) )
    #endfor
    lines.append( '' )
    
    # units:
    lines.append( '  ! tracer units:' )
    lines.append( '  character(len=*), parameter  ::  specunit(nspec_all) = (/ &' )
    # maximum length:
    mlen = max( map(len,specs.units) )
    # loop over species:
    sep = ', &'
    ispec = 0
    for itr in range(specs.nspec) :
        # enabled ?
        if not specs.enabled[itr] : continue
        # increase counter:
        ispec = ispec + 1
        # last one ?
        if ispec == nspec_all : sep = ' /)'
        # set units:
        units = specs.units[itr]
        # add line:
        fmt = "          '%-"+str(mlen)+"s'"+sep+"    ! %s"
        lines.append( fmt % (units,specs.name[itr]) )
    #endfor
    lines.append( '' )
    
    # formula:
    lines.append( '  ! tracer mole mass:' )
    lines.append( '  real, parameter  ::  specmolm(nspec_all) = (/ &' )
    # loop over species:
    sep = ', &'
    ispec = 0
    for itr in range(specs.nspec) :
        # enabled ?
        if not specs.enabled[itr] : continue
        # increase counter:
        ispec = ispec + 1
        # last one ?
        if ispec == nspec_all : sep = ' /)'
        # compose molemass from formula:
        molm = ''
        atoms = specs.formula[itr].keys()
        if 'R' in atoms :
            molm = molm+' xm_dummy'
        else :
            # init counter for parts to be added, used to add '+' signs:
            ipart = 0
            # loop over atoms:
            for atom in atoms :
                # number of times this atom occures in formula:
                n = specs.formula[itr][atom]
                # trap zeros, since this is here it is probably needed:
                if n == 0 : continue
                # increase counter, add '+' sign if necessary:
                ipart = ipart + 1
                if ipart > 1 : molm = molm+' +'
                # add contribution of this atom to mole mass,
                # first trap special, otherwise a standard atom is assumed:
                if atom.startswith('mm') :
                    # translate 'mm150' to a mole-mase of 150e-3 kg/mole:
                    molm = molm+' %s*1.0e-3' % atom.lstrip('mm')
                else :
                    # use predefined parameter with known atom masses:
                    molm = molm+' xm_%s' % atom
                #endif
                # multiply with number of this atoms in molecule:
                if n != 1 : molm = molm+' * '+str(n)
            #endfor
        #endif
        # add line:
        lines.append( '      %s%s     ! %s' % (molm,sep,specs.name[itr]) )
    #endfor
    lines.append( '' )
    # aerosol densities:
    lines.append( '  ! aerosol densities:' )
    lines.append( '  real, parameter  ::  aerdens(nspec_all) = (/ &' )
    # maximum length:
    mlen = max( map(len,specs.units) )
    # loop over species:
    sep = ', &'
    ispec = 0
    for itr in range(specs.nspec) :
        # enabled ?
        if not specs.enabled[itr] : continue
        # increase counter:
        ispec = ispec + 1
        # last one ?
        if ispec == nspec_all : sep = ' /)'
        # set density:
        aerdens = specs.aerdens[itr]
        # add line:
        fmt = "          %-"+str(mlen)+"s"+sep+"    ! %s"
        lines.append( fmt % (aerdens,specs.name[itr]) )
    #endfor
    lines.append( '' )
    
    # aerosol hygroities:
    lines.append( '  ! aerosol hygroscopicities:' )
    lines.append( '  real, parameter  ::  aerhygro(nspec_all) = (/ &' )
    # maximum length:
    mlen = max( map(len,specs.units) )
    # loop over species:
    sep = ', &'
    ispec = 0
    for itr in range(specs.nspec) :
        # enabled ?
        if not specs.enabled[itr] : continue
        # increase counter:
        ispec = ispec + 1
        # last one ?
        if ispec == nspec_all : sep = ' /)'
        # set hygroity:
        aerhygro = specs.aerhygro[itr]
        # add line:
        fmt = "          %-"+str(mlen)+"s"+sep+"    ! %s"
        lines.append( fmt % (aerhygro,specs.name[itr]) )
    #endfor 
      
    lines.append( '' ) 
    # aerosol mode:
    lines.append( '  ! aerosol modes:' )
    lines.append( '  integer, parameter              ::  NO_AEROSOL_MODE     = 0' )
    lines.append( '  integer, parameter              ::  AEROSOL_NUCL_MODE   = 1' )
    lines.append( '  integer, parameter              ::  AEROSOL_AITKEN_MODE = 2' )
    lines.append( '  integer, parameter              ::  AEROSOL_ACCUM_MODE  = 3' )
    lines.append('   integer, parameter              ::  AEROSOL_FF_MODES    = 4' )
    lines.append( '  integer, parameter              ::  AEROSOL_FINE_MODES  = 5' )
    lines.append( '  integer, parameter              ::  AEROSOL_CCC_MODE    = 6' )
    lines.append( '  integer, parameter              ::  AEROSOL_CC_MODE     = 7' )
    lines.append( '  integer, parameter              ::  AEROSOL_COARSE_MODE = 8' )
    lines.append( '  integer, parameter              ::  AEROSOL_ALL_MODES   = 9' )
    lines.append( '  integer, parameter              ::  AEROSOL_POLLEN_MODE = 10' )
    lines.append( '  integer, parameter              ::  AEROSOL_ULTRA_FINE_MODE      = 11' )
    lines.append( '  integer, parameter              ::  AEROSOL_ULTRA_FINE_FINE_MODE = 12' )
    lines.append( '  integer, parameter              ::  AEROSOL_FINE_MEDIUM_MODE     = 13' )
    lines.append( '  integer, parameter              ::  AEROSOL_MEDIUM_COARSE_MODE   = 14' )    
    lines.append( '  integer, parameter              ::  N_AEROSOL_MODES     = 14' )
    lines.append( '  ! per tracer:' )
    lines.append( '  integer, parameter  ::  specmode(nspec_all) = (/ &' )
    # loop over species:
    sep = ', &'
    ispec = 0
    for itr in range(specs.nspec) :
        # enabled ?
        if not specs.enabled[itr] : continue
        # increase counter:
        ispec = ispec + 1
        # last one ?
        if ispec == nspec_all : sep = ' /)'
        # set mode string:
        if 'aerosol' in specs.props[itr] :
            if 'nucl_mode' in specs.props[itr] :
                mode = 'AEROSOL_NUCL_MODE  '
            elif 'aitken_mode' in specs.props[itr] :
                mode = 'AEROSOL_AITKEN_MODE'
            elif 'accum_mode' in specs.props[itr] :
                mode = 'AEROSOL_ACCUM_MODE '
            elif 'ff_modes' in specs.props[itr] :
                mode = 'AEROSOL_FF_MODES '
            elif 'fine_modes' in specs.props[itr] :
                mode = 'AEROSOL_FINE_MODES '
            elif 'ccc_mode' in specs.props[itr] :
                mode = 'AEROSOL_ccc_MODE'
            elif 'cc_mode' in specs.props[itr] :
                mode = 'AEROSOL_CC_MODE'                           
            elif 'coarse_mode' in specs.props[itr] :
                mode = 'AEROSOL_COARSE_MODE'
            elif 'all_modes' in specs.props[itr] :
                mode = 'AEROSOL_All_MODES  '
            elif 'pollen_mode' in specs.props[itr] :
                mode = 'AEROSOL_POLLEN_MODE'    
                
            elif 'ultra_fine_mode' in specs.props[itr] :    
                mode = 'AEROSOL_ULTRA_FINE_MODE'
            elif 'ultra_fine_fine_mode' in specs.props[itr] :
                mode = 'AEROSOL_ULTRA_FINE_FINE_MODE'
            elif 'fine_medium_mode' in specs.props[itr] :
                mode = 'AEROSOL_FINE_MEDIUM_MODE'
            elif 'medium_coarse_mode' in specs.props[itr] :
                mode = 'AEROSOL_MEDIUM_COARSE_MODE'
            else :
                logger.error( 'no aerosol mode in "%s" properties : %s' % (specs.name[itr],specs.props[itr]) )
                raise Exception
            #endif
        else :
            # no aerosol ...
            mode = 'NO_AEROSOL_MODE '
        #endif
        # add line:
        lines.append( '      %s%s    ! %s' % (mode,sep,specs.name[itr]) )
    #endfor
    lines.append( '' )
    
    # property selections:
    for prop in indgroup_props :
        # count enabled and not-special tracers with this property:
        names = []
        for itr in range(specs.nspec) :
            # not enabled ? then skip:
            if not specs.enabled[itr] : continue
            # special ? then skip:
            if 'accum'    in specs.props[itr] : continue
            #if 'biascorr' in specs.props[itr] : continue
            # for advection, check on radicals:
            if prop == 'advected' :
                if ('radical' in specs.props[itr]) and (not advect_radicals) : continue
                if 'biascorr' in specs.props[itr] : continue
            else :
                if prop not in specs.props[itr] : continue
            #endif
            # add:
            names.append( specs.name[itr] )
        #endfor
        # count:
        nname = len(names)
        # add counter:
        lines.append( '  ! number of %s tracers:' % prop )
        lines.append( '  integer, parameter  ::  n_%s = %i' % (prop,nname) )
        # collected indices:
        if nname == 0 :
            lines.append( '  ! dummy list of tracer indices:' )
            lines.append( '  integer, parameter  ::  ispecs_%s(1) = (/ -1 /)' % prop )
        else :
            lines.append( '  ! list of tracer indices:' )
            lines.append( '  integer, parameter  ::  ispecs_%s(n_%s) = (/ &' % (prop,prop) )
            sep = ', &'
            for iname in range(nname) :
                if iname == nname-1 : sep = ' /)'
                lines.append( '      ispec_%-8s%s   ! %s' % (names[iname],sep,names[iname]) )
            #endfor
        #endfor
        lines.append( '' )
    #endfor
    
    # flag lists:
    for prop in flaglist_props :
        # start array:
        lines.append( '  ! flags to check if tracer has "%s" property:' % prop )
        lines.append( '  logical, parameter  ::  tracer_is_%s(nspec_all) = (/ &' % prop )
        # spec index in model:
        ispec = 0
        # loop over all tracers:
        for itr in range(specs.nspec) :
            # skip if not enabled:
            if not specs.enabled[itr] : continue
            # increase index:
            ispec = ispec + 1
            # tracer has this property? set fortran label:
            if prop in specs.props[itr] :
                value = '.true.'
            else :
                value = '.false.'
            #endif
            # end of line mark:
            if ispec < nspec_all :
                eol = ', &'
            else :
                eol = ' /)'
            #endif
            # write:
            lines.append( '      %-8s %s   ! %s' % (value,eol,specs.name[itr]) )
        #endfor
        lines.append( '' )

        # also mapping from tracer to subset index:
        lines.append( '  ! index of ispec in "ispecs_%s" list:' % prop )
        lines.append( '  integer, parameter  ::  ispec2%s(nspec) = (/ &' % prop )
        # spec index in model:
        ispec = 0
        # mapping index:
        iprop = 0
        # loop over all tracers:
        for itr in range(specs.nspec) :
            # skip if not enabled or special:
            if not specs.enabled[itr] : continue
            if 'accum'    in specs.props[itr] : continue
            #if 'biascorr' in specs.props[itr] : continue
            # tracer has this property? set fortran label:
            if prop in specs.props[itr] :
                iprop = iprop + 1
                value = iprop
            else :
                value = -999
            #endif
            # increase index:
            ispec = ispec + 1
            # end of line mark:
            if ispec < nspec :
                eol = ', &'
            else :
                eol = ' /)'
            #endif
            # write:
            lines.append( '      %6i %s   ! %s' % (value,eol,specs.name[itr]) )
        #endfor
        lines.append( '' )
    #endfor
    
    ## update counter:
    #lines.append( '      ! update:' )
    #lines.append( '      i_in = i_in + %i' % ngas )
    #lines.append( '' )

    # add newlines:
    xlines = []
    for line in lines : xlines.append( line+'\n' )
    
    # update file if necessary ...
    utopya.UpdateTextFile( indices_inc, xlines )
    
    ## testing ..
    #raise Exception

#enddef  #  WriteSource_Indices
   
    
# ***


def WriteSource_Work( chemistry_module, rcf, specs, reacs, labeling, labelled_specs, with_vbs ) :

    """
    Write module with routines for reaction rates and iteration step.
    """
    
    # modules:
    import os
    
    # tools:
    import utopya
    
    # first write to temporary file, replace original only if different
    # to avoid re-make every time the file is generated:
    tempfile = '%s.tmp' % chemistry_module

    # fill:
    lines = []
    lines.append( '!#######################################################################' )
    lines.append( '!' )
    lines.append( '! Gas-phase chemistry routines.' )
    lines.append( '!' )
    lines.append( '! Do not edit, automatically generated ...' )
    lines.append( '!' )
    lines.append( '!###############################################################################' )
    lines.append( '!' )
    lines.append( '#define TRACEBACK write (gol,\'("in ",a," (",a,", line",i5,")")\') rname, __FILE__, __LINE__; call goErr' )
    lines.append( '#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if' )
    lines.append( '#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if' )
    lines.append( '!' )
    lines.append( '#include "le.inc"' )
    lines.append( '!' )
    lines.append( '!###############################################################################' )
    lines.append( '' )
    lines.append( 'module LE_Chem_Work' )
    lines.append( '' )
    lines.append( '  use GO, only : gol, goPr, goErr' )
    lines.append( '' )
    lines.append( '  implicit none' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  ! --- in/out -----------------------------------' )
    lines.append( '' )
    lines.append( '  private' )
    lines.append( '' )
    lines.append( '  public  ::  nreac' )
#    lines.append( '!  public  ::  ireac0_het' )
    lines.append( '' )
    # make indices of aux array public:
    lines.append( '  public  ::  naux' )
    for nam in reacs.aux_name :
        lines.append( '  public  ::  iaux_%s' % nam )
    #endfor
    lines.append( '' )
    lines.append( '  public  ::  cldadj' )
    lines.append( '  public  ::  mw_vbs' )
    lines.append( '' )
    # public procudures:
    lines.append( '  public  ::  LE_Chem_Work_Init, LE_Chem_Work_Done' )
    lines.append( '  public  ::  LE_Chem_Work_Rates' )
    lines.append( '  public  ::  LE_Chem_Work_Iter' )
    if ( labeling ) :
        lines.append( '  public  ::  n_nonzeros' )
    #endif        
    lines.append( '' )
    lines.append( '  ! reaction indices needed for M7:' )
    lines.append( '  public  ::  ireac_M7_R71, ireac_M7_R72' )
    lines.append( '' )
    lines.append( '  ! reaction indices needed for VBS:' )
    lines.append( '  public  ::  ireac_VBS_R63, ireac_VBS_R65, ireac_VBS_R66' )
    lines.append( '' )
    lines.append( '  ! --- const ------------------------------------' )
    lines.append( '' )
    lines.append( '  ! module name:' )
    lines.append( '  character(len=*), parameter  ::  mname = \'LE_Chem_Work\'' )
    lines.append( '' )
    lines.append( '  ! number of gas-phase reactions:' )
    lines.append( '  integer, parameter  ::  nreac = %i' % reacs.nreac )
    lines.append( '' )
    lines.append( '  ! auxilary reaction rate parameters :' )
    lines.append( '  integer, parameter  ::  naux = %i' % reacs.naux )
    for iaux in range(reacs.naux) :
        lines.append( '  integer, parameter  ::  iaux_%-8s = %2i  ! %s' % (reacs.aux_name[iaux],iaux+1,reacs.aux_desc[iaux]) )
    #endfor
    lines.append( '' )
    lines.append( '  ! cldadj(11) array of k1 adjustments for cloud cover (0-10 clear to overcast)' )
    lines.append( '  ! source: maul,1980' )
    lines.append( '  real, parameter :: cldadj(0:10) = &' )
    lines.append( '        (/1.0,1.00,1.00,0.79,0.75,0.72,0.68,0.62,0.53,0.41,0.35/)' )
    lines.append( '' )
    lines.append( '  ! molecular weights needed for VBS:' )
    lines.append( '  real  ::  mw_vbs(4) = (/ 250e6, 250e6, 150e6, 180e6 /) ! mw"s for poa, sisoa, asoa, bsoa, resp.')    
    lines.append( '' )
    lines.append( '  ! reaction indices:' )
    for ireac in range(reacs.nreac) :
        # by default parametes, except for the specials which are now used
        # in explicitly implemented chemistry instead of generated:
        attr = ', parameter'
        if reacs.label[ireac].startswith( 'M7_') : attr = '           '
        if reacs.label[ireac].startswith('VBS_') : attr = '           '
        # write:
        lines.append( '  integer%s  ::  ireac_%-8s = %2i  ! %s  ->  %s' % \
                (attr,reacs.label[ireac],ireac+1,reacs.reactants[ireac],reacs.products[ireac].strip()) )
    #endfor
    # need some parameters for M7, define a dummy if necessary:
    for lab in ['VBS_R63','VBS_R65','VBS_R66','M7_R71','M7_R72'] :
        if lab not in reacs.label :
            lines.append( '  integer             ::  ireac_%-8s = -1  ! dummy for hardcoded chemistry' % lab )
        #endif
    #endfor
    #
    lines.append( '' )
    if (labeling ) :
        counter = 0
        for ispec in range(specs.nspec) :
            for irspec in range(specs.nspec) :
                if ( len(reacs.array_prod[irspec][ispec]) > 0 ) :
                    counter = counter + 1
                    lines.append( '  integer,parameter  :: %s_to_%s = %2i ' % (specs.name[irspec], specs.name[ispec], counter ) )
                #endif
            #endfor
        #endfor
        lines.append( '' )
        lines.append( ' integer, parameter      ::  n_nonzeros = %i ' %(counter) )
        lines.append( ' integer, allocatable    ::  array_column_indices(:) ' )
        lines.append( ' integer, allocatable    ::  array_row_indices(:) ' )
    #endif
    lines.append( '' )
    lines.append( '' )
    lines.append( 'contains' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  ! ====================================================================' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  subroutine LE_Chem_Work_Init( status )' )
    if ( labeling ) :
        lines.append( '    use Indices ' )
        lines.append( '    use SA_Labeling, only : SA_Chem_Gas_Setup ' )
    #endif
    lines.append( '' )
    lines.append( '    ! --- in/out ---------------------------------' )
    lines.append( '' )
    lines.append( '    integer, intent(out)                ::  status' )
    lines.append( '' )
    lines.append( '    ! --- const ----------------------------------' )
    lines.append( '' )
    lines.append( '    character(len=*), parameter ::  rname = mname//\'/LE_Chem_Work_Init\'' )
    lines.append( '' )
    lines.append( '    ! --- local ----------------------------------' )
    lines.append( '' )
    lines.append( '    ! --- begin ----------------------------------' )
    lines.append( '' )
    if ( labeling ) :
        lines.append( '    allocate(array_column_indices(n_nonzeros), stat=status) ')
        lines.append( '    IF_NOTOK_RETURN(status=1) ' )
        lines.append( '' )
        lines.append( '    allocate(array_row_indices(n_nonzeros), stat=status) ')
        lines.append( '    IF_NOTOK_RETURN(status=1) ')
        lines.append( '' )
        for ispec in range(specs.nspec) :
            for irspec in range(specs.nspec) :
                if (len(reacs.array_prod[irspec][ispec] ) > 0 ) :
                    lines.append( '    array_column_indices(%s_to_%s) = i_%s ' % (specs.name[irspec], specs.name[ispec], specs.name[irspec] ) )
                    lines.append( '    array_row_indices(%s_to_%s) = i_%s ' % (specs.name[irspec], specs.name[ispec], specs.name[ispec] ) )
                #endif
            #endfor
        #endfor
        lines.append( '' )
        lines.append( '    call SA_Chem_Gas_Setup(array_column_indices, array_row_indices, n_nonzeros, status) ' )
        lines.append( '    IF_NOTOK_RETURN(status=1) ' )
        lines.append( '' )
    #endif
    lines.append( '    ! ok' )
    lines.append( '    status = 0' )
    lines.append( '' )
    lines.append( '  end subroutine LE_Chem_Work_Init' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  ! ***' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  subroutine LE_Chem_Work_Done( status )' )
    lines.append( '' )
    lines.append( '    ! --- in/out ---------------------------------' )
    lines.append( '' )
    lines.append( '    integer, intent(out)                ::  status' )
    lines.append( '' )
    lines.append( '    ! --- const ----------------------------------' )
    lines.append( '' )
    lines.append( '    character(len=*), parameter ::  rname = mname//\'/LE_Chem_Work_Done\'' )
    lines.append( '' )
    lines.append( '    ! --- begin ----------------------------------' )
    lines.append( '' )
    if (labeling ) :
        lines.append( '    deallocate(array_column_indices) ' )
        lines.append( '    deallocate(array_row_indices) ' )
    #endif
    lines.append( '    ! ok' )
    lines.append( '    status = 0' )
    lines.append( '' )
    lines.append( '  end subroutine LE_Chem_Work_Done' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  ! ***' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  !' )
    lines.append( '  ! Reaction rate expressions.' )
    lines.append( '  !' )
    lines.append( '  !  A x <J_label>  k = A * j' )
    lines.append( '  !' )
    lines.append( '  !  A  @E          k = A * exp(-E/T)' )
    lines.append( '  !' )
    lines.append( '  !  A^B@E          k = A x (T/300)B x exp(-E/T)' )
    lines.append( '  !' )
    lines.append( '  !  k1 & k2 & F & n' )
    lines.append( '  !                         k1 [M]                          log10(k1[M]/k2) +2  -1' )
    lines.append( '  !                 k = [------------] F^G   ,   G = [ 1 + (---------------)  ]' )
    lines.append( '  !                      1 + k1[M]/k2                              n' )
    lines.append( '  !' )
    lines.append( '  !     NOTE: in original document, the formula for G included a power \'-2\' ;' )
    lines.append( '  !     Comparison with the TROE function used in CB99, and other references' )
    lines.append( '  !     (e.g. http://jpldataeval.jpl. nasa.gov/pdf/JPL_15_AllInOne.pdf ),' )
    lines.append( '  !     show that this power should be \'+2\' .' )
    lines.append( '  !     Ferd Sauter, RIVM, july 2009' )
    lines.append( '  !' )
    lines.append( '  !     NOTE: the \'log\' in the original document should be interpreted' )
    lines.append( '  !     as the 10-log, not the e-log .' )
    lines.append( '  !       Ferd Sauter, RIVM, july 2009' )
    lines.append( '  !     However, this seems only true if \'n=1\' is used.' )
    lines.append( '  !     With \'n=log(10)\' as seen in the original LOTOS-EUROS code' )
    lines.append( '  !     it should be \'log\' .' )
    lines.append( '  !       Arjo Segers, TNO, sept 2011' )
    lines.append( '  !' )
    lines.append( '  !  %2 k1 & k2 & k3' )
    lines.append( '  !                             k3[M]' )
    lines.append( '  !                 k = k1 + ------------' )
    lines.append( '  !                          1 + k3[M]/k2' )
    lines.append( '  !' )
    lines.append( '  !  %3 k1 & k2     k = k1 + k2[M]' )
    lines.append( '  !' )
    lines.append( '' )
    lines.append( '  elemental function rate1( k1, k2, F, n, M )' )
    lines.append( '' )
    lines.append( '    ! --- in/out ---------------------------------' )
    lines.append( '' )
    lines.append( '    real                      ::  rate1' )
    lines.append( '    real, intent(in)          ::  k1' )
    lines.append( '    real, intent(in)          ::  k2' )
    lines.append( '    real, intent(in)          ::  F' )
    lines.append( '    real, intent(in)          ::  n' )
    lines.append( '    real, intent(in)          ::  M' )
    lines.append( '' )
    lines.append( '    ! --- local ----------------------------------' )
    lines.append( '' )
    lines.append( '    real         ::  G' )
    lines.append( '' )
    lines.append( '    ! --- begin ----------------------------------' )
    lines.append( '' )
    lines.append( '    !! CB05 definition, usually n=1 :' )
    lines.append( '    !G = 1.0 / ( 1.0 + ( log10(k1*M/k2)/n )**2 )' )
    lines.append( '' )
    lines.append( '    ! LOTOS-EUROS definition, usually n=log(10) :' )
    lines.append( '    G = 1.0 / ( 1.0 + ( log(k1*M/k2)/n )**2 )' )
    lines.append( '' )
    lines.append( '    rate1 = ( (k1*M)/(1.0+k1*M/k2) ) * F**G' )
    lines.append( '' )
    lines.append( '  end function rate1' )
    lines.append( '' )
    lines.append( '  ! *' )
    lines.append( '' )
    lines.append( '  elemental function rate2( k1, k2, k3, M )' )
    lines.append( '' )
    lines.append( '    ! --- in/out ---------------------------------' )
    lines.append( '' )
    lines.append( '    real                      ::  rate2' )
    lines.append( '    real, intent(in)          ::  k1' )
    lines.append( '    real, intent(in)          ::  k2' )
    lines.append( '    real, intent(in)          ::  k3' )
    lines.append( '    real, intent(in)          ::  M' )
    lines.append( '' )
    lines.append( '    ! --- begin ----------------------------------' )
    lines.append( '' )
    lines.append( '    rate2 = k1 + (k3*M)/(1.0+k3*M/k2)' )
    lines.append( '' )
    lines.append( '  end function rate2' )
    lines.append( '' )
    lines.append( '  ! *' )
    lines.append( '' )
    lines.append( '  elemental function rate3( k1, k2, M )' )
    lines.append( '' )
    lines.append( '    ! --- in/out ---------------------------------' )
    lines.append( '' )
    lines.append( '    real                      ::  rate3' )
    lines.append( '    real, intent(in)          ::  k1' )
    lines.append( '    real, intent(in)          ::  k2' )
    lines.append( '    real, intent(in)          ::  M' )
    lines.append( '' )
    lines.append( '    ! --- begin ----------------------------------' )
    lines.append( '' )
    lines.append( '    rate3 = k1 + k2*M' )
    lines.append( '' )
    lines.append( '  end function rate3' )
    lines.append( '' )
    lines.append( '  ! *' )
    lines.append( '' )
    lines.append( '  elemental function phux( x,y,z, zen )' )
    lines.append( '' )
    lines.append( '    ! --- in/out ---------------------------------' )
    lines.append( '' )
    lines.append( '    real                      ::  phux' )
    lines.append( '    real, intent(in)          ::  x, y, z' )
    lines.append( '    real, intent(in)          ::  zen' )
    lines.append( '' )
    lines.append( '    ! --- local ----------------------------------' )
    lines.append( '' )
    lines.append( '    real    ::  arg' )
    lines.append( '    real    ::  fac' )
    lines.append( '' )
    lines.append( '    ! --- begin ----------------------------------' )
    lines.append( '' )
    lines.append( '    arg = zen*z' )
    lines.append( '    if ( arg < 1.57 ) then' )
    lines.append( '      fac = max( -10.0, y*(1.0-1.0/cos(arg)) )' )
    lines.append( '      phux = x*exp(fac)' )
    lines.append( '    else' )
    lines.append( '      phux = 0.0' )
    lines.append( '    end if' )
    lines.append( '' )
    lines.append( '  end function phux' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  ! ====================================================================' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  ! computation of the photochemical reaction rates' )
    lines.append( '  ! and the temperature dependent rates' )
    lines.append( '  ! unit: ppb**(-n) min**(-1)' )
    lines.append( '  !   with \'n\' the order of reaction, thus n=1,2, ..' )
    lines.append( '' )
    lines.append( '  subroutine LE_Chem_Work_Rates( rk, aux, rk_het )' )
    lines.append( '' )
#    lines.append( '    use Binas, only : Avog, Rgas' )
#    lines.append( '    use dims, only : nspec, runF' )
    lines.append( '    use LE_Chem_Hetro, only : nreac_het' )
    lines.append( '    use LE_Chem_Hetro, only : ireac_N2O5_NH4HSO4a_f' )
    lines.append( '    use LE_Chem_Hetro, only : ireac_N2O5_ss_f' )
    lines.append( '    use LE_Chem_Hetro, only : ireac_N2O5_ss_c' )
    lines.append( '    use LE_Chem_Hetro, only : ireac_HNO3_ss_f' )
    lines.append( '    use LE_Chem_Hetro, only : ireac_HNO3_ss_c' )
    lines.append( '' )
    lines.append( '    ! --- in/out ---------------------------------' )
    lines.append( '' )
    lines.append( '    real(8), intent(out)       ::  rk(nreac)         ! reaction rates (1/ppb**n/min)' )
    lines.append( '    real(8), intent(inout)     ::  aux(naux)         ! auxilary values' )
    lines.append( '    real(8), intent(in)        ::  rk_het(nreac_het) ! hetr. reac. rates  (1/ppb**n/min)' )
    lines.append( '' )
    lines.append( '    ! --- const ----------------------------------' )
    lines.append( '' )
    lines.append( '    ! Constant used for 10-based logarithms in reaction rates.' )
    lines.append( '    ! Should be replaced by the 10-log available in Fortran:' )
    lines.append( '    !    log(x)/ln10 = log(x)/log(10) = log10(x)' )
    lines.append( '    real, parameter   ::  ln10 = 2.3025851' )
    lines.append( '' )
    lines.append( '    ! --- local ----------------------------------' )
    lines.append( '' )
    lines.append( '    ! auxilary values' )
    lines.append( '    real              ::  air' )
    lines.append( '    real              ::  ch2o' )
    lines.append( '    real              ::  T' )
    lines.append( '    real              ::  zen' )
    lines.append( '    real              ::  cldfac' )
    lines.append( '    real              ::  cldsulf' )
    lines.append( '    real              ::  rh' )
    lines.append( '    real              ::  invT' )
    lines.append( '    real              ::  Tdiv300' )
    lines.append( '' )
    lines.append( '    ! --- begin ----------------------------------' )
    lines.append( '' )
    lines.append( '    ! auxilary concentrations:' )
    lines.append( '    air     = aux(iaux_air)      ! (mlc/cm3)' )
    lines.append( '    ch2o    = aux(iaux_H2O)      ! water concentration' )
    lines.append( '' )
    lines.append( '    ! extract auxilary values:' )
    lines.append( '    T       = aux(iaux_T     )   ! temperature (K)' )
    lines.append( '    zen     = aux(iaux_zen   )   ! solar zenith angle ?' )
    lines.append( '    cldfac  = aux(iaux_cldfac)   ! cloud cover factor applied to photolysis rates' )
    lines.append( '    cldsulf = aux(iaux_cldsulf)  ! cloud cover factor for sulf reaction' )
    lines.append( '    rh      = aux(iaux_rh    )   ! relative humidity at 2m' )
    lines.append( '' )
    lines.append( '    ! precompute inversions as a much faster multiplication ...' )
    lines.append( '    invT    = 1.0/T' )
    lines.append( '    Tdiv300 = T/300.0' )

    # loop over reactions (0,..,nreac-1) :
    ireac = 0
    while ireac < reacs.nreac :

        # current:
        rate_expr = reacs.rate[ireac].strip()
        
        lines.append( '' )                       
        # add reaction as comment:
        lines.append( '    ! %s : %s  ->  %s' % (reacs.label[ireac],reacs.reactants[ireac],reacs.products[ireac].strip()) )
        # add rate expression too:
        lines.append( '    ! rate expr. : %s' % rate_expr )

        #
        # ~ rate expr. :       8.3E-5    $      2.0               $              0.1   $    90.0
        #   expand to  :   k = 8.3E-5 * (1.0 +  2.0 * cldsulf) * max(1.0, 1.0 +  0.1 *(rh - 90.0) )
        if '$' in rate_expr :
            # split:
            A,B,C,D = rate_expr.split('$')
            # clean up :
            A = A.strip()
            B = B.strip()
            C = C.strip()
            D = D.strip()
            # form expression:
            rk = '%s * (1.0 + %s * cldsulf) * max(1.0, 1.0 + %s*(Rh - %s) )' % (A,B,C,D)

        #
        # ~ rate expr.   :      k1  &  A x <J_label>
        #   expand to    :  k = k1  +  A * J
        #
        elif ('&' in rate_expr) and ('x' in rate_expr) :
            # thermal + photolysis rate
            T,P = rate_expr.split('&')
            # photolysis rate
            A,J = P.split('x')
            # clean up :
            T = T.strip()
            A = A.strip()
            J = J.strip()
            # test content of first part:
            if '@' in T :
                # default temperature dependency:
                rk = cb05_rk( T )
                # specials ...
                if J == '<N2O5_IUPAC05>' :
                    rk = rk + ' + phux( 3.79e-5, 1.70537, 0.80153, zen ) * cldfac'
                else :
                    logger.error( 'unsupported photolysis label "%s"' % J )
                    logger.error( '  reaction label : %s' % reacs.label[ireac] )
                    logger.error( '  rate expr.     : %s' % rate_expr )
                    raise Exception
                #endif
            else :
                logger.error( 'do not know what to do with first part "%s"' % T )
                logger.error( '  reaction label : %s' % reacs.label[ireac] )
                logger.error( '  rate expr.     : %s' % rate_expr )
                raise Exception
            #endif

        #
        # ~ rate expr   :   k1 & k2 &  F  & n &  C   /  kd
        #   defaults                  0.3   1 & 1.0
        #   expand to   :   rate1(k1,k2,F,n) * C / kd
        #
        elif rate_expr.find('/') >=0 :    
            O, R = rate_expr.split('/')
            # cleanup
            O.strip()
            R.strip()
            # first part:
            rk = ''
            if '&' in O :
                expr = O.split('&')
                rk1 = cb05_rk(expr[0])
                rk2 = cb05_rk(expr[1])
                rF = '0.3'
                if len(expr) >= 3 : rF = expr[2].strip()
                rn = '1.0'
                if len(expr) >= 4 : rn = expr[3].strip()
                rk = 'rate1( %s, %s, %s, %s, air )' % (rk1,rk2,rF,rn) 
                if len(expr) >= 5 : 
                    rp = expr[4].strip()
                    rk = rp + '*' + rk
                #endif
            elif O.find('*') :
                rk = O
            else :
                logger.error( 'could not resolve first expression "%s" from "%s"' % (O,rate_expr) )
                raise Exception
            #endif
            # second part:
            if '@' in R :
                rk = rk + ' / &\n             (' + cb05_rk(R) + ')'
            else :
                logger.error( 'could not resolve second expression "%s" from "%s"' % (R,rate_expr) )
                raise Exception
            #endif
                       
        # ~ rate expr   :   A x <J_label>
        #   exapand to  :   k = A * J
        elif rate_expr.find('x') >= 0 :
            # photolysis rate
            A,J = rate_expr.split('x')
            # cleanup:
            A = A.strip()
            J = J.strip()
            # temporary ...
            #rk = '%s * J_%s' % (A,J)
            # hardcoded values ...
            rk = ''
            if J == '<NO2_SAPRC99>'         : rk = str(A) + ' * phux( 1.07e-2, 1.01319, 0.83330, zen ) * cldfac'
            if J == '<O3_O3P_IUPAC05>'      : rk = 'phux( 5.36e-4, 0.34764, 0.91030, zen ) * cldfac'
            if J == '<O3_O1D_IUPAC05>'      : rk = 'phux( 3.22e-5, 4.45037, 0.78028, zen ) * cldfac'
            if J == '<NO3NO2_SAPRC99>'      : rk = 'phux( 2.73e-1, 0.29327, 0.92401, zen ) * cldfac'
            if J == '<NO3NO_SAPRC99>'       : rk = 'phux( 2.74e-2, 0.26226, 0.92849, zen ) * cldfac'
            if J == '<HNO2_IUPAC05>'        : rk = 'phux( 8.96e-4, 0.99438, 0.83295, zen ) * cldfac'
            if J == '<H2O2_SAPRC99>'        : rk = 'phux( 7.78e-6, 1.91463, 0.79810, zen ) * cldfac'
            if J == '<HO2NO2_IUPAC05>'      : rk = '0.0 !; stop \'no photolysis rate for HO2->NO2\''
            if J == '<HNO3_IUPAC05>'        : rk = 'phux( 5.48e-7, 2.86922, 0.79561, zen ) * cldfac'
            if J == '<N2O5_IUPAC05>'        : rk = 'phux( 3.79e-5, 1.70537, 0.80153, zen ) * cldfac'
            if J == '<NTR_IUPAC05>'         : rk = '0.0 !; stop \'no photolysis rate for NTR\''
            if J == '<COOH_SAPRC99>'        : rk = '0.0 !; stop \'no photolysis rate for COOH\''
            if J == '<HCHO_R_SAPRC99>'      : rk = str(A) + ' * phux( 4.05e-5, 2.06917, 0.80267, zen ) * cldfac'
            if J == '<HCHO_M_SAPRC99>'      : rk = 'phux( 4.92e-5, 1.60973, 0.80184, zen ) * cldfac'
            if J == '<CCHO_R_SAPRC99>'      : rk = 'phux( 5.40e-6, 2.52915, 0.79722, zen ) * cldfac'
            if J == '<PAN_IUPAC05>'         : rk = '0.0 !; stop \'no photolysis rate for PAN\''
            if J == '<C2CHO_SAPRC99>'       : rk = '0.0 !; stop \'no photolysis rate for C2CHO\''
            #
            #if J == '<MGLY_IUPAC05>'        : rk = '0.0 ; stop \'no photolysis rate for MGLY; 0.02*J[NO2] ?\''
            # JW: J[mgly] you can use 5.5*J[bCH2O] channel as in TM.
            # For the moment, use the 0.02*J[NO2] ass suggested:
            #if J == '<MGLY_IUPAC05>'        : rk = str(A) + ' * phux( 1.07e-2, 1.01319, 0.83330, zen ) * cldfac'
            #
            #if J == '<ACROLEIN_SAPRC99>'    : rk = '0.0 !; stop \'no photolysis rate for ACROLEIN\''
            # 1.7e-5 times JNO2 (J== <NO2_SAPRC99>
            #if J == '<ACROLEIN_SAPRC99>'    : rk = str(A) + ' * phux( 1.07e-2, 1.01319, 0.83330, zen ) * cldfac'  
            #
            ## adhoc ...
            #rk = 'aux(iaux_Jfac) * '+rk
            # check ..
            if len(rk) == 0 :
                print( 'ERROR - unsupported photolysis expression : ', J )
                raise
            #endif

        #
        #    k1 & k2 & F & n
        #                           k1 [M]                          log10(k1[M]/k2) +2  -1
        #                   k = [------------] F^G   ,   G = [ 1 + (---------------)   ]
        #                        1 + k1[M]/k2                               n         
        #  
        #        NOTE: in the original document, the formula for G included a power '-2' ;
        #        Comparison with the TROE function used in CB99, and other references
        #        (e.g. http://jpldataeval.jpl. nasa.gov/pdf/JPL_15_AllInOne.pdf ), 
        #        show that this power should be '+2' .
        #          Ferd Sauter, RIVM, july 2009
        #
        #        NOTE: the 'log' in the original document should be interpreted
        #        as the 10-log, not the e-log .
        #          Ferd Sauter, RIVM, july 2009
        #        However, this seems only true if 'n=1' is used.
        #        With 'n=log(10)' as seen in the original LOTOS-EUROS code
        #        it should be 'log' .
        #
        #    %2 k1 & k2 & k3
        #                               k3[M]
        #                   k = k1 + ------------
        #                            1 + k3[M]/k2
        #  
        #    %3 k1 & k2     k = k1 + k2[M]
        #
        elif rate_expr.find('&') >= 0 :
            # select form:
            if rate_expr.startswith('%3') :
                expr = rate_expr[2:].split('&')
                rk1 = cb05_rk(expr[0])
                rk2 = cb05_rk(expr[1])
                rk = 'rate3( %s, %s, air )' % (rk1,rk2)
            elif rate_expr.startswith('%2') :
                expr = rate_expr[2:].split('&')
                rk1 = cb05_rk(expr[0])
                rk2 = cb05_rk(expr[1])
                rk3 = cb05_rk(expr[2])
                rk = 'rate2( %s, %s, %s, air )' % (rk1,rk2,rk3)
            else :
                expr = rate_expr.split('&')
                rk1 = cb05_rk(expr[0])
                rk2 = cb05_rk(expr[1])
                rF = '0.3'
                if len(expr) >= 3 : rF = expr[2].strip()
                rn = '1.0'
                if len(expr) >= 4 : rn = expr[3].strip()
                rk = 'rate1( %s, %s, %s, %s, air )' % (rk1,rk2,rF,rn)
                if len(expr) >= 5 : 
                    rp = expr[4].strip()
                    rk = rp + '*' + rk
                #endif
            #endif

        #
        # ~ rate expr  :   A^B@E
        #   expand to  :   k = A x (T/300)B x exp(-E/T)
        #  
        else :
            # default expression:
            rk = cb05_rk( rate_expr )

        #endif

        # code to set reactions rate:
        lines.append( '    rk(ireac_%s) = %s' % (reacs.label[ireac],rk) )
        
        ireac = ireac + 1
        # testing ...
        #break
        
    #endwhile
    
    # add dummy assignment ?
    if ireac == 0 :
        lines.append( '' )
        lines.append( '    ! no reactions, dummy assignment:' )
        lines.append( '    rk = 0.0' )
    #endif
    
    # write closing:
    lines.append( '' )
    lines.append( '    !' )
    lines.append( '    ! Set some extra auxilary values, these are results of lumping reactions' )
    lines.append( '    !' )
    lines.append( '    ! For example the reaction:' )
    lines.append( '    !       O + H2O   --> 2 OH' )
    lines.append( '    ! is involved in the photolysis of ozone:' )
    lines.append( '    !       O3        --> O2 + O' )
    lines.append( '    !' )
    lines.append( '    ! Further the reactions:' )
    lines.append( '    !       ROR       --> 1.10 ALD + 0.96 XO2 + 0.94 HO2 +' )
    lines.append( '    !                     0.04 XO2N + 0.02 ROR - 2.10 PAR' )
    lines.append( '    !' )
    lines.append( '    !       ROR       --> HO2' )
    lines.append( '    !' )
    lines.append( '    !       ROR + NO2 -->' )
    lines.append( '    !' )
    lines.append( '    ! are involved in the reaction:' )
    lines.append( '    !       PAR + OH  --> 0.87 XO2 + 0.13 XO2N + 0.11 HO2 +' )
    lines.append( '    !                     0.11 ALD + 0.76 ROR  - 0.11 PAR' )
    lines.append( '    !' )
    lines.append( '    ! The fractions below are the result of the different reaction rates' )
    lines.append( '    ! (ratios depends on temperature):' )
    lines.append( '    !' )
    lines.append( '    ! rat 10 must be in order of mlc/cm3 thus multiplied with factor ppm_to_mlccm3  ')
    lines.append( '    ! This factor is 2.4621e13 with st temp (298K ) and st pres (1.013hPa)' )
    lines.append( '    aux(iaux_rat10  ) = 42900.0 * 1.e+06 * aux(iaux_ppm_to_mlccm3)' )
    lines.append( '    aux(iaux_alpha8 ) = cH2O * 326000.0 / ( cH2O * 326000.0 + aux(iaux_rat10) )' )
    lines.append( '    aux(iaux_rat75  ) = 90000.0 * exp( 7000.0*(1/298.0-1.0/T) )' )
    lines.append( '    aux(iaux_alpha10) = aux(iaux_rat75) / ( aux(iaux_rat75) + 390000.0 )' )
    lines.append( '' )
    lines.append( '  end subroutine LE_Chem_Work_Rates' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  ! ***' )
    lines.append( '' )
    lines.append( '' )
    lines.append( '  pure subroutine LE_Chem_Work_Iter( n, y, ysum, gdt, yp, yl, Q, &' )
    if (labeling ) :
        lines.append( '                                       n_nonzeros, array_prod, &')
    #endif
    if (with_vbs ) :
        lines.append( '                                       drog_loss, drog_prod_neg, &')
    #endif
    lines.append( '                                       nreac, rk, naux, aux )')
    lines.append( '' )
    lines.append( '    use Indices' )
    lines.append( '' )
    lines.append( '    ! --- in/out ---------------------------------' )
    lines.append( '' )
    lines.append( '    integer, intent(in)     ::  n             ! number of concentrations' )
    lines.append( '    real, intent(inout)     ::  y(1:n)        ! concentrations  (conc)' )
    lines.append( '    real, intent(in)        ::  ysum(1:n)     ! work array' )
    lines.append( '    real, intent(in)        ::  gdt           ! sec' )
    lines.append( '    real, intent(out)       ::  yp(1:n)       ! production rate (conc/s)' )
    lines.append( '    real, intent(out)       ::  yl(1:n)       ! loss rate (conc/s)' )
    lines.append( '    real, intent(in)        ::  Q(1:n)        ! source term (conc/s)' )
    if ( labeling ) :
        lines.append( '    integer, intent(in)     ::  n_nonzeros' )
        lines.append( '    real, intent(out)       ::  array_prod(n_nonzeros) ' )
    #endif
    if ( with_vbs ) :
        lines.append( '    real, intent(out)       ::  drog_loss(1:n)     ! loss rate of soa precursors ')
        lines.append( '    real, intent(out)       ::  drog_prod_neg(1:n) ! negative production rate of soa precursors ')
    #endif
    lines.append( '    integer, intent(in)     ::  nreac         ! number of reactions' )
    lines.append( '    real(8), intent(in)     ::  rk(nreac)     ! reaction rate coeff.' )
    lines.append( '    integer, intent(in)     ::  naux          ! number of auxilary values' )
    lines.append( '    real(8), intent(in)     ::  aux(naux)     ! auxilary values' )
    lines.append( '' )
    lines.append( '    ! --- local ----------------------------------' )
    lines.append( '' )
    lines.append( '    real                    :: alpha8' )
    lines.append( '    real                    :: alpha9' )
    lines.append( '    real                    :: alpha10' )
#    lines.append( '    real                    :: rat10' )
    lines.append( '    real                    :: rat75' )
    lines.append( '' )
    ispec = 0    
    while ispec < (specs.nspec) :
        if (labeling) :
            irspec = 0            
            while irspec < (specs.nspec)  :
                if ( 'ftmp' in reacs.array_prod[irspec][ispec] or 'ftmp' in reacs.yp[ispec] or 'ftmp' in reacs.yl[ispec] ) :
                    lines.append( '    real                    :: ftmp_%s ' %(specs.name[ispec] ) )
                    irspec = specs.nspec
                #endif
                irspec = irspec + 1
            #endwhile
        else :
            if ('ftmp' in reacs.yp[ispec] ) :
                lines.append( '    real                    ::  ftmp_%s ' %(specs.name[ispec]) )
            elif ('ftmp' in reacs.yl[ispec] ) :
                lines.append( '    real                    ::  ftmp_%s ' %(specs.name[ispec]) )
            #endif
        #endif               
        ispec = ispec + 1
    #endwhile
    lines.append( '' )
    if ( with_vbs ) :
        #lines.append( '    integer                 ::  n_vbs_bins' ) # get from Indices!
        lines.append( '    integer                 ::  s,v' )
        lines.append( '    real, parameter  ::  conv_par_alkanes = 10.0' )
        lines.append( '    real, parameter  ::  conv_ole_alkenes = 10.0' )
        lines.append( '    real, parameter  ::  mw_soa_prec(n_soa_prec) = (/ 30e6, 28e6, 92e6, 106e6, 68e6, 136e6 /) ')
        lines.append( '    real       ::  alpha_par(n_vbs_bins)' )
        lines.append( '    real       ::  alpha_ole(n_vbs_bins)' )
        lines.append( '    real       ::  alpha_tol(n_vbs_bins)' )
        lines.append( '    real       ::  alpha_xyl(n_vbs_bins)' )
        lines.append( '    real       ::  alpha_iso(n_vbs_bins)' )
        lines.append( '    real       ::  alpha_terp(n_vbs_bins)' )
        lines.append( '    real       ::  soa_prec_reacted' )
        lines.append( '    real(8)    ::  rk_ro2_no                   ! reaction rate for RO2 radicals with NO' )
        lines.append( '    real(8)    ::  rk_ro2_ho2                  ! reaction rate for RO2 radicals with HO2' )
        lines.append( '    real(8)    ::  rk_ro2_ro2(n_soa_prec)      ! reaction rate for RO2 radicals with themselves' )
        lines.append( '    real(8)    ::  rk_ro2_tot(n_soa_prec)      ! sum of reaction rates for RO2 radicals' )
        lines.append( '    real       ::  branchratio(n_soa_prec)' )
        lines.append( '    real       ::  mass_yields_eff(n_soa_prec, n_vbs_bins) ! 2D array containing effective (i.e. NOx dependent) mass yields for each SOA precursor' )
        lines.append( '' )    
        lines.append( '    real, parameter     ::  mass_yields_highnox(n_soa_prec, n_vbs_bins) = &  ! 2D array containing mass yields for each soa precursor' )
        lines.append( '    !rows:    par     ole     tol     xyl    iso     terp  ' )   
        lines.append( '    reshape((/0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.0  , & ! column vbs 1 ' )  
        lines.append( '              0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.0  , & ! column vbs 2 ' )
        lines.append( '              0.0   , 0.002 , 0.003 , 0.002, 0.001 , 0.012, & ! column vbs 3 ' )
        lines.append( '              0.094 , 0.0155, 0.165 , 0.195, 0.023 , 0.122, & ! column vbs 4 ' )
        lines.append( '              0.0   , 0.0605, 0.300 , 0.300, 0.015 , 0.201, & ! column vbs 5 ' )
        lines.append( '              0.0   , 0.210 , 0.435 , 0.435, 0.0   , 0.500, & ! column vbs 6 ' )
        lines.append( '              0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.0  , & ! column vbs 7 ' )
        lines.append( '              0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.0  , & ! column vbs 8 ' )
        lines.append( '              0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.0  /), & ! column vbs 9 ' )
        lines.append( '    shape(mass_yields_highnox)) ' )
        lines.append( '' )    
        lines.append( '    real, parameter     ::  mass_yields_lownox(n_soa_prec, n_vbs_bins) = &  ! 2D array containing mass yields for each soa precursor' )
        lines.append( '    !rows:    par     ole     tol     xyl    iso     terp   ' )         
        lines.append( '    reshape((/0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.0  , & ! column vbs 1' )        
        lines.append( '              0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.000, & ! column vbs 2' )     
        lines.append( '              0.0   , 0.014 , 0.075 , 0.075, 0.009 , 0.107, & ! column vbs 3' )  
        lines.append( '              0.1875, 0.0265, 0.225 , 0.300, 0.030 , 0.092, & ! column vbs 4' )
        lines.append( '              0.0   , 0.0945, 0.375 , 0.375, 0.015 , 0.359, & ! column vbs 5' )
        lines.append( '              0.0   , 0.300 , 0.525 , 0.525, 0.0   , 0.600, & ! column vbs 6' )
        lines.append( '              0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.0  , & ! column vbs 7' )
        lines.append( '              0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.0  , & ! column vbs 8' )
        lines.append( '              0.0   , 0.0   , 0.0   , 0.0  , 0.0   , 0.0  /), & ! column vbs 9' )
        lines.append( '    shape(mass_yields_lownox))    ' )

    #endif
    
    lines.append( '' )
    lines.append( '    ! --- begin ----------------------------------' )
    lines.append( '' )
    if ( labeling ) :
        lines.append( '    array_prod = 0.0 ' )
        lines.append( '' )
    #endif
    #
    if ( with_vbs ) :
        lines.append( '    drog_loss = 0.0 ' )
        lines.append( '    drog_prod_neg = 0.0 ' )
        lines.append( '    mass_yields_eff = 0.0 ' )
        lines.append( '    rk_ro2_no   = rk(ireac_VBS_R63)*y(i_NO)' )
        lines.append( '    rk_ro2_ho2  = rk(ireac_VBS_R66)*y(i_HO2)' )
        lines.append( '    ' )
        lines.append( '    !for spec in vbs_precursor_species !order prec in drog: potolxiterp' )
        lines.append( '    do s=1, n_soa_prec' )
        lines.append( '    ' )
        lines.append( '        !rk_ro2_ro2(s)  = rk(ireac_VBS_R65)*soa_prec_reacted(s)' )
        lines.append( '        rk_ro2_ro2(s)  = 0.0 ! for now assumed to be negligible, as in Pye et al., 2010' )
        lines.append( '        rk_ro2_tot(s)  = rk_ro2_no + rk_ro2_ro2(s) + rk_ro2_ho2' )
        lines.append( '    ' )
        lines.append( '        if (rk_ro2_tot(s) > 0.0)  then' )
        lines.append( '            branchratio(s) = real( rk_ro2_no / rk_ro2_tot(s) )' )
        lines.append( '        else' )
        lines.append( '            branchratio(s) = 0.5' )
        lines.append( '        endif' )
        lines.append( '    ' )
        lines.append( '        mass_yields_eff(s,:) = branchratio(s) * mass_yields_highnox(s,:) + (1.0 - branchratio(s)) * mass_yields_lownox(s,:) ! 2nd dimension is n_vbs_bins for all of which branching ratio is equal' )
        lines.append( '    enddo' )
        lines.append( '    ' )
        lines.append( '    ! convert mass yield to ppb yield and assign to precursor species' )
        lines.append( '    do v=1, n_vbs_bins' )
        lines.append( '      alpha_par(v) = mass_yields_eff(1,v) * (mw_soa_prec(1) / mw_vbs(3)) / conv_par_alkanes' )
        lines.append( '      alpha_ole(v) = mass_yields_eff(2,v) * (mw_soa_prec(2) / mw_vbs(3)) / conv_ole_alkenes' )
        lines.append( '      alpha_tol(v) = mass_yields_eff(3,v) * (mw_soa_prec(3) / mw_vbs(3)) ' )
        lines.append( '      alpha_xyl(v) = mass_yields_eff(4,v) * (mw_soa_prec(4) / mw_vbs(3))' )
        lines.append( '      alpha_iso(v) = mass_yields_eff(5,v) * (mw_soa_prec(5) / mw_vbs(4))' )
        lines.append( '      alpha_terp(v) = mass_yields_eff(6,v) * (mw_soa_prec(6) / mw_vbs(4))' )
        lines.append( '    enddo' )

    #endif
    #
    # factors which are need for lumping the reactions, at this moment
    # the reaction O + H2O --> 2OH is set into the ozone photlysis reaction
    # and further the reaction of ROR is set into the reaction PAR + OH -->
    # which has ROR as a product.
    # For now these factors are hardcoded is this generating program    
    
    # create lists with length number of specs and filled with blanks
    alpha8  = []
    alpha9  = []
    alpha10 = []
    rat75   = []
    ftmp    = []
    
    for ispec in range(specs.nspec) :
        alpha8.append('')
        alpha9.append('')
        alpha10.append('')
        rat75.append('')
        ftmp.append('')
    #endfor
    
    for ispec in range(specs.nspec) :
        if specs.name[ispec] == 'NO2' :
            rat75[ispec] = 'aux(iaux_rat75)'           
            alpha9[ispec] = '22.*aux(iaux_ppb_to_mlccm3)/(y(ispec_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))'
            ftmp[ispec] = '0.8*alpha9'
            
        elif specs.name[ispec] == 'O3' :
            alpha8[ispec] = 'aux(iaux_alpha8)'
            ftmp[ispec] = '1.-alpha8'
            
        elif specs.name[ispec] == 'ALD' :
            rat75[ispec] = 'aux(iaux_rat75)'
            alpha9[ispec] = 'y(i_NO2)*22./(y(i_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))'
            alpha10[ispec] = 'aux(iaux_alpha10)'
            ftmp[ispec] = '0.13 + (0.8*(1.-alpha9)*(0.185+alpha10*1.115))'
            
        elif specs.name[ispec] == 'PAR' :
            rat75[ispec] = 'aux(iaux_rat75)'
            alpha9[ispec] = 'y(i_NO2)*22./(y(i_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))'
            alpha10[ispec] = 'aux(iaux_alpha10)'
            ftmp[ispec] = '0.13+(0.8*(1.-alpha9)*(0.415+alpha10*1.885))'
            
        elif specs.name[ispec] == 'OH' :
            alpha8[ispec] = 'aux(iaux_alpha8)'
            ftmp[ispec] = '2.*alpha8'
            
        elif specs.name[ispec] == 'HO2' :
            alpha9[ispec] = 'y(i_NO2)*22./(y(i_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))'
            ftmp[ispec] = '0.13 + (0.8*(1.-alpha9))'
            
        elif specs.name[ispec] == 'XO2' :
            rat75[ispec] = 'aux(iaux_rat75)' 
            alpha9[ispec] = 'y(i_NO2)*22./(y(i_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))'
            alpha10[ispec] = 'aux(iaux_alpha10)'
            ftmp[ispec] = '0.93 + (0.8*(1.-alpha9)*(0.615+alpha10*0.385))'
            
        #endif
    #endfor
    
    # flag to see if any reactions were added:
    with_reacs = False
            
    # loop over all tracers:
    for ispec in range(specs.nspec) :

        # tracer not enabled ? then skip:
        if not specs.enabled[ispec] : continue

        # no reactions for this tracer ? then skip:
        if len(reacs.yir[ispec]) == 0 : continue
        
        # reset flag:
        with_reacs = True
    
        # current name:
        specname = specs.name[ispec]
        
        # intro:
        lines.append( '' )
        
        #if 'soa' in specs.props[ispec] :
        #    lines.append( '   if (runF%do_secorg) then' )
        ##endif
        lines.append( '    !' )
        lines.append( '    ! %s' % specname )
        lines.append( '    !' )
        
        # add reactions as comment:
        lines.append( '' )
        for ir in reacs.yir[ispec] :
            lines.append( '    ! %s : %s  ->  %s' % 
               (reacs.label[ir],reacs.reactants[ir],reacs.products[ir].strip()) )
        #endfor
        
        if (labeling and specs.name[ispec].strip().lower() in labelled_specs ) :
            ftmp_bool = False
            irspec = 0
            while (irspec < specs.nspec ) :
                if ('ftmp' in reacs.array_prod[irspec][ispec]) or ('ftmp' in reacs.yl[ispec]) or ('ftmp' in reacs.yp[ispec] ) :
                    ftmp_bool = True
                    irspec = specs.nspec
                #endif
                irspec = irspec + 1
            #endwhile
            
            if (ftmp_bool) :

                lines.append('')
                if rat75[ispec] != '' :
                    lines.append('    rat75 = %s' %(rat75[ispec]))
                #endif
                if alpha8[ispec] != '' :
                    lines.append('    alpha8 = %s' %(alpha8[ispec]))
                #endif
                if alpha9[ispec] != '' :
                    lines.append('    alpha9 = %s' %(alpha9[ispec]))
                #endif
                if alpha10[ispec] != '':
                    lines.append('    alpha10 = %s' %(alpha10[ispec]))
                #endif
                if ftmp[ispec] != '':
                    lines.append( '    ftmp_%s = %s ' %(specname, ftmp[ispec]) )
                #endif
            #endif
            lines.append( '' )
            lines.append( '    yl(ispec_%s) = 0.0 %s ' %(specname, reacs.yl[ispec]) )
            lines.append( '' )
            counter = 0

            for irspec in range(specs.nspec) :
                if ( len(reacs.array_prod[irspec][ispec] ) > 0 ) :
                    lines.append( '    array_prod(%s_to_%s) = ' %( specs.name[irspec], specname) + reacs.array_prod[irspec][ispec] )
                    counter = counter + 1
                #endif
            #endfor
            lines.append( '' )
            if (counter == 0 ) :
                lines.append( '    yp(ispec_%s) = Q(ispec_%s) ' %(specname, specname) )
            else :
                lines.append( '    yp(ispec_%s) = Q(ispec_%s) %s & ' %(specname, specname, reacs.yp[ispec]) )
            #endif

            counter_2 = 0
            for irspec in range(specs.nspec) :
                if ( len(reacs.array_prod[irspec][ispec] ) > 0 ) :
                    if ( counter_2 < counter -1 ) :
                        lines.append( '       + array_prod(%s_to_%s) & ' %(specs.name[irspec], specname ) )
                    else :
                        lines.append( '       + array_prod(%s_to_%s)   ' %(specs.name[irspec], specname ) )
                    #endif
                    counter_2 = counter_2 + 1
                #endif
            #endfor
            lines.append( '' )
            lines.append( '    y(ispec_%s) = max( 0.0, (ysum(ispec_%s) + gdt*yp(ispec_%s) ) / ( 1.0 + gdt*yl(ispec_%s) ) )' %(specname,specname,specname,specname) )
        else :
            if ( ('ftmp' in reacs.yp[ispec]) or ('ftmp' in reacs.yl[ispec]) ) :
                lines.append('')
                if rat75[ispec] != '' :
                    lines.append('    rat75 = %s' %(rat75[ispec]))
                #endif
                if alpha8[ispec] != '' :
                    lines.append('    alpha8 = %s' %(alpha8[ispec]))
                #endif
                if alpha9[ispec] != '' :
                    lines.append('    alpha9 = %s' %(alpha9[ispec]))
                #endif
                if alpha10[ispec] != '':
                    lines.append('    alpha10 = %s' %(alpha10[ispec]))
                #endif
                if ftmp[ispec] != '':
                    lines.append( '    ftmp_%s = %s ' %(specname, ftmp[ispec]) )
                #endif
            #endif

            # add iteration if necessary:
            if (len(reacs.yl[ispec]) > 0) or (len(reacs.yp[ispec]) > 0) :
                lines.append( '' )
                lines.append( '    yl(ispec_%s) = 0.0 %s' % (specname,reacs.yl[ispec]) )
                lines.append( '' )
                lines.append( '    yp(ispec_%s) = Q(ispec_%s) %s' % (specname,specname,reacs.yp[ispec]) )
                lines.append( '' )
                lines.append( '    y(ispec_%s) = max( 0.0, ( ysum(ispec_%s) + gdt*yp(ispec_%s) ) / ( 1.0 + gdt*yl(ispec_%s) ) )' % (specname,specname,specname,specname) )
#                lines.append( '    y(ispec_%s) = ( ysum(ispec_%s) + gdt*yp(ispec_%s) ) / ( 1.0 + gdt*yl(ispec_%s) ) ' % (specname,specname,specname,specname) )
            #endif
            # soa precursor ?
            if ( with_vbs and ('soa_prec' in specs.props[ispec]) ) :
                if (len(reacs.drog_loss[ispec]) > 0 ) :
                    lines.append( '' )
                    lines.append( '    drog_loss(ispec_%s) = 0.0 + %s ' %( specname, reacs.drog_loss[ispec]) )
                    lines.append( '' )
                #endif
                if (len(reacs.drog_prod_neg[ispec]) > 0 ) :
                    lines.append( '' )
                    lines.append( '    drog_prod_neg(ispec_%s) = 0.0 + %s ' %(specname, reacs.drog_prod_neg[ispec]) )
                    lines.append( '' )
                #endif
            #endif
            
        #endif #if labeling
        #if 'soa' in specs.props[ispec] :
        #    lines.append( '   endif' )
        #    lines.append( '' )
        ##endif

    #endfor # specs loop
    
    # add something to prevent compilation errors ?
    if not with_reacs :
        lines.append( '    ! no production and loss reactions found:' )
        lines.append( '    yp = 0.0' )
        lines.append( '    yl = 0.0' )
        lines.append( '    ! copy input to output:' )
        lines.append( '    y = ysum' )
    #endif
    
    lines.append( '' )
    lines.append( '  end subroutine LE_Chem_Work_Iter' )
    lines.append( '' )
    lines.append( '' )
    lines.append( 'end module LE_Chem_Work' )

    # add newlines, split lines first at existing '\n' characters:
    xlines = []
    for line in lines :
        for subline in line.split('\n') :
            xlines.append( subline+'\n' )
        #endfor
    #endfor
    
    # replace the target file if the content has changed:
    utopya.UpdateTextFile( chemistry_module, xlines )
    
#enddef  # WriteSource_Work


# ***


def cb05_rk( rate_expr ) :

    """
    Form reaction rate out of expression:

       A ^ B @ E  : k = A (T/300)^B exp(-E/T)

    Some of A, B, and E might be missing.
    """

    # split into 'A^B' and 'E', later might be empty:
    #AB,sep,E = rate_expr.partition('@')
    ss = rate_expr.split('@')
    AB = ss[0]
    E = ''
    if len(ss) > 1 : E = ss[1]
    # split into 'A' and 'B', later might be empty:
    #A,sep,B = AB.partition('^')
    ss = AB.split('^')
    A = ss[0]
    B = ''
    if len(ss) > 1 : B = ss[1]
    # init with A part:
    A = A.strip()
    rk = A
    # add B part:
    B = B.strip()
    if len(B) > 0 : rk = '(' + rk + ' * (Tdiv300)**(%s))' % B
    # add E part; remove minus for negative E:
    E = E.strip()
    if len(E) > 0 :
        if E.startswith('-') :
            #rk = rk + ' * exp(%s*invT)' % E[1:]
            #rk = rk + ' * exp(%s/tk)' % E[1:]
            rk = rk + ' * exp(%s/T)' % E[1:]
        else :
            #rk = rk + ' * exp(-%s*invT)' % E
            #rk = rk + ' * exp(-%s/tk)' % E
            rk = rk + ' * exp(-%s/T)' % E
        #endif
    #endif  

    # ok
    return rk
    
#enddef  # cb05_rk
        

# ------------------------------------------------
# end
# ------------------------------------------------



