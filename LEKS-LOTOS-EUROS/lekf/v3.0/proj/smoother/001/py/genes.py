#! /usr/bin/env python

# ------------------------------------------------
# help
# ------------------------------------------------

"""
GeneS - Generate Sources

"""

# ------------------------------------------------
# common
# ------------------------------------------------

# get common stuff:
#   - logger instance
from genes_common import *


# ------------------------------------------------
# main class
# ------------------------------------------------

class GeneS( object ) :

    """
    Main class to generate a chemistry code.
    """
    
    def __init__( self, rcfile, basekey ) :
    
        """
        Create a chemistry code using settings in rcfile.
        Rcfile keys start with the provided basekey.
        """
        
        # modules:
        import os
        
        # tools:
        import rc
        import genes_tracers
        import genes_reactions
        import genes_model_le
        
        #
        # settings
        #
        
        # read settings:
        rcf = rc.RcFile( rcfile )
        
        ## enable labeling code ?
        #with_labelling = rcf.get( basekey+'.labelling', 'bool' )
        
        # storage for group lists:
        groups = { 'supported' : [], 'selected' : [] }
        # list with all supported grouperties:
        fields = rcf.get( basekey+'.group.supported' ).split()
        for field in fields : groups['supported'].append(field.strip())
        # grouperties to select tracers and reactions:
        fields = rcf.get( basekey+'.group.selected' ).split()
        for field in fields : groups['selected'].append(field.strip())

        # explicit cloud chemistry ?
        with_chem_cloud = rcf.get('le.chem.cloud','bool')
        # if not, then implicit:
        if not with_chem_cloud : groups['selected'].append('in-cloud')
        
        # check if selected groups are supported:
        for group in groups['selected'] :
            if group not in groups['supported'] :
                logger.error( 'selected group "%s" not one of supported : %s' % (group,str(groups['supported'])) )
                raise Exception
            #endif
        #endfor
        
        # storage for property lists:
        props = { 'supported' : [] }
        # list with all supported grouperties:
        fields = rcf.get( basekey+'.prop.supported' ).split()
        for field in fields : props['supported'].append(field.strip())
        
        # skip some tracers ?
        fields = rcf.get( basekey+'.tracers.skip' ).split()
        # cleanup and convert to lower case:
        skip_tracers = []
        for field in fields : skip_tracers.append(field.strip().lower())
        
        # labeling ?
        with_labeling = rcf.get(basekey+'.labeling', 'bool', default=False)
        labelled_specs = []
        if with_labeling :
            logger.info( 'Labeling is enabled ...' )
            # read labelled specs
            fields = rcf.get( 'labels.labelled.specs' ).split()
            for field in fields : labelled_specs.append(field.strip().lower())
            logger.info( ' Labelled species : %s ' %(str(labelled_specs)) )
            
        else :
            logger.info( ' Labeling is not enabled... ' )
        #endif
                    
        # with vbs?
        with_vbs = rcf.get(basekey+'.vbs', 'bool')
        
        if ( with_vbs and with_labeling ) :
            logger.error( ' ################################################## ' )
            logger.error( ' VBS code does not work together with labeling code ' )
            logger.error( ' ################################################## ' )
            raise Exception
        elif with_vbs :
            logger.info( 'VBS scheme enabled ... ' )
        else :
            logger.info( 'VBS scheme not enabled ... ' )
        #endif
        
        #
        # read tracers
        #

        # name of file with tracer info:
        species_file = rcf.get( basekey+'.tracers.file' )
        
        # info ...
        logger.info( '      read tracers from "%s" ...' % species_file )

        # read:
        specs = genes_tracers.GeneSTracers( species_file, groups, props, skip_tracers )

        #
        # read reactions
        #

        # name of file with reaction info:
        reactions_file = rcf.get( basekey+'.reactions.file' )
        
        # info ...
        logger.info( '      read reactions from "%s" ...' % reactions_file )

        # read reactions:
        reacs = genes_reactions.GeneSReactions( reactions_file, specs, groups, labelled_specs )    

        #
        # write source code
        #
        
        # info ...
        logger.info( '      write source files ...' )

        # write code:
        genes_model_le.WriteSource( rcf, basekey, specs, reacs, with_labeling, labelled_specs, with_vbs, os.curdir )
    
    #enddef   # __init__

#endclass    # GeneS


# ------------------------------------------------
# end
# ------------------------------------------------



