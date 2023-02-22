#! /usr/bin/env python

# ------------------------------------------------
# help
# ------------------------------------------------

"""
GeneS - Generate Sources

Class to store info on reactions

"""

# ------------------------------------------------
# common
# ------------------------------------------------

# get common stuff:
#   - logger instance
from genes_common import *


# ------------------------------------------------
# class
# ------------------------------------------------

class GeneSReactions( object ) :

    """
    Class to store info on reactions.
    """
    
    def __init__( self, txtfile, specs, groups, labelled_specs =[] ) :

        """
        Read CBm4 reactions from table.

        Example of input file:

          # Label, Reactants, Products         , Rate Expression          , groups 
          R1  , NO2         , NO + O           , 1.0 x <NO2_SAPRC99>      , cbm4
          R2  , O+O2+M      , O3 + M           , 6.0E-34^-2.4             , cbm4

          R.. , TERP + OH   ,                  , 60.*1.2E-11*2.7E10 @ 732 , cbm4 soa  

        Provide specs instance of GeneSTracers class.
        
        Provide list with selected groups, reactions are only included
        if all their groups are in the "groups['selected']" list.
        """

        # external:
        import os
        import numpy

        # file parameters:
        fcomment = '#'
        fcontinue = '...'
        sep = ','

        # found ?
        if not os.path.exists(txtfile) :
            logger.error( 'input file not found:' )
            logger.error( '  %s' % txtfile )
            raise Exception
        #endif

        # dimensions:    
        nspec = specs.nspec

        # supported auxilary parameters:
        aux = {}
        aux['air'    ] = { 'desc' : 'air (concentration)' }
        aux['H2O'    ] = { 'desc' : 'water vapour (concentration)' }
        aux['O2'     ] = { 'desc' : 'oxygen molecule (concentration)' }
        aux['H2'     ] = { 'desc' : 'hydrogen (concentration)' }
        aux['p'      ] = { 'desc' : 'air pressure (Pa)' }
        aux['T'      ] = { 'desc' : 'temperature (K)' }
        aux['zen'    ] = { 'desc' : 'solar zenith angle ?' }
        aux['cloud'  ] = { 'desc' : 'cloud fraction [0-1]' }
        aux['cldfac' ] = { 'desc' : 'cloud fraction applied to photolysis rates [0-1]' }
        aux['rh'     ] = { 'desc' : 'relative humidity at 2m' }
        aux['cldsulf'] = { 'desc' : 'cloud factor for sulf reaction' }
        aux['hetn2o5'] = { 'desc' : 'heterogene chemistry' }
        aux['nox0'   ] = { 'desc' : 'chem_solve_nox_value' }
        aux['rat75'  ] = { 'desc' : 'ratio between reaction rates for lumping reactions' }
        aux['rat10'  ] = { 'desc' : 'ratio between reaction rates for lumping reactions' }
        aux['alpha8' ] = { 'desc' : 'factor used for reverse reactions of ozone photolysis' }
        aux['alpha10'] = { 'desc' : 'factor used for lumping reactions' }
        aux['ppm_to_mlccm3'] = { 'desc' : 'factor used for conversion from ppm to mlc/cm3' }
        aux['ppb_to_mlccm3'] = { 'desc' : 'factor used for conversion from ppb to mlc/cm3' }
        #aux['SO4a_f' ] = { 'desc' : 'sulfate aerosol fine mode' }
        #aux['Na_f'   ] = { 'desc' : 'natrium aerosol fine mode' }
        #aux['Na_c'   ] = { 'desc' : 'natrium aerosol coarse mode' }
        # collect names:
        aux_names = aux.keys()
        # actually used names:
        aux_used = ['air','H2O','p','T','zen','cloud','cldfac','rh','cldsulf',
                    'hetn2o5', 'nox0',
                    'rat75','rat10','alpha8','alpha10', 'ppm_to_mlccm3', 'ppb_to_mlccm3']

        # init loop:
        ireac = 0   # reaction number : 1, 2, ...
        self.label     = []
        self.reactants = []
        self.products  = []
        self.rate      = []
        self.groups    = []
        self.yp  = ['' for i in range(nspec)]
        self.yl  = ['' for i in range(nspec)]
        self.yir = [[] for i in range(nspec)]  # reaction indices: 0, 1, ...
        self.drog_loss     = ['' for i in range(nspec) ]
        self.drog_prod_neg = ['' for i in range(nspec) ]
        self.array_prod    = numpy.zeros( (nspec,nspec), dtype = 'S512' )
        self.array_prod    = [ ['' for i in range(nspec) ] for i in range(nspec) ]

        # open file:
        f = open( txtfile, mode='r' )
        
        # no headers yet ...
        headers = []
        
        # read lines until eof
        while True:

            # read line:
            line = f.readline()

            # zero length ? then eof (empty line has always a cr/lf):
            if len(line) == 0 : break

            # skip comment:
            if line.startswith(fcomment) : continue

            # removes leading and ending white space:
            line = line.strip()
            
            # empty line remained ? then next:
            if len(line) == 0 : continue

            # header line ?
            if len(headers) == 0 :
                # split:
                fields = line.split(sep)
                for field in fields : headers.append(field.strip())
                # check ...
                headers_expected = ['label','reactants','products','rate expression','groups']
                if headers != headers_expected :
                    logger.error( 'wrong headers ...' )
                    logger.error( '  found    : %s' % headers )
                    logger.error( '  expected : %s' % headers_expected )
                    raise Exception
                #endif
                # next:
                continue
            #endif

            # continuation ?
            while line.endswith(fcontinue) :
                # read next line:
                line2 = f.readline()
                # add, remove continuation mark:
                line = line[0:-len(fcontinue)] + line2.strip()
            #endwhile

            # split at comma's :
            try :
                label,reactants,products,rate,grouplist = line.split(sep)
            except :
                logger.error( 'could not split line into (label,reactants,products,rate,groups) :' )
                logger.error( line )
                raise Exception
            #endtry
            
            # convert into list:
            groupnames = []
            for group in grouplist.split() : groupnames.append(group.strip())
            
            # check ...
            for group in groupnames :
                if group not in groups['supported'] :           
                    logger.error( 'reaction group "%s" not one of supported : %s' % (group,str(groups['supported'])) )
                    raise Exception
                #endif
            #endfor

            # all groupnames present ?
            all_present = True
            for group in groupnames :
                all_present = all_present and (group in groups['selected'])
            #endfor
            # all should be present ...
            if not all_present : continue
            
            # accepted; increase counter:
            ireac = ireac + 1

            ## debug ...
            #print( line.strip() )

            # cleanup:
            label = label.strip()

            ## info:
            #logger.info( '' )
            #logger.info( label )
            #logger.info( '  reactants: %s' % reactants )
            #logger.info( '  products : %s' % products )
            #logger.info( '  rate     : %s' % rate )
            
            # store:
            self.label.append(label)
            self.reactants.append(reactants)
            self.products.append(products)
            self.rate.append(rate)
            self.groups.append(groupnames)
            
            # VBS loss reaction rates should not be included in loss/prod terms
            # the equivalent reactions are already in, these rates are only needed (hardcoded) to determine 
            # branch ration for nox regime in the VBS chemistry
            if label.__contains__('VBS_R') : continue

            # 
            # Reactants
            #
            #   A + B + C + ...  ->
            #
            # For each reactant, add the loss rate to the loss function:
            #   yl(A) = .. + k * B * C
            #   yl(B) = .. + k * A * C
            #      :
            #      
            # reactants are seperated by plus signs:
            rcts = reactants.split('+')

            # 
            # Products
            #
            #   A + B + C + ...  ->  p * P + q * Q + r * R ...
            #
            # For each product, add the production rate to the production function:
            #   yp(P) = .. + k * A * B * C
            #   yp(Q) = .. + k * A * B * C
            #      :
            #                
            # for negative productions, change form '- 0.2*PAR' into '+ -0.2 * PAR'
            products = products.replace( '-', '+ -' )
            # products are seperated by plus signs:
            prds = products.split('+') 
            
            ## info ...
            #logger.info( '  create loss terms ...' )

            # loop over reactants to make loss terms
            for rct in rcts :
            
                ## info ...
                #logger.info( '    rct = %s' % rct )
            
                # cleanup:
                spec = rct.strip()
                # no multiple body reactions yet ...
                if '*' in spec :
                    logger.error( 'no factors in source terms supported : %s' % spec )
                    logger.error( 'do not use "2*A -> B" but "A + A -> B" instead' )
                    raise Exception
                #endif
                # search in tracer list:
                if spec not in specs.name :
                    # molucele 'M' is air ...
                    if spec == 'M' : spec = 'air'
                    # should be an auxilary parameter ...
                    if spec not in aux_names :
                        logger.error( 'Reactant is neither a specie or an aux : %s' % spec )
                        raise Exception
                    #endif
                    # no loss expression for an auxilary reactant ...
                    continue
                #endtry
                # index:
                ispec = specs.name.index(spec)
                # add reaction index:
                if ireac-1 not in self.yir[ispec] : self.yir[ispec].append(ireac-1)

                # check if the reactant is also produced by the reaction and thus available 
                # on the product side of the reaction.
                # if this is the case and the production term equals one, then the specie is a Catalysis
                # and it will not be used in production or loss rates
                Catalysis = False
                for prd in prds :
                    # cleanup:
                    prd = prd.strip()
                    # no supported species produced ? then skip:
                    if len(prd) == 0 : continue
                    # exact match ? then leave: Catalysis!
                    if prd.strip() == rct.strip() :
                        Catalysis = True
                        break
                    #endif
                    # split in factor and specie:
                    ind = prd.find('*')
                    if ind > 0 :
                        # if a factor is available before the product, this factor is changed negative
                        #fac,sep,spec = prd.partition('*')
                        ss = prd.split('*') ; 
                        if ss[0].strip()[0] == '-' :  # negative factor changed to positive
                            fac = ss[0].strip()[1:] ;
                            prdspec = ss[1]
                        else :                        # positive factor changed to negative
                            fac = '-' + ss[0]
                            prdspec = ss[1]
                        #endif    

                    else :
                        # If no factor is available (this means factor is 1 or minus 1)

                        if prd[0].strip() == '-' : # negative factor changed to positive
                            fac = '1.0'
                            prdspec = prd[1:]
                        else :                     # positive factor changed to negative
                            fac = '-1.0'
                            prdspec = prd
                        #endif
                    #endif
                    # cleanup:
                    fac = fac.strip()
                    prdspec = prdspec.strip()  
                    
                #endfor

                # Catalysis ?
                if Catalysis : continue

                # init contribution to loss line:
                self.yl[ispec] = self.yl[ispec] + ' &\n        + rk(ireac_%s)' % label
                # loop over other reactants:
                for orct in rcts :
                    # skip self:
                    if orct == rct : continue
                    # get representation:
                    ospecval,aux_used = self.get_spec_repr( orct, specs, aux_names, aux_used )
                    ## info ...
                    #logger.info( '      other reactant : %s' % orct )
                    #logger.info( '        representation : %s' % ospecval )
                    # add multiplication with other concentration:
                    self.yl[ispec] = self.yl[ispec] + ' * %s' % ospecval       
                    
                    # soa precursor ?
                    if ( 'soa_prec' in specs.props[ispec] ) and ( orct.strip().lower() in ['o3', 'no3', 'oh' ] ) :
                        self.drog_loss[ispec] = self.drog_loss[ispec] + '&\n         + rk(ireac_%s)' %(label)
                        self.drog_loss[ispec] = self.drog_loss[ispec] + '* %s' %(ospecval)
                    #endif
                #endfor

                # info ...
                #print( 'add to yl(',specs.name[ispec],') : ', yl[ispec] )
            #endfor 

            ## info ...
            #logger.info( '  create production terms ...' )

            # loop over products, to make production terms
            for prd in prds :            
                # cleanup:
                prd = prd.strip()

                # Catalysis?
                Catalysis= False
                for rct in rcts :
                    # exact match ?
                    if rct.strip() == prd.strip() : 
                        Catalysis = True
                    #endif
                if Catalysis : continue

                ## info ...
                #logger.info( '    prd = %s' % prd )

                # no supported species produced ? then skip:
                if len(prd) == 0 : continue
                # split in factor and specie:
                ind = prd.find('*')
                if ind > 0 :
                    #fac,sep,spec = prd.partition('*')
                    ss = prd.split('*') ; fac = ss[0] ; spec = ss[1]
                else :
                    # might be a lonely '-' sign ...
                    if prd[0] == '-' :
                        fac = '-1.0'
                        spec = prd[1:]
                    else :
                        fac = ''
                        spec = prd
                    #endif
                #endif                     

                # cleanup:
                fac = fac.strip()
                spec = spec.strip()


                # search in tracer list:
                try :
                    ispec = specs.name.index(spec)
                except :
                    # molucele 'M' is air ...
                    if spec == 'M' : spec = 'air'
                    # should be an auxilary conc ...
                    if spec not in aux_names :
                        logger.error( 'Product is neither a specie or an aux : %s', spec )
                        logger.error( '  line : %s' % line )
                        logger.error( '  file : %s' % txtfile )
                        logger.warning( '   continue ...' )
                        #raise Exception
                    #endif
                    # no loss expression for an auxilary reactant ...
                    continue
                #endtry

                # add reaction number:
                if ireac-1 not in self.yir[ispec] : self.yir[ispec].append(ireac-1)            

                # labelled specie?
                if (spec.lower() in labelled_specs) :

                    # state properties of produced specie    
                    Carbon = False
                    Nitrogen = False              
                    Sulphur = False

                    if specs.formula[ispec]['C'] > 0 : Carbon = True
                    if specs.formula[ispec]['N'] > 0 : Nitrogen = True
                    if specs.formula[ispec]['S'] > 0 : Sulphur = True                        

                    # init number of specific sources in reactants
                    ncarbon_rct = 0
                    nnitrogen_rct = 0
                    nsulphur_rct = 0

                    nrspec_labelled = 0 # number of labelled reactants

                    # count number of specific sources in reactants                                                          
                    for rct in rcts :
                        rspec = rct.strip()

                        if rspec.lower() not in labelled_specs :
                            continue
                        else :
                            nrspec_labelled = nrspec_labelled + 1
                        #endif

                        try :
                            irspec = specs.name.index(rspec)
                            if specs.formula[irspec]['C'] > 0 : ncarbon_rct = ncarbon_rct + 1
                            if specs.formula[irspec]['N'] > 0 : nnitrogen_rct = nnitrogen_rct + 1
                            if specs.formula[irspec]['S'] > 0 : nsulphur_rct = nsulphur_rct + 1
                        except :
                            # molucele 'M' is air ...
                            if rspec == 'M' : rspec = 'air'
                            # should be an auxilary conc ...
                            if rspec not in aux_names :
                                logger.error( 'Product is neither a specie or an aux : %s', rspec )
                                logger.error( '  line : %s' % line )
                                logger.error( '  file : %s' % txtfile )
                                logger.warning( '   continue ...' )
                                #raise Exception
                            #endif
                            # no loss expression for an auxilary reactant ...
                            continue
                        #endtry
                    #endfor

                    # if none of the reactants is labelled, add to the 'normal' production term (only possible for negative production of NO2 out of Carbon hydroxation)
                    if nrspec_labelled == 0 :

                        # init line:
                        self.yp[ispec] = self.yp[ispec] + ' &\n       + '

                        # add reaction rate:
                        self.yp[ispec] = self.yp[ispec] + ('rk(ireac_%s)' % label)

                        # negative produced soa precursor ?
                        if len(fac) > 0 :
                            if (fac.strip()[0] == '-' ) and ( 'soa_prec' in specs.props[ispec] ) :
                                self.drog_prod_neg[ispec] = self.drog_prod_neg[ispec] + ' &\n       + '
                                self.drog_prod_neg[ispec] = self.drog_prod_neg[ispec] + ( 'rk(ireac_%s)' %(label) )
                                self.drog_prod_neg[ispec] = self.drog_prod_neg[ispec] + ' * (' + fac + ')'
                            #endif
                        #endif

                        # add production factor ?
                        if len(fac) > 0 : 
                            # enclose if necessary:
                            if fac[0] == '-' : fac = '('+fac+')'
                            # add:
                            self.yp[ispec] = self.yp[ispec] + ' * ' + fac
                        #endif

                        # loop over reactants:
                        for rct in rcts :
                            # cleanup:
                            rspec = rct.strip()
                            # get representation:
                            rspecval,aux_used = self.get_spec_repr( rspec, specs, aux_names, aux_used )
                            # add multiplication with other concentration:
                            self.yp[ispec] = self.yp[ispec] + ' * %s' % rspecval

                            # negative produces soa precursor ?
                            if len(fac) > 0 :
                                if (fac[0] == '-' ) and ( 'soa_prec' in specs.props[ispec] ) :
                                    self.drog_prod_neg[ispec] = self.drog_prod_neg[ispec] + '* %s' %(rspecval)
                                #endif
                            #endif
                        #endfor

                    else :
                        rspecval_label = ''                        

                        if (Carbon and ncarbon_rct > 1) : rspecval_label = rspecval_label + ' * %1.2f' %(float(1.0/ncarbon_rct))
                        if (Nitrogen and nnitrogen_rct > 1) : rspecval_label = rspecval_label + ' * %1.2f' %(float(1.0/nnitrogen_rct))                                
                        if (Sulphur and nsulphur_rct > 1) : rspecval_label = rspecval_label + ' * %1.2f' %(float(1.0/nsulphur_rct))   
                        if (Carbon and Nitrogen and ncarbon_rct >= 1 and nnitrogen_rct >= 1) : rspecval_label = ' * 0.50 ' + rspecval_label                             

                        if (len(fac) > 0 ) :      
                            # enclose if necessary:
                            if fac[0] == '-' : fac = '('+fac+')'
                            # add:
                            rspecval_label = rspecval_label + ' * ' + fac                                                        
                        #endif

                        # loop over reactants
                        for rct in rcts :

                            rspec = rct.strip()
                            rspecval,aux_used = self.get_spec_repr( rspec, specs, aux_names, aux_used )
                            rspecval_label = rspecval_label + ' * ' + rspecval                             

                        #endfor

                        # loop over reactants
                        for rct in rcts :
                            rspec = rct.strip()
                            try:
                                irspec = specs.name.index(rspec)

                            except :
                                # molucele 'M' is air ...
                                if rspec == 'M' : rspec = 'air'
                                # should be an auxilary conc ...
                                if rspec not in aux_names :
                                    logger.error( 'Product is neither a specie or an aux : %s', rspec )
                                    logger.error( '  line : %s' % line )
                                    logger.error( '  file : %s' % txtfile )
                                    logger.warning( '   continue ...' )
                                    #raise Exception
                                #endif
                                # no loss expression for an auxilary reactant ...
                                continue
                            #endtry
                            if (Carbon and not Nitrogen and specs.formula[irspec]['C'] > 0 and rspec.lower() not in labelled_specs) :
                                logger.error( prd + ' is labelled, while origin ' + rspec + ' is not labelled')
                                raise Exception
                            elif (Carbon and not Nitrogen and specs.formula[irspec]['C'] > 0 ):
                                self.array_prod[irspec][ispec] = self.array_prod[irspec][ispec] + (' &\n        + rk(ireac_%s)%s ' %(label,rspecval_label) )

                            elif (Nitrogen and not Carbon and specs.formula[irspec]['N'] > 0 and rspec.lower() not in labelled_specs):
                                logger.error( prd + ' is labelled, while origin ' + rspec + ' is not labelled')
                                raise Exception
                            elif (Nitrogen and not Carbon and specs.formula[irspec]['N'] > 0 ) :
                                self.array_prod[irspec][ispec] = self.array_prod[irspec][ispec] + (' &\n        + rk(ireac_%s)%s ' %(label,rspecval_label) )

                            elif (Sulphur and specs.formula[irspec]['S'] > 0 and rspec.lower() not in labelled_specs):
                                logger.error( prd + ' is labelled, while origin ' + rspec + ' is not labelled')
                                raise Exception
                            elif (Sulphur and specs.formula[irspec]['S'] > 0 ):
                                self.array_prod[irspec][ispec] = self.array_prod[irspec][ispec] + (' &\n        + rk(ireac_%s)%s ' %(label,rspecval_label) )

                            elif ( (Carbon and Nitrogen) and ( (specs.formula[irspec]['C'] > 0 and rspec.lower() in labelled_specs) or (specs.formula[irspec]['N'] > 0  and rspec.lower() in labelled_specs) ) ) :
                                self.array_prod[irspec][ispec] = self.array_prod[irspec][ispec] + (' &\n        + rk(ireac_%s)%s ' %(label,rspecval_label) )

                            else :
                                continue
                            #endif

                        #endfor # reactants
                    #endif # labelled reactants?


                else :


                    # init line:
                    self.yp[ispec] = self.yp[ispec] + ' &\n        + '

                    # add reaction rate:
                    self.yp[ispec] = self.yp[ispec] + ('rk(ireac_%s)' % label)

                    # negative produced soa precursor ?
                    if len(fac) > 0 :
                        if (fac.strip()[0] == '-' ) and ( 'soa_prec' in specs.props[ispec] ) :
                            self.drog_prod_neg[ispec] = self.drog_prod_neg[ispec] + ' &\n       + '
                            self.drog_prod_neg[ispec] = self.drog_prod_neg[ispec] + ( 'rk(ireac_%s)' %(label) )
                            self.drog_prod_neg[ispec] = self.drog_prod_neg[ispec] + ' * (' + fac + ')'
                        #endif
                    #endif

                    # add production factor ?
                    if len(fac) > 0 :
                        # enclose if necessary:
                        if fac[0] == '-' : fac = '('+fac+')'
                        # add:
                        self.yp[ispec] = self.yp[ispec] + ' * ' + fac
                    #endif

                    # loop over reactants:
                    for rct in rcts :
                        # cleanup:
                        rspec = rct.strip()
                        # get representation:
                        rspecval,aux_used = self.get_spec_repr( rspec, specs, aux_names, aux_used )
                        # add multiplication with other concentration:
                        self.yp[ispec] = self.yp[ispec] + ' * %s' % rspecval

                        # negative produces soa precursor ?
                        if len(fac) > 0 :
                            if (fac[0] == '-' ) and ( 'soa_prec' in specs.props[ispec] ) :
                                self.drog_prod_neg[ispec] = self.drog_prod_neg[ispec] + '* %s' %(rspecval)
                            #endif
                        #endif
                    #endfor
                #endif

                # info ...
                #print( 'add to yp(',specs.name[ispec],') : ', yp[ispec] )
            #endfor

            # testing ...
            #if ireac == 2 : break
            #if label.strip() == 'RH2f' : raise Exception

        #endwhile

        # close:
        f.close()
        
        # info:
        #print( '' )

        # set counters:
        self.nreac = ireac

        # lists with actually used auxilary parameters:
        self.aux_name = [] ; self.aux_desc = []
        for name in aux_used :
            self.aux_name.append( name)
            self.aux_desc.append( aux[name]['desc'] )
        #endfor
        # counter:
        self.naux = len(self.aux_name)

    #enddef  # __init__
    
    # *
    
    def get_spec_repr( self, orct, specs, aux_names, aux_used ) :
    
        """
        Get represenation for reactant.
        """
    
        # cleanup:
        ospec = orct.strip()
        # multiple body reaction ?
        if '*' in ospec :
            # split:
            ofac,ospec = ospec.split('*')
            # cleanup:
            ofac = ofac.strip()
            ospec = ospec.strip()
        else :
            ofac = ''
        #endif
        # search in tracer list:
        if ospec in specs.name :
            # representation of the value:
            ospecval = 'y(ispec_%s)' % ospec
        else :
            # molucele 'M' is air ...
            if ospec == 'M' : ospec = 'air'
            # auxilary conc ?
            if ospec in aux_names :
                # representation of the value:
                ospecval = 'aux(iaux_%s)' % ospec
                # add to list:
                if ospec not in aux_used : aux_used.append(ospec)
            else :
                print( 'reactant is neither a specie or an aux : ', ospec )
                raise
            #endif
        #endif
        # add factor if necessary:
        if len(ofac) > 0 : ospecval = ofac+'*'+ospecval
        
        # ok
        return ospecval,aux_used
        
    #enddef   # get_spec_repr


#endclass   # GeneSReactions


# ------------------------------------------------
# end
# ------------------------------------------------

