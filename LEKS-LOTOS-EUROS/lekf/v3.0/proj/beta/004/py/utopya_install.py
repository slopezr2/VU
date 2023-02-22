
#-------------------------------------------------
# help
#-------------------------------------------------

"""
*************************
``utopya_install`` module
*************************

The :py:mod:`utopya_install` module provides classes to copy files to the local system,
eventually from remote archives. This is typically used to install input data.
All classes are accesible through the top level :py:mod:`utopya` module.

* The :py:class:`.UtopyaInstall` class creates an instance of a top level
  object that, based on rcfile settings, performs all installation tasks
  requested in the settings.
  
The configuration of the objects is done using rcfile settings.


Example script
==============

An example script and settings, as well as test data is available
in the 'test' directory:

* `test/UtopyaInstaller <../../../test/UtopyaInstall>`_
* `test/UtopyaInstaller.rc <../../../test/UtopyaInstall.rc>`_


Class hierarchy
===============

The classes provided by this module have been derived with the following hierchy:

* :py:class:`.UtopyaRc`

  * :py:class:`.UtopyaInstall`

Classes
=======

"""


#-------------------------------------------------
# modules
#-------------------------------------------------

import utopya_rc


#-------------------------------------------------
# builder
#-------------------------------------------------

class UtopyaInstall( utopya_rc.UtopyaRc ) :

    """
    Base class for objects that install input files.
    
    The initialization arguments look like::
    
        UtopyaInstall( 'settings.rc', rcbase='model.install', 
                          env={ 'WORKDIR' : '/scratch/yourname/work' } )
                          
    The first argument defines the name of the rcfile with settings,
    and the second is the base name of the settings that will be used.
    The optional 'env' dictionairy could be used to pass environment variables
    that are expanded in the rcfile.
    
    The first configuration in the settings should be the definition of a 
    set list. This is just a list of keywords with nicknames for
    data sets to be installed::
    
        ! list of data sets to be installed:
        model.install.sets       :  meteo obs
        
    For each of these sets, the details should be specified with keys
    that start with ``<rcbase>.<setname>``. Here we show an example for
    the set named 'meteo', which should install a timeseries of meteo files.
    
    First specify the target directory of the installation::
    
        ! installation directory; empty for current:
        model.install.meteo.dir        : /scratch/run/meteo
        
    Specify the archive directory where the input files are located.
    This might be a remote path,
    in which case the specification should follow the formatting
    supported by the :py:mod:`gss` module::
    
        ! archive path with files to be installed:
        model.install.meteo.arch     :   ec:/you/data/meteo
        
    Eventually specify a list of archive directories;
    these will be scanned in the specified order until the source
    file is found::
    
        ! archive paths with files to be installed:
        model.install.meteo.arch     :   ec:/you/data/meteo \\
                                         /data/other-meteo
                                         
    Specifiy a list of files to be copied from the archive;
    an error will be raised if a file is not found::
    
        ! files to be installed:
        model.install.meteo.files    :   t.nc q.nc u.nc v.nc
        
    If a file has the extension that is recoqnized as a package
    ('.tar','.tar.gz','.tgz', or '.zip'), 
    then the file is automatically unpacked and removed::
    
        ! package files, unpacked automatically:
        model.install.meteo.files    :   t.tar q.tar
        
    Files that only need to be installed if present in the archive
    should be enclosed by square brackets::
    
        ! package files, unpacked automatically:
        model.install.meteo.files    :   t.tar q.tar [extra.nc]
        
    To install a timeseries, first specify a time range and time step size,
    and include templates for year, month, day etc in the filenames.
    The templates should be accepted by the :py:meth:`datetime.datetime.strftime` method::
    
        ! copy files per day:
        model.install.meteo.timerange  :   2000-01-01 to 2010-12-31 by 1 day
        ! filenames:
        model.install.meteo.files      :   t_%Y%m%d.nc q_%Y%m%d.nc

    The general form of the time range description is::
        
        <start-time> [up]to <end-time> by <time-step>
            
    where most common forms for times and steps are supported,
    for example::
    
        2000-01-01 00:00 to 2006-01-01 00:00 by 1.5 hours
        2000-01-01       to 2006-01-01       by 1   days
        2000-01-01       to 2006-01-01       by 1   months
        
    With seperator "to" no file will be enstalled for the end time,
    for example because this is included in the last file already.
    If a file should be retrieved for the end time, use "upto".

    Set the following flag if files should always be renewed;
    by default, files that already exists are kept::
    
        ! renew existing files?
        model.install.meteo.renew      :   True
        
    Another optional flag could be set to reset the modification time
    of each installed file; this is useful on scratch file systems
    that automatically remove files that have not been changed for a while::
    
        ! touch installed or already present files?
        model.install.meteo.touch      :   True

    Enable the following optional flag to have messages printed while
    getting files; this is for example useful to debug the retrieval;
    default is False::

        ! print messages while getting files?
        model.install.meteo.verbose    :   True
        
    Note that package files (.tar, .zip) are always renewed even if the
    'renew' flag is 'False' (default), since they are removed after unpacking.
    To only renew a package file if one of the included files is not present yet,
    use the 'from' construction to specify that a file should be extracted from
    the specified package; note that now a sepeation character ';' should be used::
    
        ! daily files, extract from monthly package; use ';' as seperation:
        model.install.meteo.files      :   t_%Y%m%d.nc from t_%Y%m.tar ; \\
                                           q_%Y%m%d.nc from q_%Y%m.tar
        
    With this, the archive file is only retrieved and (completely) unpacked
    if one of the included files is not present yet.
    """
    
    def __init__( self, rcfile, rcbase='', env={} ) :
    
        """
        Initialize archiver object.
        """
        
        # modules:
        import os
        
        # init base object:
        utopya_rc.UtopyaRc.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # info ...
        self.logger.info( 'install ...' )
        self.logger.info( '  settings from : %s' % rcfile )
        self.logger.info( '  base key      : %s' % rcbase )
        if len(env) > 0 :
            self.logger.info( '  environment:' )
            for key in env.keys() :
                self.logger.info( '    %s = %s' % (key,env[key]) )
            #endfor
        #endif
        
        # read set list:
        setnames = self.GetSetting( 'sets' ).split()
        # info ...
        self.logger.info( '  loop over datasets ...' )
        # loop:
        for setname in setnames :

            # info ...
            self.logger.info( '    data set "%s" ...' % setname )
            
            # call installer:
            self.InstallSet( setname, indent='      ' )

        #endfor # set
        
    #enddef __init__
    
    # *
    
    def InstallSet( self, setname, indent='' ) :
      
        """
        Install dataset following settings specified by name.
        """
        
        # modules:
        import os

        # target directory to install data:
        tdir = self.GetSetting( setname+'.dir' )
        # change base dir?
        if len(tdir) > 0 :
            # not present yet?
            if not os.path.isdir(tdir) :
                # info ...
                self.logger.info( indent+'create %s ...' % tdir )
                # create:
                os.makedirs( tdir )
            #endif
            # info ...
            self.logger.info( indent+'change to %s ...' % tdir )
            # save current:
            owd = os.getcwd()
            # change:
            os.chdir( tdir)
        #endif # basedir defined?
        
        # archive directory:
        archdirs = self.GetSetting( setname+'.arch' ).split()
        # count:
        narch = len(archdirs)
        # info ...
        self.logger.info( indent+'archive location(s):' )
        for archdir in archdirs :
            self.logger.info( indent+'  '+archdir )
        #endfor
        
        # display progress when getting files?
        verbose = self.GetSetting( setname+'.verbose', totype='bool', default=False )
        # renew existing files?
        renew = self.GetSetting( setname+'.renew', totype='bool', default=False )
        # touch all extracted or present files?
        touch = self.GetSetting( setname+'.touch', totype='bool', default=False )
        
        # files to be installed, should be present:
        line = self.GetSetting( setname+'.files', default='' )
        fnames = []
        if len(line) > 0 :
            sep = None
            if (';' in line) or (' from ' in line) : sep = ';'
            for f in line.split(sep) : fnames.append( f.strip() )
        #endif
            
        # time range, might not be defined:
        trange = self.GetSetting( setname+'.timerange', default='None' )
        # defined?
        if trange != 'None' :
            # split timerange definition in start time, end time, and step;
            # return list of time values:
            tt = self.GetTimes( trange )
        else :
            # dummy ...
            tt = [None]
        #endif

        # time loop:
        for t in tt :
            # info ..
            if t is None :
                self.logger.info( indent+'install ...' )
            else :
                self.logger.info( indent+'install for time %s ...' % t.strftime('%Y-%m-%d %H:%M') )
            #endif

            # loop over all files to retrieve:
            for fname_template in fnames :

                # optional?
                if fname_template.startswith('[') and fname_template.endswith(']') :
                    # set flag:
                    optional = True
                    # remove enclosing '[..]':
                    fname_template = fname_template[1:-1]
                else :
                    # should be present:
                    optional = False
                #endif
                # for log messages ...
                omsg = ''
                if optional : omsg = ' (if present)'
                
                # replace time templates?
                if t is None :
                    fname = fname_template
                else :
                    fname = self.ExpandTime( fname_template, t )
                #endif
                
                # source file; either same, or package:
                if ' from ' in fname :
                    # extract source from description:
                    fname,source = fname.split(' from ')
                    # cleanup:
                    fname = fname.strip()
                    source = source.strip()
                else :
                    # no specific source file:
                    source = fname
                #endif
                
                # renew?
                if renew or (not os.path.isfile(fname)) :
                    # info ..
                    if os.path.isfile(fname) :
                        self.logger.info( indent+'  renew %s%s ...' % (fname,omsg) )
                    else :
                        self.logger.info( indent+'  get   %s%s ...' % (fname,omsg) )
                    #endif
                    if source != fname :
                        self.logger.info( indent+'    source %s%s ...' % (source,omsg) )
                    #endif
                    # loop over possible archives:
                    narch = len(archdirs)
                    for iarch in range(narch) :
                        # current:
                        archdir = archdirs[iarch]
                        # file is optional if more than one archive is specified,
                        # and this is not the last attempt yet:
                        trynext = (narch > 1) and (iarch < len(archdirs)-1)
                        # get file, eventually unpack:
                        self.GetFile( archdir, source, optional=(optional or trynext), 
                                        touch=touch, verbose=verbose, indent=indent+'    ' )
                        # found?
                        if os.path.isfile(fname) : break
                    #endfor
                    # file from package ? then check if present ..
                    if (fname != source) and (not os.path.isfile(fname)) and (not optional) :
                        self.logger.error( 'file "%s" not found after installation of "%s"' % (fname,source) )
                        raise Exception
                    #endif
                else :
                    # info ...
                    self.logger.info( indent+'  keep  %s ...' % fname )
                #endif
                
                # touch file?
                if touch and os.path.isfile(fname) : os.utime( fname, None )

            #endfor # files

        #endfor # times

        # change base dir?
        if len(tdir) > 0 :
            # back:
            os.chdir( owd )
        #endif # basedir defined?
        
    #enddef InstallSet
    
    # *
    
    def ExpandTime( self, filename_template, t ) :
    
        """
        Expand time values in filename template:
           input_%Y%m%d.nc
           
        Special for templates with two different times:
           %Y%m%d/endtime_%(%Y%m%d + 1 day).nc
        """
        
        # modules:
        import datetime

        # special template present?
        while '%(' in filename_template :
            # start:
            i0 = filename_template.index('%(')
            # end:
            i1 = i0 + filename_template[i0:].index(')')+1
            # remove markup:
            line = filename_template[i0+2:i1-1]
            # split:
            fields = line.split()
            # extract first elements:
            tfmt = fields[0]
            oper = fields[1]
            # switch:
            if oper == '+' :
                # increment fields:
                n    = int(fields[2])
                step = fields[3]
                # time step:
                if step == 'day' :
                    dt = datetime.timedelta(n)
                else :
                    self.logger.error( 'unsupported step size "%s"' % step )
                    raise Exception
                #endif
                # alternative time:
                tx = t + dt
                # expand:
                value = tx.strftime( tfmt )
                # insert:
                filename_template = filename_template[:i0]+value+filename_template[i1:]
            else :
                # not yet ...
                self.logger.error( 'unsupported operator "%s" in "%s" in template "%s"' % (oper,line,filename_template) )
                raise Exception
            #endif
        #endif

        # expand:
        return t.strftime( filename_template )

    #enddef ExpandTime

    # *
    
    def GetFile( self, path, filename, optional=False, touch=False, verbose=False, indent='' ) :
      
        """
        Get specified file from path.
        
        The path could denote a remote locations, in which case it should
        be supported by the :py:mod:`gss` module.
        
        If the file is a package (tar or zipfile), then it will be unpacked
        by the :py:meth:`.Unpack` method and removed.
        """
        
        # modules:
        import os
        
        # tools:
        import gss
        
        # full name of sourcefile:
        sourcefile = os.path.join(path,filename)
        
        # check presence:
        exists = gss.IsFile( sourcefile )
        # not present?
        if not exists :
            # might be ok ...
            if optional :
                # info ..
                self.logger.warning( indent+'file not found: %s' % sourcefile )
                # ok
                return
            else :
                self.logger.error( indent+'file not found: %s' % sourcefile )
                raise Exception
            #endif
        #endif
        
        # get file:
        gss.Copy( sourcefile, filename, verbose=verbose )
        
        # unpack if necssary, based on filename extension ;
        # limitted to target file if specified:
        self.Unpack( filename, touch=touch, indent=indent )
        
    #enddef GetFile
    
    # *
    
    def Unpack( self, filename, touch=False, indent='' ) :
      
        """
        Unpack file if it has the exension for a tar or zip file.
        After unpacking, the package file is removed.
        
        Unpacked files are uncompressed if they have the extension '.gz'.
        """
        
        # modules:
        import os
        
        # no unpacking by default:
        command = None
        # extract extension:
        bname,ext = os.path.splitext(filename)
        # unpack?
        if ext in ['.tar','.tar.gz','.tgz'] :
            # unpack with tar:
            command = 'tar x -v -f %s' % filename
        #
        elif ext in ['.zip'] :
            # unpack with zip:
            command = 'unzip %s' % filename
        #endif
        # run?
        if command is not None :
            # info ...
            self.logger.info( indent+'unpack ...' )
            # call subprocess, trap errors:
            p = self.Call( command, shell=True, indent=indent+'  ' )
            # loop over extracted files:
            for outfile in p.stdout :
                # extracted file, might be replaced by unpacked version:
                xoutfile = outfile
                # unzip if necessary:
                if outfile.endswith('.gz') :
                    # info ...
                    self.logger.info( indent+'  gunzip %s ...' % outfile )
                    # unzip command, force overwrite if already present:
                    command = 'gunzip -f %s' % outfile
                    # run:
                    self.Call( command, shell=True, indent=indent+'  ' )
                    # reset name:
                    xoutfile = outfile.replace('.gz','')
                #endif
                # reset modification time?
                if touch : os.utime( xoutfile, None )
            #endfor # extracted files
            # info ...
            self.logger.info( indent+'remove ...' )
            # remove:
            os.remove( filename )
        #endif
    
    #enddef Unpack
    
    # *
    
    def GetTimes( self, trange ) :
      
        """
        Split time range description into start time, end time, and step.        
        Return value is a list of :py:class:`datetime.datetime` objects.
        """
        
        # modules:
        import datetime
        import dateutil.relativedelta
        
        # check ..
        if (' to ' not in trange) and (' upto ' not in trange) :
            self.logger.error( 'timerange description should be "<starttime> [up]to <endtime> by <step>", found "%s"' % trange )
            raise Exception
        #endif
        # split:
        if ' upto ' in trange :
            include2 = True
            stime1,trange2 = trange.split(' upto ')
        else :
            include2 = False
            stime1,trange2 = trange.split(' to ')
        #endif
        # check ..
        if ' by ' not in trange2 :
            self.logger.error( 'timerange description should be "<starttime> to <endtime> by <step>", found "%s"' % trange2 )
            raise Exception
        #endif
        # split:
        stime2,ststep = trange2.split(' by ')
        
        # remove seperation characters:
        for c in ['/','-',':'] :
            stime1 = stime1.replace(c,' ')
            stime2 = stime2.replace(c,' ')
        #endfor
        # extract start time:
        try :
            t1 = datetime.datetime( *map(int,stime1.split()) )
        except :
            self.logger.error( 'could not extract starttime from "%s" from description "%s"' % (stime1,trange) )
            raise
        #endtry
        
        # extract end time:
        try :
            t2 = datetime.datetime( *map(int,stime2.split()) )
        except :
            self.logger.error( 'could not extract endtime from "%s" from description "%s"' % (stime2,trange) )
            raise
        #endtry
        
        # split step:
        parts = ststep.split()
        # check ..
        if len(parts) != 2 :
            self.logger.error( 'format of timestep should be "<number> <units>", found "%s"' % ststep )
            raise
        #endif
        # copy:
        sn,units = parts
        
        # extract number:
        try :
            if '.' in sn :
                n = float(sn)
            else :
                n = int(sn)
            #endif
        except :
            self.logger.error( 'could not extract step number from "%s" from description "%s"' % (sn,trange) )
            raise Exception
        #endtry
        
        # translate units "year" to "years" etc for use in 'relativedelta',
        # otherwise these are interpreted as absolute numbers ...
        if units in ['year','month','day','hour','minute','second'] :
            units = units+'s'
        #endif
        # check ..
        if units not in ['years','months','days','hours','minutes','seconds'] :
            self.logger.error( 'unsupported step units "%s" from description "%s"' % (units,trange) )
            raise Exception
        #endif
        
        # time step:
        dt = dateutil.relativedelta.relativedelta( **{ units : n } )
        
        # init result:
        tt = []
        # fill:
        t = t1
        while t < t2 :
            # add:
            tt.append( t )
            # next:
            t = t + dt
        #endwhile
        # add end?
        if include2 : tt.append( t2 )
        
        # ok
        return tt
      
    #enddef GetTimes


#endclass UtopyaInstall


#-------------------------------------------------
# end
#-------------------------------------------------

