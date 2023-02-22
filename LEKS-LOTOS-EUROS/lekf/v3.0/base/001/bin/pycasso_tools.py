

"""
PYCASSO tools
"""


#-------------------------------------------------
# routines
#-------------------------------------------------


def write_text_file( fname, text ) :

    """
    Write a text file.
    Arguments:
      fname   : target file name
      text    : list of strings, a line should end with '\n'
    """
    
    # external:
    import os
    import logging
    
    # info ...
    logging.debug( '    write %s ...' % fname )

    # write new text:
    f = open( fname, 'w' )
    f.writelines( text )
    f.close()

    # ok
    return
    
#enddef


# ***


def diff_text_files( fname1, fname2 ) :

    """
    Return True if the files are different, False if they are the same.
    """
    
    # external:
    import os
    import logging
    
    # read files:
    f = open( fname1, 'r' )
    lines1 = f.readlines()
    f.close()
    f = open( fname2, 'r' )
    lines2 = f.readlines()
    f.close()
    
    # compare:
    return lines1 != lines2
    
#enddef


def replace_text_file_if_different( targetfile, newfile ) :

    """
    Replace target file with new file if the content is different.
    """
    
    # external:
    import os
    
    # target file not present yet ?
    if not os.path.isfile(targetfile) :
        # rename:
        os.rename( newfile, targetfile )
    else :
        # check if different ...
        different = diff_text_files( targetfile, newfile )
        # different ?
        if different :
            # different files, replace target file by new generated version:
            os.rename( newfile, targetfile )
        else :
            # files are the same, remove the temporary version:
            os.remove( newfile )
        #endif
    #endif            

#enddef


# ***


def update_text_file( fname, newtext ) :

    """
    Replace a file by a new text if the later differs
    from the current content.
    Arguments:
      fname      : target file name
      newtext    : list of strings, a line should end with '\n'
    """
    
    # external:
    import os
    import logging
    
    ## info ...
    #logging.info( '  update %s ...' % fname )
    # file exists alread?
    if os.path.exists(fname) :
        # read current content:
        f = open( fname, 'r' )
        oldtext = f.readlines()
        f.close()
        # differences ?
        rewrite = newtext != oldtext
        ## info ...
        #for iline in range(len(oldtext)) :
        #    if iline < len(newtext) :
        #        if oldtext[iline] != newtext[iline] :
        #            logging.info( '    first different lines:' )
        #            logging.info( '      old: %s' % oldtext[iline] )
        #            logging.info( '      new: %s' % newtext[iline] )
        #            break
        #        #endif
        #    #endif
        ##endfor
    else :
        # no file yet, always rewrite:
        rewrite = True
    #endif
    
    # write file ?
    if rewrite :
        # for info message:
        stat = 'replace'
        # write new text:
        f = open( fname, 'w' )
        f.writelines( newtext )
        f.close()
    else :
        # for info message:
        stat = 'keep'
    #endif
    
    # info ...
    logging.debug( '    %-8s %-40s' % (stat,fname) )

    # ok
    return
    
#enddef


# ***


def modify_text_file( fname, key, value ) :

    """
    Modify a text file by replacing a key by a new value.
    Arguments:
      fname      : text file name
      key        : value to be replaced
      value      : replacement value
    """

    # modules:
    import os
    import logging

    # check ...
    if not os.path.isfile(fname) :
        logging.error( 'no file or file not found : %s' % fname )
        raise IOError
    #endif
    
    # read file into list of lines:
    f = open( fname, 'r' )
    lines = f.readlines()
    f.close()

    # copy while replacing key with value:
    newlines = []
    for line in lines :
        # replace key with value:
        line = line.replace( key, value )
        # add:
        newlines.append(line)
    #endfor
    
    # write again:
    f = open( fname, 'w' )
    f.writelines(newlines)
    f.close()

    # ok
    return
    
#enddef


# ***


def CreateDirs( dirname, forceclean=False ) :

    """
    Create a directory and report success.
    """

    # external:
    import os
    import shutil
    import logging
    
    # already present ?
    if os.path.isdir( dirname ) :
        # remove existing directory ?
        if forceclean:
            # info ...
            logging.info( 'remove existing %s ...' % dirname )
            # use the chainsaw:
            shutil.rmtree( dirname )
        #endif
    #endif
    
    # not present (anymore) ?
    if not os.path.isdir( dirname ) :
        # info ...
        logging.info( 'create new directory %s ...' % dirname )
        # create directory including parents if necesary:
        os.makedirs(dirname)
    #endif
    
    # ok
    return None

#enddef

            
#-------------------------------------------------
# end
#-------------------------------------------------
