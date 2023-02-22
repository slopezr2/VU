LEKF - LOTOS-EUROS Kalman Filter
================================

Ensemble Kalman filter around LOTOS-EUROS model.


Quick clone
-----------

The easiest way is to use an ssh url to clone the repository.
See [How to use the repository](https://ci.tno.nl/gitlab/lotos-euros/lekf/-/wikis/How-to-use-the-repository) 
on the wiki for how to add a new ssh private key to your account:

Then use the following command to clone the repository:

    git  clone  git@ci.tno.nl:lotos-euros/lekf.git


Documentation
-------------

A browsable documentation is included in the source tree:

    doc/build/html/index.html

Eventually first re-create the documentation from the its source files:

    make docu

Note that the documenation requies 'sphinx', 
which is usually part of your python distribution.




Quick start
-----------

Choose the patch version; '000' is the first patch of this release.

To start a run with default settings for this patch:

  ./base/004/bin/setup-lekf  base/004/rc/lekf.rc


Patches
-------

000
  Parallel over modes, supports LE v2.1.001.

001
  Parallel over modes for LE v2.1.004 on single domain.

002
  Parallel over domains using LE v2.1.004 .
  
003
  Changes for CAMS50 operational version.
  
004  
  Changes following coupling to LOTOS-EUROS v2.2.001.


Configuration
-------------

A template for the configuration of the Kalman Filter run is:

  base/004/rc/lekf.rc
  
This file includes the settings for the model.
If the model version is 'v2.2.001', then the included settings are:

  base/004/rc/lotos-euros-v2.2.001.rc
  

Soure code
----------

The source code consists of the following files:

 1. First the LEKF base source is copied:

        base/004/
        
    Note that this base is valid for a particular version
    of the model only; see '4' for source files valid for
    other versions.
       
 2. A complete LOTOS-EUROS code, which consists of a base version plus
    eventually some replacements for specific projects.
    The base- and project-directories are specified in the 
    lotos-euros.rc file included in the lekf.rc .

 3. Changes specific to a particular combination of LEKF 
    and LOTOS-EUROS versions.
    These are stored in a path including the patch number of the LEKF code
    and the version number of the model:

        proj/lotos-euros/004/v2.2.001/

    An example of a change is the modification of the 'le_driver'
    to let deposition velocities be computed by the filter timestep
    initialization, to allow adding of uncertainty to this fields.


