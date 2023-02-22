LOTOS-EUROS Kalman Filter - Usage
=================================



Quick start
-----------

Choose the patch version; '000' is the first patch of this release.

To start a run with default settings for this patch:

  ./base/000/bin/setup_lekf  base/000/rc/lekf-v1.0.000.rc


Patches
-------

000
  Parallel over modes, supports LE v2.1.001.

001
  Parallel over modes for LE v2.1.004 on single domain.

002
  Parallel over domains using LE v2.1.004 .


Configuration
-------------

A template for the configuration of the Kalman Filter run is:

  base/000/rc/lekf-v3.0.000.rc
  
This file includes the settings for the model.
If the model version is 'v2.0.001', then the included settings are:

  base/000/rc/lekf-lotos-euros-v2.0.001.rc
  

Soure code
----------

The source code consists of the following files:

 1. First the LEKF base source is copied:

        base/000/
        
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

        proj/lotos-euros/000/v2.0.001/

    An example of a change is the modification of the 'le_driver'
    to let deposition velocities be computed by the filter timestep
    initialization, to allow adding of uncertainty to this fields.


