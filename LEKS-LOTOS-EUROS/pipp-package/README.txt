
LOTOS-EUROS code Santiago
----------------

Open source version v2.2.001, but also including:
    - M7 aerosol scheme
 - radiance computations

Some specific fixes that will be part of v2.2.002 are copied to:
  proj/pipp/001/


Guides
------

LOTOS-EUROS_User-Guide_v2.2.000.pdf
  Latest version of User Guide (how to run the model).
  
LOTOS-EUROS-Reference-Guide-v2.2.000.pdf
  Latest version of Reference Guide (what is in the model).
  
LETR-2016-001-aerosol-modes3.pdf
  Research notes on aerosol sizes and optical properties.


  
About aerosol radius
--------------------

The aerosol sizes are defined in:
  base/001/data/aerosol-radius.csv

The name of this file is used in the settings:
  base/001/rc/lotos-euros-expert.rc
     genes.radius.file             :  ../data/aerosol-radius.csv
This is then used to create the include file with tracer properties;
check the result in the build directory:
  build/src/le_indices.inc


Work by Jianbing Jin
--------------------

The "aerosol-radius.csv" file was created by Jianbing when 
testing optical properties of dust plumes.
Related to "AEROSOURCE" project.

Draft paper is on Overleaf.
To be shared by Jianbing?

Some usefule stuff might be on TNO project folders in:
  ${PROJECTS}/SMO/Space/2019/AEROSOURCE/users/jinj/
Content:
  Angstrom_Assimilation_paper/
  LE_M2_M7_comparison/
  LOTOS-EUROS/
  v2.2_M7/
  v2.2_Radius/
  
