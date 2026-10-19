
--------------------------------------------------------------------------------

This directory contains the pointing history, efficiency functions for each
block, and the list of detected objects for the mal and mah blocks of Mike
Alexandersen's thesis (Alexandersen et al. 2016, AJ), observed in the
MegaPrime r filter.

The content of this directory is:

MA
   \---------- README.ma (this file)    Describes the content of the directory
    \--------- pointings.list           List of pointings and times of the
     \                                  survey
      \------- ma11la.eff               Efficiency function for the mal block
       \------ ma12ha.eff               Efficiency function for the mah block
        \----- MA.detections            List of objects really detected in
                                        mal and mah, with their dynamical class
                                        if known, magnitude at detection, rate
                                        of motion, H mag, and more.

See CFEPS_g/README.formats for the format of a survey "block".

--------------------------------------------------------------------------------
