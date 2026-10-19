
ObsSummary README for the OSSOSv11 release 

-------------------------------------------------------------------
Contents

Each subdirectory holds exactly one survey: its pointings.list, the .eff
efficiency file of every block named in pointings.list, and its detections.

CFEPS_g         subdir  The CFEPS blocks acquired in the g filter.
                        This includes the 6 objects detected in CFEPS and
                        (blindly, not targeted) re-observed by OSSOS. They are
                        listed under their CFEPS designations, with their CFEPS
                        detection circumstances. They also appear in OSSOS_v11
                        under their OSSOS designations, with their OSSOS
                        detection circumstances; the correspondence is in the
                        top-level Duplicates file. A combined
                        CFEPS_g + OSSOS_v11 run counts them twice.
CFEPS_r         subdir  The CFEPS L3h block, acquired in MegaPrime r filter.
CFEPS_presurvey subdir  The CFEPS presurvey block (Jones et al. 2006), R filter.
HiLat           subdir  The High Latitude survey, r filter.
MA              subdir  The mal and mah blocks (Alexandersen et al. 2016),
                        r filter.
OSSOS_v11       subdir  All detections from OSSOS and only OSSOS, r filter.
jwst_sample_a   subdir  JWST GO 1568 Sample A (Eduardo et al. 2026); a
                        multi-epoch survey root, not part of the v11 release.
N26             subdir  HST GO-9433 heliostack (Napier et al. 2026), F606W
                        STMAG (filter t); not part of the v11 release. See
                        N26/README.n26.

Detection files for the characterized samples can be found in the above
directories. File names are <sample>/<sample>.detections and, where
available, <sample>/<sample>.CDS. As an example, the OSSOS files are in
OSSOS_v11 as:

OSSOS_v11.detections    Orbital elements and classes. This is a
                        space-separated value file, with name of the column
                        given on the last comment line (starting with #) at the
                        start of the file. See header of CDS file for a
                        description of the columns.
OSSOS_v11.CDS           Same as OSSOS_v11.detections in CDS format

-------------------------------------------------------------------
File Formats

-The .detections file gives an averaged magnitude (and band, and uncertainty)
    for each object DURING THE DISCOVERY NIGHT'S TRIPLE (this is the only
    magnitude that matters for the characterization).  The 'surmised' absolute
    magnitude (that is, the H_r magnitude surmised using the average m_r and
    the discovery geometry) is also given (it has the same error). The dynamical
    classification of the objects is given in the first several columns; see the
    wiki page and file headers for details.

-------------------------------------------------------------------
Caveats/clarifications for the v11 release

- There are 49/37/84/67/147/87/67/55/105/146 characterized objects in
    E/O/L/H/P/M/S/T/C/D block, respectively.  
    The unclassified uo3eXX, uo3oXX, uo3lXX, uo4hXX, uo5pXXX, uo5mXX, uo5sXX,
    uo5tXX, uo5cXXX, uo5dXXX objects are not part of the release but can be
    obtained by request to coreossos@taos.asiaa.sinica.edu.tw  although
    they are listed in the file in the top level directory
- Orbital classifications are 'secured' according to SSBN08 nomenclature. 
  The detections file gives an 'S'=secure or 'I'=insecure rating
- Orbital classifications for less than 10% of the objects is still insecure.
- Note that usage of the CFEPS survey will require the user make choices related
    to what they believe the g-r colour distributions are...
- All products for CFEPS, Hilat and MA (abbrevation for the two survey blocks
    mal and mah done in Mike Alexandersen's thesis) were previously published,
    and are being provided in this release to the OSSOS team for convenience. 
