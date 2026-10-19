c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c Generic model routines for Survey Simulator, version 2.0 for OSSOS
c
c Calling sequence in SurveySimulator.f (survey simulator driver) is:
c
c Loop (some condition):
c     call GiMeObj(arg_list_1)
c     Check model ended:
c         set exit condition
c     call Detos1(arg_list_2)
c     Check detection and tracking:
c         store results
c
c where arg_list_1 is
c (filena, seed, a, e, inc, node, peri, M, epoch, h, color, gb, ph,
c  period, amp, comp, ierr)
c with:
c
c INPUT
c     filena: name of containing description of model, read in by the
c             model subroutine "modname" (CH)
c     seed  : Random number generator seed (I4)
c
c OUTPUT
c     a     : semimajor axis (R8)
c     e     : eccentricity (R8)
c     inc   : Inclination [rad] (R8)
c     node  : Longitude of node [rad] (R8)
c     peri  : Argument of perihelion [rad] (R8)
c     M     : Mean anomaly [rad] (R8)
c     epoch : epoch of the orbital elements, in Julian Day (R8)
c     h     : absolute magnitude of object in band filter "x" (R8)
c     color : array of colors "y-x", where the index of "y" is as
c             described in detos1 (10*R8)
c                color(1) : g-x
c                color(2) : r-x
c                color(3) : i-x
c                color(4) : z-x
c                color(5) : u-x
c                color(6) : V-x
c                color(7) : B-x
c                color(8) : R-x
c                color(9) : I-x
c     gb    : opposition surge factor, Bowell formalism (R8)
c     ph    : phase of lightcurve at epoch [rad] (R8)
c     period: period of lightcurve [day] (R8)
c     amp   : amplitude of lightcurve [mag] (R8)
c     commen: user specified string containing whatever the user wants (CH*100)
c     nchar : number of characters in the comment string that should be
c             printed out in output files if the object is detected;
c             maximum of 100 (I4)
c     ierr  : return code
c                  0 : nominal run, things are good
c                100 : end of model, exit after checking this object
c                -10 : could not get all orbital elements, skip object
c                -20 : something went grossly wrong, should quit
c
c The model subroutines can access files using logical unit numbers from
c 10 to 15. This range in reseved for them and won't be used by the
c drivers nor SurveySubs routines.
c
c It is good practice that when first started, the GiMeObj routine
c writes a file describing the model used, the versino and the date of
c the routine.
c
c Since this routine is called once for every object created, it needs
c to get all the required parameters once when it is called the first
c time, then save these values for future use.
c
c The following routine gives a working example of a model routine. It
c is probably worth reading it through.
c
c The survey simulator expects orbital elements with respect to ecliptic
c reference frame.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c File generated on 2013-07-01
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c First define variables so they are accessible from a Python wrapper
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      block data modeldata

      real*8
     $  lambdaN, lambdaU, epoch_m, emaxTROJ

      common /com_time/ epoch_m, lambdaN, lambdaU
      common /com_11/ emaxTROJ

      data
c One cannot use the form "sh /15.d0*drad/", so one needs to perform
c the computation before hand and give the result here, or do the
c multiplication in the program.
c
c Values of time and planet positions for all models
c
     $  lambdaN /5.489d0/,
c     $  lambdaN /5.77165d0/,
     $  lambdaU /0.06744d0/,
c     $  epoch_m /2455859.5d0/
     $  epoch_m /2453157.5d0/
      end

c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c Now define working routines and functions
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

c  *********************************************************************
c  *********************************************************************
c  Subroutine res_3_2
c
c  currently generates objects in the 3:2 mean motion resonance
c  Based on what was in plut.f from kozai_plutinos
c
c  Set up to generate a from a uniform +/- about the res center
c  e from a gaussian distribution (truncated at e=0)
c  i from a sin(i)*gaussian distribution
c  resonant amplitude from a triangle shaped thing
c
c  Should be replaced/modified as the user wishes
c      
c  Kat Volk May 2013
c
c  2013-05-14: J.-M. Petit adapted to new GiMeObj routine
c
c  2013-07-01: J.-M. Petit modified API to fit new Drivers.f
c
c  2013-09-20: J.-M. Petit changed H distributino to include divot
c
c  2014-04-25: M. Alexandersen, fixed lib-amp and a-bounds
c  2014-06-07: M. Alexandersen, fixed a-bounds again
c
c  2014-10-03: M. Alexandersen, transform to Trojans
c
c  *********************************************************************
c  *********************************************************************

      subroutine GiMeObj (filena, seed, a, e, inc, node, peri, M,
     $  epoch, h, color, gb, ph, period, amp, commen, nchar, ierr)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c This routine generates an object from a model stored in a file.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c J-M. Petit  Observatoire de Besancon
c Version 1 : May 2013
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c INPUT
c     filena: name of containing description of model, read in by the
c             model subroutine "modname" (CH)
c     seed  : Random number generator seed (I4)
c
c OUTPUT
c     a     : semimajor axis (R8)
c     e     : eccentricity (R8)
c     inc   : Inclination [rad] (R8)
c     node  : Longitude of node [rad] (R8)
c     peri  : Argument of perihelion [rad] (R8)
c     M     : Mean anomaly [rad] (R8)
c     epoch : Time of elements [JD] (R8)
c     h     : Absolute magnitude of object in 'x' band, what ever this is (R8)
c     color : Array of colors (10*R8)
c                colors(1) : g-x
c                colors(2) : r-x
c                colors(3) : i-x
c                colors(4) : z-x
c                colors(5) : u-x
c                colors(6) : V-x
c                colors(7) : B-x
c                colors(8) : R-x
c                colors(9) : I-x
c     gb    : opposition surge factor, Bowell formalism (R8)
c     ph    : phase of lightcurve at epoch [rad] (R8)
c     period: period of lightcurve [day] (R8)
c     amp   : amplitude of lightcurve [mag] (R8)
c     commen: user specified string containing whatever the user wants (CH*100)
c     nchar : number of characters in the comment string that should be
c             printed out in output files if the object is detected;
c             maximum of 100 (I4)
c     ierr  : return code (I4)
c                  0 : nominal run, things are good
c                100 : end of model, exit after checking this object
c                -10 : could not get all orbital elements, skip object
c                -20 : something went grossly wrong, should quit
c
c The user can fill the 100-character 'commen' string any way they
c wish; this comment string will be printed in the driver on the output
c line of each detection.  Examples of the comment might be resonance name
c and libration amplitude, or the name of a component in the GiMeObj model
c that the object responds to.  The nchar variable (passed back to Driver)
c allows the user to pring only the first nchar characters of this string.
c
c This routine uses logical unit 10 to access the file containing the model.
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c Set of F2PY directive to create a Python module
c
Cf2py intent(in) filena
Cf2py intent(in,out) seed
Cf2py intent(out) a
Cf2py intent(out) e
Cf2py intent(out) inc
Cf2py intent(out) node
Cf2py intent(out) peri
Cf2py intent(out) M
Cf2py intent(out) epoch
Cf2py intent(out) h
Cf2py intent(out) color
Cf2py intent(out) gb
Cf2py intent(out) ph
Cf2py intent(out) period
Cf2py intent(out) amp
Cf2py intent(out) commen
Cf2py intent(out) nchar
Cf2py intent(out) ierr
  
      implicit none

c Calling arguments
      integer*4 ierr, seed, nchar
      real*8 a, e, inc, node, peri, M, epoch, h, color(*), gb, ph,
     $  period, amp, colour(10) !MA
      character filena*(*), commen*100, distri_file*80, 
     $  objflag*5, filename*80

c Some values better set up as parameters
      integer*4
     $  lun_m,                  ! Logical unit number for model description file
     $  lun_ll                  ! Logical unit number for logging
      real*8
     $  Pi,                     ! Pi
     $  TwoPi,                  ! 2*Pi
     $  drad                    ! Degree to radian convertion: Pi/180
c Set the values
      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.0d0*Pi, drad = Pi/180.0d0,
     $  lun_m = 20, lun_ll = 21)

c Internal storage
      integer*4
     $  i,                      ! Dummy index
     $  j,                      ! Dummy index
     $  kozflag,                ! Flag for kozai resonance
     $  nparam                  ! Number of parameters for distribution
      real*8
     $  resamp,                 ! Resonance amplitude
     $  phi11,                  ! Resonant angle
     $  size_dist_one,          ! Returns H according to exponential distrib
     $  h_params(6),            ! Parameters for the H distribution
     $  params(6),              ! Parameters for some distribution
     $  incmin,                 ! Lower limit of inclination distribution
     $  incmax,                 ! Upper limit of inclination distribution
     $  emin,                   ! Lower limit of eccentricity distribution
     $  emax,                   ! Upper limit of eccentricity distribution
     $  lmin,                   ! Lower limit of libration distribution
     $  lmax,                   ! Upper limit of libration distribution
     $  epoch_m,                ! Epoch of elements [JD]
     $  random,                 ! Random number
     $  ht,                     ! Random number stored for later
     $  ran3,                   ! Random number generator
     $  gasdev,                 ! Returns a random number with gaussian
                                ! distribution of given width and center
     $  q,                      ! Perihelion distance
     $  lambdaN,                ! Neptune's mean longitude at reference epoch
     $  lambdaU,                ! Uranus' mean longitude at some reference epoch
     $  alphaU11, alphaN11, alphaso, ! frac of Uranian, Neptunian, scattering
     $  betaUH, betaUT,         ! Fraction of Uranian Horseshoes and Trojans
     $  betaNH, betaNT,         ! Fraction of Neptunian Horseshoes and Trojans
     $  ecent,                  ! mean eccentricity
     $  ew,                     ! standard deviation of eccentricity
     $  colcent,                !MA mean colour         
     $  colw,                   !MA standard deviation of colour
     $  emaxTROJ,               ! Maximum eccentricity of plutinos
     $  sg2deg,                 ! standard deviation in inclination (degree)
     $  sg2,                    ! standard deviation in inclination (radian)
     $  libampc,                ! width of the libration amp dist
     $  h0,                     ! lowest value for h distribution
     $  h1,                     ! knee value for h distribution
     $  h2,                     ! highest value for h distribution
     $  slope,                  ! slope of size distribution, before divot
     $  slope2,                 ! slope of size distribution, past divot
     $  ct,                     ! constrast factor at divot
     $  ckc,                    ! nb/(nb+ns)
     $  nb,                     ! number of ojbect brighter than divot
     $  ns,                     ! number of ojbect fainter than divot
     $  h0s10,                  ! normalizing factor at low end of H-distrib
     $  h1s10,                  ! normalizing factor at high end of H-distrib
     $  r,                      ! Distance of object to Sun
     $  x,                      ! X cartesian coordinate of object
     $  y,                      ! Y cartesian coordinate of object
     $  z,                      ! Z cartesian coordinate of object
     $  rayleigh,               ! Name of function for other distributions
     $  onecomp                 ! Name of function for inclination distribution
      logical
     $  first,                  ! Tells if first call to routine
     $  log                     ! True if we want to log data
      integer*4  nlines, n_elements

c Centaur orb. elems. to be read in from file
      real*8
     $  am(306000), em(306000), im(306000), fbin(306000), fcum(306000)
      integer*2
     $  ac(306000), ec(306000), ic(306000)

c Place some variables in common block so they can be accessed directly
c by a Python program.
      common /com_11/ emaxTROJ
      common /com_time/ epoch_m, lambdaN, lambdaU

c Ran3, size_dist_one and gasdev are functions defined outside
c this routine.
      external ran3, size_dist_one, gasdev, onecomp, onecompjmp, 
     $  rayleigh

c Sets initial values
      data first /.true./

      save

c This is the first call
      if (first) then
c Define values for lightcurve and opposition surge effect
         gb = -0.12d0
         ph = 0.d0
         period = 0.6d0
         amp = 0.d0
c Reads in other parameters describing the model.
         open(unit=lun_m, file=filena, status="unknown", err=1000)
         read(lun_m,*) sg2deg
         sg2=sg2deg*drad
         read(lun_m,*) h0
         read(lun_m,*) h1
         read(lun_m,*) h2
         read(lun_m,*) slope
         read(lun_m,*) slope2
         h_params(1) = h0
         h_params(2) = h1
         h_params(3) = h2
         h_params(4) = slope
         h_params(5) = slope2
         read(lun_m,*) ct
         nb = (10.d0**(slope*h1) - 10.d0**(slope*h0))/slope
         ns =10.d0**((slope-slope2)*h1)/slope2
     $     *(10.d0**(slope2*h2) - 10.d0**(slope2*h1))/ct
         ckc = nb/(nb + ns)
         read(lun_m,*) ew
cMA         read (lun_m, *) (color(i),i=1,10)
         read (lun_m, *) alphaU11, alphaN11, alphaso
         read (lun_m, *) betaUH, betaUT
         read (lun_m, *) betaNH, betaNT
         read (lun_m, *) n_elements
         read (lun_m, *) filename
         read (lun_m, *) (colour(i),i=1,10)
         read (lun_m, *) log
         close(lun_m)

c Read in the residence file to draw SO's from.
         open (unit=2730, file=filename, status='old', err=1000)
         do j=1, n_elements
          read(2730,*) ac(j), ec(j), ic(j),
     $               am(j), em(j), im(j),
     $               fbin(j), fcum(j)
         enddo
         close(2730)

c Writes a file describing the model that was used.
         open (unit=lun_ll, file='ModelUsed', access='append',
     $     status='unknown')
         write (lun_ll, '(a)')
     $     '# Model of Scattering objects and temporary 1:1 resonance, 
     $     version 2.0, 2014-10-08'
         write (lun_ll, '(''#'')')
         write (lun_ll, '(a,5(1x,f6.2))')
     $        '# H limits and slopes:', (h_params(i),i=1,5)
         write (lun_ll, '(a,3(1x,f6.4))')
     $        '# Scattering, Uranian, Neptunian fraction:', 
     $        alphaso, alphaU11, alphaN11
         write (lun_ll, '(a,3(1x,f6.4))')
     $        '# Uranian Horseshoe, Trojan, Quasi-sat fraction:', 
     $        betaUH, betaUT, 1-betaUH-betaUT
         write (lun_ll, '(a,3(1x,f6.4))')
     $        '# Neptunian Horseshoe, Trojan, Quasi-sat fraction:', 
     $        betaNH, betaNT, 1-betaNH-betaNT
         write (lun_ll, '(a,1(1x,f6.2),a)')
     $        '# Width of inc:              ', sg2deg
         write (lun_ll, '(a,1(1x,f6.3),a)')
     $        '# Width of ecc:               ', ew
         write (lun_ll, '(''#'')')
         close (lun_ll)
c Change "first" so this is not called anymore
         first = .false.
      end if

cMA Very simple randomising of colours with a bimodal distribution 
cMA as seen from centaurs and scatterers in Piexinho et al 2012. 
cMA Implemented by M. Alexandersen, only appropriate for g-r
         colw=0.1
         random=ran3(seed)
         if (random .lt. 0.5) then
           colcent=0.55
           color(1)=gasdev(colcent,colw,seed)
         else
           colcent=0.93
           color(1)=gasdev(colcent,colw,seed)
         endif
         do j = 2,10
            color(j) = colour(j)
         end do

c Pick object type.
 1040 continue
      objflag='??'
      random=ran3(seed)
      if (random .lt. alphaN11) goto 1050           ! Object is in N 1:1
      if (random .lt. alphaN11+alphaU11) goto 1060  ! Object is in U 1:1
      if (random .lt. alphaN11+alphaU11+alphaso) goto 1070  ! ScatteringObj
      if (alphaN11+alphaU11+alphaso .gt. 1.0d0) then
              write(6,*) 'sum >= 1, stopping'
              stop
      endif
      goto 1040 ! This shouldn't happen

c Neptunian Coorbital section
 1050 continue
c Pull an eccentricity      
      random=ran3(seed)
      params(1)=ew
      nparam=1
      emin=(0.0d0)
      emax=(0.30d0)
      call incdism (seed,nparam,params,emin,emax,e,2,ierr,onecompjmp)

c pick 'a'
      random=ran3(seed)
      a= 30.10 + (random-0.5)*2.0 * 0.20

c cutoff if Uranus-approaching. (won't do for now)
      q = a*(1.0d0-e)
c      if (q .lt. 22.d0) goto 1050

c Pull an inclination
      random = ran3(seed)
      params(1) = sg2
      nparam = 1
      incmin = (0.0d0*drad)
      incmax = (90.0d0*drad)
      call incdism (seed, nparam, params, incmin, incmax, inc,
     $  1, ierr, onecomp)

c node and M picked randomly.  
      random=ran3(seed)
      node = random*TwoPi
      random=ran3(seed)
      M = random*TwoPi

c resonant argument and libration amplitude are different for horse/troj
         ht=ran3(seed)
         if (ht .lt. betaNT) then
           objflag='84 NT'
c Need random resonant argument, with libration amplitude
 1054      continue
           resamp = ran3(seed)*70.0d0
           random = ran3(seed)
           if (random .gt. resamp/70.0d0) goto 1054
c Now choose phi11
           random=ran3(seed)
           phi11 = 60.0d0*drad +
     $              ((resamp/10)**1.8)*drad +
     $              sin(TwoPi*random)*resamp*drad
c Second term in the phi11 is to offset the centre.
           random=ran3(seed)
           if (random .ge. 0.5) then 
             objflag='85 NT'
             phi11 = 360.0d0*drad - phi11 ! object trailing
           endif
         else if (ht .lt. betaNT+betaNH) then
           objflag='83 NH'
 1055      continue
           resamp = ran3(seed)*25.0d0
           random = ran3(seed)
           if (random .gt. resamp/25.0d0) goto 1055
           resamp=resamp+135.0d0
           random=ran3(seed)
           phi11 = 180.0d0*drad + sin(TwoPi*random)*resamp*drad
         else
           objflag='80 NQ'
 1056      continue
           resamp = ran3(seed)*45.0d0+5.0d0
           random=ran3(seed)
           phi11 = 0.0d0*drad + sin(TwoPi*random)*resamp*drad
         endif
         call zero2pi(phi11)
c Set argument of pericenter based on resonant angle
         peri = phi11 - 1.d0*m +1.d0*lambdaN - node
         call zero2pi(peri)
         goto 1999

c Uranian Coorbital section
 1060 continue
c Pull an eccentricity      
      random=ran3(seed)
      params(1)=ew
      nparam=1
      emin=(0.0d0)
      emax=(0.30d0)
      call incdism (seed,nparam,params,emin,emax,e,4,ierr,onecompjmp)

c pick 'a'
      random=ran3(seed)
      a= 19.23 + (random-0.5)*2.0 * 0.20

c cutoff if Saturn-approaching. (won't do for now)
      q = a*(1.0d0-e)
c      if (q .lt. 12.d0) goto 1061

c Pull an inclination
      random = ran3(seed)
      params(1) = sg2
      nparam = 1
      incmin = (0.0d0*drad)
      incmax = (90.0d0*drad)
      call incdism (seed, nparam, params, incmin, incmax, inc,
     $  3, ierr, onecomp)

c node and M picked randomly.  
      random=ran3(seed)
      node = random*TwoPi
      random=ran3(seed)
      M = random*TwoPi

c resonant argument and libration amplitude are different for horse/troj
         ht=ran3(seed)
         if (ht .lt. betaUT) then
           objflag='74 UT'
c Need random resonant argument, with libration amplitude
 1064      continue
           resamp = ran3(seed)*70.0d0
           random = ran3(seed)
           if (random .gt. resamp/70.0d0) goto 1064
c Now choose phi11
           random=ran3(seed)
           phi11 = 60.0d0*drad +
     $              ((resamp/10)**1.8)*drad +
     $              sin(TwoPi*random)*resamp*drad
c Second term in the phi11 is to offset the centre.
           random=ran3(seed)
           if (random .ge. 0.5) then 
             objflag='75 UT'
             phi11 = 360.0d0*drad - phi11 ! object trailing
           endif
         else if (ht .lt. betaUT+betaUH) then
           objflag='73 UH'
 1065      continue
           resamp = ran3(seed)*25.0d0
           random = ran3(seed)
           if (random .gt. resamp/25.0d0) goto 1065
           resamp=resamp+135.0d0
           random=ran3(seed)
           phi11 = 180.0d0*drad + sin(TwoPi*random)*resamp*drad
         else
           objflag='70 NQ'
 1066      continue
           resamp = ran3(seed)*45.0d0+5.0d0
           random=ran3(seed)
           phi11 = 0.0d0*drad + sin(TwoPi*random)*resamp*drad
         endif
c Set argument of pericenter based on resonant angle
         call zero2pi(phi11)
         peri = phi11 - 1.d0*m +1.d0*lambdaN - node
         call zero2pi(peri)
         goto 1999

c        Generate a scattering object/Centaur
 1070    continue
         objflag='99 SO'
         resamp=0.0d0
         phi11=0.0

c Pick a line from lookup table
         random=ran3(seed)
c This is my terrible attempt at doing in fortran
c what a "where" function can do in IDL: 
         j=1
         do while (fcum(j) .lt. random)
            j=j+1
         enddo
c Now j is the index of the cell that I want elements from
c I select the elements and fuzz them by the size of the cells
         a=am(j)+ran3(seed)*0.5d0-0.25d0
         e=em(j)+ran3(seed)*0.02d0-0.01d0
         inc=(im(j)+ran3(seed)*2.0d0-1.0d0)*drad
c  pick node, argperi & M randomly
         node=ran3(seed)*TwoPi
         peri=ran3(seed)*TwoPi
         M=ran3(seed)*TwoPi
c done
         goto 1999

 1999 continue
c
c H-mag distribution
      random = ran3(seed)
      if (random .lt. ckc) then
         slope = h_params(4)
         h0 = h_params(1)
         h1 = h_params(2)
      else
         slope = h_params(5)
         h0 = h_params(2)
         h1 = h_params(3)
      end if
      h0s10 = 10.d0**(slope*h0)
      h1s10 = 10.d0**(slope*h1)
      random = ran3(seed)
      h = log10(random*(h1s10 - h0s10) + h0s10)/slope

      epoch = epoch_m

c Stores informations on resonance
      commen = ' '
      write (commen(2:12), '(f10.5,1x)') resamp ! resamp is in degrs everywhere
      write (commen(14:24), '(f10.6,1x)') phi11/drad
      write (commen(25:30), '(1x,a)') objflag
      nchar = 35

c Stores object if user requested
      if (log) then
         call pos_cart(a,e,inc,node,peri,M,x,y,z)
         r = sqrt(x*x + y*y + z*z)
         open (unit=lun_ll, file='ModelUsed', access='append',
     $     status='old')
cMA edited the next few lines to get output more like detect and track file.
         write(lun_ll,101) a,e,inc/drad,a*(1.d0-e),r,M/drad,a/a,h,
     $     color(1), x, y, z, commen(1:nchar), node/drad, peri/drad
         close (lun_ll)
 101     format(7(f8.4,1x),f6.2,1x,f5.2,3(1x,f7.3),1x,a,2(1xf8.4))
      end if
      ierr = 0

      return

 1000 continue
      ierr = -20
      return

      end

c--------------------------------------------
c   gasdev
c********************************************************************
c  Given a center, width and seed return a value drawn from a gaussian

      FUNCTION gasdev(x0,sigma,rs)

      implicit none
      real*8 x0, sigma
      integer*4 iset, rs
      real*8 v1, v2, rsq, gasdev, gset
      real*8 fac, random, ran3
      external ran3
      SAVE iset,gset
      
      if  (iset.eq.0) then
 12      continue
             random = ran3(rs)
             v1= 2d0*random - 1d0
             random = ran3(rs)
             v2= 2d0*random - 1d0
             rsq=v1*v1+v2*v2
         if (rsq.ge.1.0) goto 12 

         fac=sqrt(-2d0*log(rsq)/rsq)
         gset=v1*fac
         iset=1
         gasdev=v2*fac*sigma+x0
         return 
      else
         iset=0
         gasdev=gset*sigma+x0
         return
      endif
      END

      real*8 function rayleigh (nparam, param, inc)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c This routine returns the unnormalized "probability"
c density of a Rayleigh function.
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c
c M. Alexandersen, UBC
c Version 1 : October 2014
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
c INPUT
c     nparam: Number of parameters (I4)
c     param : Parameters (n*R8)
c     inc   : Inclination (R8)
c
c OUPUT
c     rayleigh: Value of the inclination (R8)
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
Cf2py intent(in) nparam
Cf2py intent(in), depend(nparam) :: param
Cf2py intent(in) inc

      implicit none

      integer
     $  nparam

      real*8
     $  param(*), inc

      real*8
     $  Pi, TwoPi

      parameter
     $  (Pi = 3.141592653589793d0, TwoPi = 2.0d0*Pi)

      real*8
     $  fe, s1, angle, t1, t3

      if (nparam .ne. 1) stop
      s1 = param(1)
      t1 = 2.*s1**2
      t3 = -inc**2
      fe = exp(t3/t1)
      rayleigh = (inc)*fe

      return
      end


      include 'ModelUtils.f'
