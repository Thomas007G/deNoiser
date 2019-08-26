#!/usr/bin/env python

# v4 21 Jun 2013: now correctly deals with fibre offset and includes
#                 choice of Gaussian or Moffat profile
# v4.1 28 jan 2014: updated efficiencies and central obstruction
# v4.1.1 30 jan 2014: updated blue LR efficiencies
# v4.1.2 13 mar 2014: "clean" version to remove Kelson-specific code
# v4.1.3 27 jan 2015: updated typical image FWHM to 0.95" (includes typical
#                     seeing + PRI contribution); added default exposure time
#                     of 23 minutes (x3=69 minutes) at surveys start
# v4.1.3.1 29 jan 2015: updated default exposure time
#                     to 18 minutes (x3=55 minutes) at surveys start
# v4.1.3.2 20 may 2016: updated efficiencies following 17 nov 2015 updates;
#                     typical exposure times now 17 minutes
# v4.2     12 jan 2018: updated efficiencies following 22 dec 2017 almost-as-built throughput
#                       (WEAVE-SYS-007 v1.8)
#                       added green HR (=blueHR2) VPH grating
# v4.2.3   15 feb 2018: added blurring due to PRI image quality; "seeing" is now as expected (and not total image FWHM); corrected default seeing to 0.75"
# follows Chris Benn's signal.f to compute signal-to-noise for WEAVE

from sys import argv
import signal2noise

help="""Exposure time calculator for WEAVE.  Run as
% ./signalWEAVE <optional switches> <magnitude or surface brightness> <exposure time in seconds>

Notation: <> represents a number required on input

If no exposure time is given, the ETC assumes a single 18 minute exposure.

Switches:

* Instrument Mode *
choose ONE of
-blueLR [default]
-redLR
-blueHR
-greenHR
-redHR

* Object/observing conditions parameters *
-band <filter of magnitude or surface brightness> [default: V] (choices: UBVRI)
-X <airmass> [default: 1.2]
-sky <surface brightness of sky between the lines> [default: 22.7 mag/sq. arcsec in B band]
-skyband <band of interest for sky background> [default: B] (choices: UBVRI)
-profile <seeing profile> [default: Moffat] (choices: Gaussian or Moffat)
 -beta <beta of Moffat profile> [default: 2.5] (only if profile is Moffat)
-FWHM <full width at half maximum of seeing disk> [default: 0.8]
-readnoise <RN in e-> [default: 2.5]
-dark <dark current in e-/hour> [default: 0.]
-offset <offset of fiber center from object center in arcseconds> [default: 0.1]
-sb (switch to set input as surface brightness in mag./sq. arcsec instead of point source magnitudes)

Return values should be self-explanatory."""

if len(argv)==1:
   print help
   exit()

# new in v4
# PSF type
argv,profile=signal2noise.rdarg(argv,'-profile',None,'Moffat')
betam=None
if profile=='Moffat':
    argv,betam=signal2noise.rdarg(argv,'-betam',float,2.5)
elif profile!='Gaussian':
    print 'Profile',profile,'not yet implemented. Stopping.'
    sys.exit()

# S/N per Angstrom?  Default is S/N per resolution element
argv,ang=signal2noise.rdarg(argv,'-ang',single=1)
if ang:
    snrper='SNR'
else:
    snrper='SNRres'
# Fiber injection f/ratio
argv,fFP=signal2noise.rdarg(argv,'-fFocalPlane',float,3.2)
# Collimator f/ratio
argv,fcol=signal2noise.rdarg(argv,'-fCollimator',float,3.1)
# Camera f/ratio
argv,fcam=signal2noise.rdarg(argv,'-fCamera',float,1.8)
# resolution
argv,res=signal2noise.rdarg(argv,'-resolution',float,5000.)
# central wavelength
argv,cwave=signal2noise.rdarg(argv,'-cwave',float,4900.)

# Arm
argv,red=signal2noise.rdarg(argv,'-redArm',None,single=1)
argv,blue=signal2noise.rdarg(argv,'-blueArm',default=1,single=1)
if blue:
    blue=1
    red=None
    pfEff=0.803
else:
    blue=None
    pfEff=0.811

# efficiencies: detector, fibers, spectrograph, PFcorr - red LR defaults
QE=0.859
fiberEff=0.868
specEff=0.684

# shortcuts
argv,blueLR=signal2noise.rdarg(argv,'-blueLR',default=None,single=1)
argv,redLR=signal2noise.rdarg(argv,'-redLR',default=None,single=1)
argv,blueHR=signal2noise.rdarg(argv,'-blueHR',default=None,single=1)
argv,greenHR=signal2noise.rdarg(argv,'-greenHR',default=None,single=1)
argv,redHR=signal2noise.rdarg(argv,'-redHR',default=None,single=1)
##if not blueLR and not redLR and not blueHR and not redHR:
##   blueLR=1
if blueLR: # note: values are averages over 400-590 nm
   print 'Mode: blueLR'
   res=5750.
   blue=1
   red=None
   cwave=4900.
   pfEff=0.803
   fiberEff=0.784
   specEff=0.545
   QE=0.924
elif redLR: # note: values are averages over 610-900 nm
   print 'Mode: redLR'
   res=5750.
   blue=None
   red=1
   cwave=7950.
   pfEff=0.810
   fiberEff=0.852
   specEff=0.579
   QE=0.805
elif blueHR: # average over 413-455 nm
   print 'Mode: blueHR'
   res=21000.
   blue=1
   red=None
   cwave=4250.
   pfEff=0.799
   fiberEff=0.699
   specEff=0.363 # note that the response is *highly* peaked
   QE=0.923
elif greenHR: # average over 483-533 nm
   print 'Mode: greenHR'
   res=21000.
   green=1
   red=None
   cwave=5130.
   pfEff=0.803
   fiberEff=0.794
   specEff=0.470 # note that the response is *highly* peaked
   QE=0.927
elif redHR: # average over 608-680 nm
   print 'Mode: redHR'
   res=21000.
   blue=None
   cwave=6450.
   red=1
   pfEff=0.807
   fiberEff=0.837
   specEff=0.462
   QE=0.947


argv,pfEff=signal2noise.rdarg(argv,'-PF',float,pfEff)
argv,QE=signal2noise.rdarg(argv,'-QE',float,QE)
argv,fiberEff=signal2noise.rdarg(argv,'-fiber',float,fiberEff)
argv,specEff=signal2noise.rdarg(argv,'-spec',float,specEff)

argv,eff=signal2noise.rdarg(argv,'-eff',float,None)

# airmass
argv,airmass=signal2noise.rdarg(argv,'-X',float,1.2)

# bandpass for sky brightness
argv,skyband=signal2noise.rdarg(argv,'-skyband',None,'B')
# sky brightness, in mag/sq. arcsec
argv,skysb=signal2noise.rdarg(argv,'-sky',float,22.7)
# default is sky brightness *between lines* on dark night in B, from SIGNAL
# sky flux in photons/s/m^2/Ang/arcsec^2 from Ellis & Bland-Hawthorne (2008)
# at 900 nm... note that this is ~6.4 x brighter than the above estimate...
#sky=0.08

# band of interest
argv,band=signal2noise.rdarg(argv,'-band',None,'V')
# seeing FWHM
argv,seeing=signal2noise.rdarg(argv,'-FWHM',float,0.75)
seeing=signal2noise.totalSeeing(seeing)
# fiber diameter
argv,fiberD=signal2noise.rdarg(argv,'-diameter',float,1.3045)

# detector read noise in e-
argv,rn=signal2noise.rdarg(argv,'-readnoise',float,2.5)
# dark current/hour
argv,dark=signal2noise.rdarg(argv,'-dark',float,0.)

# surface brightness (otherwise default to point source)?
argv,sb=signal2noise.rdarg(argv,'-sb',None,single=1)

# output radial velocity error?
argv,rv=signal2noise.rdarg(argv,'-rv',None,single=1)
argv,rvscale=signal2noise.rdarg(argv,'-rvscale',float,0.6)

# offset of fiber from object center
# UPDATE 01.12.2012: now assume offset of 0.1 arcsec
argv,offset=signal2noise.rdarg(argv,'-offset',float,0.1)

# magnitude (point source: magnitude; surface brightness: magnitude/sq. arcsec)
mag=float(argv[1])
# exposure time
if len(argv)>2:
   time=float(argv[2])
else:
   # typical exposure time for single exposure, in seconds
   time=17.*60.

S=signal2noise.signal(QE,fiberEff,specEff,pfEff,res,offset,fiberD,fFP,fcol,
                      fcam,cwave=cwave,profile=profile,betam=betam)
snr=S.S2N(time,mag,band,airmass,fiberD,seeing,rn,dark,eff,sb,skysb,skyband)
s2n=snr['SNR']

print 'Spectral resolving power R = %5d' % (S.res)
print 'Average dispersion (Ang/pixel) = %8.4f' % (S.disp)
print 'Number of pixels/fiber along slit = %5.2f' % (S.fiberCCD)
print 'Resolution element (Ang) = %6.3f' % (S.fiberCCD*S.disp)
print 'Efficiency = %4.2f Effective area = %5.2f m^2' % \
      (snr['efficiency'],snr['effectivearea'])
print 'Object photons/pixel = %d sky photons (between lines)/pixel = %d' % \
      (snr['objectphotons'],snr['skyphotons'])
print '(both integrated over spatial direction)'
print 'S/N/Ang = %8.2f S/N/resolution element = %8.2f' % (s2n,snr['SNRres'])
print 'S/N/pix = %8.2f' % (snr['SNRpix'])
if rv:
    sigmarv=S.RVaccuracy(s2n,scale=rvscale)
    print 'RV error (km/s) = %8.2f' % (sigmarv)

