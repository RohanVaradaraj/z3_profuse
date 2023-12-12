evalglobal = TRUE

library(ProFuse)
library(ProSpect)
library(ProFit)
library(ProFound)
library(Rfits)
library(ParmOff)
library(magicaxis)
library(celestial)
library(imager)

redshift = 0.0447
data('BC03lr')
data('Dale_NormTot')
data('AGN_UnOb_Sparse')
data('Dale_M2L_func')
filters=c('u_VST', 'g_VST', 'r_VST', 'i_VST', 'Z_VISTA',
          'Y_VISTA', 'J_VISTA', 'H_VISTA', 'Ks_VISTA')
filtout={}
for(i in filters){filtout=c(filtout,list(approxfun(getfilt(i))))}

loc = c(1200,480)
cut = -299:300

cenwaves = cenwave[match(c('u_VST', 'g_VST', 'r_VST', 'i_VST', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'Ks_VISTA'), cenwave$filter),'cenwave']

agemax = 13.3e9 - cosdistTravelTime(z=redshift, H0 = 67.8, OmegaM = 0.308)*1e9

image_list = list(
  u = Rfits_read_image(system.file("extdata", 'MultiBand/u.fits', package="ProFound"),ext=2)$imDat[loc[1] + cut, loc[2] + cut],
  g = Rfits_read_image(system.file("extdata", 'MultiBand/g.fits', package="ProFound"),ext=2)$imDat[loc[1] + cut, loc[2] + cut],
  r = Rfits_read_image(system.file("extdata", 'MultiBand/r.fits', package="ProFound"),ext=2)$imDat[loc[1] + cut, loc[2] + cut],
  i = Rfits_read_image(system.file("extdata", 'MultiBand/i.fits', package="ProFound"),ext=2)$imDat[loc[1] + cut, loc[2] + cut],
  Z = Rfits_read_image(system.file("extdata", 'MultiBand/Z.fits', package="ProFound"),ext=2)$imDat[loc[1] + cut, loc[2] + cut],
  Y = Rfits_read_image(system.file("extdata", 'MultiBand/Y.fits', package="ProFound"),ext=2)$imDat[loc[1] + cut, loc[2] + cut],
  J = Rfits_read_image(system.file("extdata", 'MultiBand/J.fits', package="ProFound"),ext=2)$imDat[loc[1] + cut, loc[2] + cut],
  H = Rfits_read_image(system.file("extdata", 'MultiBand/H.fits', package="ProFound"),ext=2)$imDat[loc[1] + cut, loc[2] + cut],
  Ks = Rfits_read_image(system.file("extdata", 'MultiBand/Ks.fits', package="ProFound"),ext=2)$imDat[loc[1] + cut, loc[2] + cut]
)

MF2F = profuseMultiBandFound2Fit(image_list=image_list,
                                 magzero=c(0,0,0,0,30,30,30,30,30),
                                 tolerance=20,
                                 parm_global = c("sersic.xcen1", "sersic.ycen1",
                                                 "sersic.re1", "sersic.re2",
                                                 "sersic.ang2", "sersic.axrat2"),
                                 parm_ProSpect = list(
                                   mSFR_1 = 0,
                                   mpeak_1 = 10,
                                   mperiod_1 = 0.3,
                                   mskew_1 = 0,
                                   #tau_screen_1 = -0.5,
                                   #alpha_SF_screen_1 = 3,
                                   Zfinal_1 = -2,
                                   
                                   mSFR_2 = 0,
                                   mpeak_2 = 1,
                                   mperiod_2 = 0.3,
                                   mskew_2 = 0,
                                   #tau_screen_2 = -0.5,
                                   #alpha_SF_screen_2 = 3
                                   Zfinal_2 = -2
                                 ),
                                 logged_ProSpect = c(
                                   TRUE,                              #         mSFR_1 = 0,
                                   FALSE,                             #         mpeak_1 = 5,
                                   TRUE,                              #        mperiod_1 = 0.3,
                                   FALSE,                             #         mskew_1 = 0,
                                   #TRUE,                              #        tau_screen_1 = -0.5,
                                   #FALSE,                           #         alpha_SF_screen_1 = 3,
                                   TRUE,                              #        Zfinal_1 = -2,
                                   TRUE,                              #        mSFR_2 = 0,
                                   FALSE,                             #         mpeak_2 = 5,
                                   TRUE,                              #        mperiod_2 = 0.3,
                                   FALSE,                             #         mskew_2 = 0,
                                   #TRUE,                              #        tau_screen_2 = -0.5,
                                   #FALSE                              #        alpha_SF_screen_2 = 3
                                   TRUE                              #        Zfinal_2 = -2,
                                 ),
                                 intervals_ProSpect = list(
                                   lo=c(
                                     -4,                              #       mSFR_1 = 0,
                                     0,                               #       mpeak_1 = 5,
                                     -0.5,                            #          mperiod_1 = 0.3,
                                     -1,                              #        mskew_1 = 0,
                                     #-2.5,                            #          tau_screen_1 = -0.5,
                                     #0,                               #       alpha_SF_screen_1 = 3,
                                     -4,                              #        Zfinal_1 = -2,
                                     -4,                              #        mSFR_2 = 0,
                                     0,                               #       mpeak_2 = 5,
                                     -0.5,                            #          mperiod_2 = 0.3,
                                     -1,                              #        mskew_2 = 0,
                                     #-2.5,                            #          tau_screen_2 = -0.5,
                                     #0                                #      alpha_SF_screen_2 = 3
                                     -4                               #        Zfinal_2 = -2,
                                   ),
                                   hi=c(
                                     3,                               #       mSFR_1 = 0,
                                     10,                              #        mpeak_1 = 5,
                                     1,                               #       mperiod_1 = 0.3,
                                     1,                               #       mskew_1 = 0,
                                     #1,                               #       tau_screen_1 = -0.5,
                                     #4,                               #       alpha_SF_screen_1 = 3,
                                     -1.3,                              #        Zfinal_1 = -2,
                                     3,                               #       mSFR_2 = 0,
                                     10,                              #        mpeak_2 = 5,
                                     1,                               #       mperiod_2 = 0.3,
                                     1,                               #       mskew_2 = 0,
                                     #1,                               #       tau_screen_2 = -0.5,
                                     #4                                #      alpha_SF_screen_2 = 3
                                     -1.3                              #        Zfinal_2 = -2,
                                   )
                                 ),
                                 data_ProSpect = list(massfunc=massfunc_snorm_trunc,
                                                      speclib=BC03lr,
                                                      Dale=Dale_NormTot,
                                                      filtout=filtout,
                                                      z=redshift,
                                                      Z=Zfunc_massmap_lin,
                                                      agemax=agemax
                                 )
)

highfit = profuseMultiBandDoFit(MF2F=MF2F)

magcurve(massfunc_snorm_trunc(age=x,mSFR=10^highfit$parm["mSFR_1"],mpeak=highfit$parm["mpeak_1"],mperiod=10^highfit$parm["mperiod_1"],
                              mskew=highfit$parm["mskew_1"], magemax=agemax/1e9),0,13e9,add=FALSE,col='red', ylim=c(0,10),xlab='Age (Yr)', ylab='SFR (Msol / Yr)')
magcurve(massfunc_snorm_trunc(age=x,mSFR=10^highfit$parm["mSFR_2"],mpeak=highfit$parm["mpeak_2"],mperiod=10^highfit$parm["mperiod_2"],
                              mskew=highfit$parm["mskew_2"], magemax=agemax/1e9),0,13e9,add=TRUE,col='blue')
legend('topright', legend=c('Bulge', 'Disk'), col=c('red', 'blue'), lty=1)

profitLikeModel(highfit$parm, MF2F, makeplots = TRUE)

