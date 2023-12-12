evalglobal = TRUE

library(ProFuse)
library(ProSpect)
library(ProFit)
library(ProFound)
library(Rfits)
library(ParmOff)
library(magicaxis)
library(celestial)

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

print('OK')

