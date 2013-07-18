#!/bin/csh
set name='HD163296.CO32.regridded.cen15'

# Fourier inversion / make a dirty image
rm -rf $name.map $name.bm
invert vis=$name.vis map=$name.map beam=$name.bm \
       imsize=1024 cell=0.05 robust=2

# deconvolution / clean
rm -rf $name.clean
clean map=$name.map beam=$name.bm out=$name.clean \
      niters=100000 cutoff=0.04 region='arcsec,box(-7,-7,7,7)'

# restore with synthesized beam
rm -rf $name.cm
restor model=$name.clean beam=$name.bm map=$name.map out=$name.cm 

# output synthesized model image to fits
rm -rf $name.fits
fits in=$name.cm op=xyout out=$name.fits

clean:
rm -rf $name.map $name.bm

