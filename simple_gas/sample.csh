#! /bin/csh -f
# syntax: ./sample.csh <model_filename> <data_filename>
# where <model_filename> is the prefix of the full-resolution model cube
# 	<data_filename> is the prefix for the data cube
#	<nchan+1> is the # of velocity channels you want output (+1)

start:

# SAMPLING / VISIBILITY CALCULATIONS
# ----------------------------------

# convert model image to MIRIAD format
rm -rf $1.im
fits  in=$1.fits op=xyin out=$1.im

# update model image header
puthd in=$1.im/epoch  value=`gethd in=$2.cm/epoch`
puthd in=$1.im/crval1 value=`gethd in=$2.cm/crval1`
puthd in=$1.im/crval2 value=`gethd in=$2.cm/crval2`
puthd in=$1.im/restfreq value=345.79599 type=double

# print warnings
echo 'BE VIGILANT: setting rest frequency by hand!'

# resample model in the velocity plane
echo 'BE VIGILANT: skipping any regridding in velocity plane'
#rm -rf sky.vis
#regrid in=$1.im tin=$2.cm out=sky.vis axes=3

# resample model onto observed (u,v) tracks
rm -rf $1.model.vis $1.model.vis.fits
uvmodel vis=$2.vis model=$1.im out=$1.model.vis options=replace
fits in=$1.model.vis op=uvout out=$1.model.vis.fits

# skip imaging?
#goto end

# IMAGING
# -------

# MODEL IMAGE

# Fourier inversion / make a dirty image
rm -rf $1.model.map $1.model.bm
invert vis=$1.model.vis map=$1.model.map beam=$1.model.bm \
       imsize=1024 cell=0.05 robust=2

# deconvolution / clean
rm -rf $1.model.clean
clean map=$1.model.map beam=$1.model.bm out=$1.model.clean \
      niters=100000 cutoff=0.04 region='arcsec,box(-7,-7,7,7)'

# restore with synthesized beam
rm -rf $1.model.cm
restor model=$1.model.clean beam=$1.model.bm map=$1.model.map out=$1.model.cm 

# output synthesized model image to fits
rm -rf $1.model.fits
fits in=$1.model.cm op=xyout out=$1.model.fits

clean:
rm -rf $1.model.map $1.model.bm
rm -rf $1.model.clean
rm -rf $1.model.cm

end:
