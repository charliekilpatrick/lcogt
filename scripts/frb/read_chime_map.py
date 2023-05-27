import numpy as np
from astropy.io import fits
import h5py
from astropy_healpix import HEALPix
import healpy as hp
import os

filename = 'FRB20190617B_localization.h5'

f = h5py.File(filename,'r')

nside = f['healpix'].attrs['nside']
size = hp.nside2npix(nside)

data = np.zeros(size)

ipix = f['healpix']['ipix']
cl = f['healpix']['CL']

print(ipix)
print(cl)

for i,c in zip(ipix, cl):
    # flip cl so all values with 1-0.XX are outside the XX% probability interval
    data[i]=1.0-c

if os.path.exists('test.fits'):
    os.remove('test.fits')

hp.fitsfunc.write_map('test.fits', data, nest=True, coord='C',
    overwrite=True)
