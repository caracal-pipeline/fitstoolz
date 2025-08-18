import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import unittest
from fitstoolz import LOG as log
import os
import uuid
from fitstoolz.reader import FitsData


class TestExapndFits(unittest.TestCase):
    def setUp(self):
        
        self.ra0 = 0.0
        self.dec0 = np.deg2rad(-30.0)
        self.ncorr = 2
        self.start_freq = 1.3e9
        self.nchan = 5
        self.dfreq = 1e6

        # image parameters
        self.img_size = 256
        self.cell_size = 3e-6  # arcsec
        
        
        self.test_files = []
        # Set up logging level
        self.original_log_level = log.level
        
    def tearDown(self):
        """Clean up after each test method runs."""
        # Remove any temporary files created
        for file in self.test_files:
            if os.path.exists(file):
                os.remove(file)
        
        # Reset logging level
        log.setLevel(self.original_log_level)


    def test_fits_predicting_all_stokes_linear_basis(self):
        """
        Test visibility prediction from FITS images of all Stokes parameters, ncorr = 4
        Validates:
            - Output shape of visibilities
            - XX = I + Q
            - XY = U + iV
            - YX = U - iV
            - YY = I - Q
        """
        self.ncorr = 4
        # the numbers below are unphysical—they are just for testing the computation
        stokes_params = [('I', 1.0), ('Q', 2.0), ('U', 3.0), ('V', 4.0)]

        filenames = []        
        for stokes in stokes_params:
            wcs = WCS(naxis=3)
            wcs.wcs.ctype = ['RA---SIN', 'DEC--SIN', 'FREQ']
            wcs.wcs.cdelt = np.array([-self.cell_size/3600, self.cell_size/3600, self.dfreq]) # pixel scale in deg
            wcs.wcs.crpix = [self.img_size//2, self.img_size//2, 1] # reference pixel
            wcs.wcs.crval = [np.rad2deg(self.ra0), np.rad2deg(self.dec0), self.start_freq] # reference pixel RA and Dec in deg
            
            
            # make header
            header = wcs.to_header()
            header['BUNIT'] = 'Jy'
            
            # make image
            image = np.zeros((self.nchan, self.img_size, self.img_size))
            image[:, self.img_size//2, self.img_size//2] = stokes[1] # put a point source at the center
            # write to FITS file
            hdu = fits.PrimaryHDU(image, header=header)
            test_filename = f'test_{uuid.uuid4()}_{stokes[0]}.fits'
            self.test_files.append(test_filename)
            hdu.writeto(test_filename, overwrite=True)
            
        
            filenames.append(test_filename)
        
        fname = filenames[0]
        fnames_quv = filenames[1:]
        myfits = FitsData(fname)
        myfits.add_axis("STOKES", 0, "stokes", axis_grid=[1], attrs=dict(ref_pixel=0, pixel_size=1, units="Jy", dim="stokes"))
        myfits.expand_along_axis_from_files("STOKES", fnames_quv)
        
