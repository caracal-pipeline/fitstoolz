import logging
import os
import uuid

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from fitstoolz import LOGGER
from fitstoolz.reader import FitsData

log = logging.getLogger(LOGGER)

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class InitTests:
    def __init__(self):
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

    def __del__(self):
        """Clean up after each test method runs."""
        # Remove any temporary files created
        for file in self.test_files:
            if os.path.exists(file):
                os.remove(file)


@pytest.fixture
def params():
    return InitTests()


class TestExpandAlongDim:
    def test_fits_predicting_all_stokes_linear_basis(self, params):
        """
        Test visibility prediction from FITS images of all Stokes parameters, ncorr = 4
        Validates:
            - Output shape of visibilities
            - XX = I + Q
            - XY = U + iV
            - YX = U - iV
            - YY = I - Q
        """
        params.ncorr = 4
        # the numbers below are unphysical—they are just for testing the computation
        stokes_params = [("I", 1.0), ("Q", 2.0), ("U", 3.0), ("V", 4.0)]

        filenames = []
        for stokes in stokes_params:
            wcs = WCS(naxis=3)
            wcs.wcs.ctype = ["RA---SIN", "DEC--SIN", "FREQ"]
            wcs.wcs.cdelt = np.array([-params.cell_size / 3600, params.cell_size / 3600, params.dfreq])
            wcs.wcs.crpix = [params.img_size // 2, params.img_size // 2, 1]
            wcs.wcs.crval = [np.rad2deg(params.ra0), np.rad2deg(params.dec0), params.start_freq]

            # make header
            header = wcs.to_header()
            header["BUNIT"] = "Jy"

            # make image
            image = np.zeros((params.nchan, params.img_size, params.img_size))
            image[:, params.img_size // 2, params.img_size // 2] = stokes[1]  # put a point source at the center
            # write to FITS file
            hdu = fits.PrimaryHDU(image, header=header)
            hdul = fits.HDUList([hdu])
            test_filename = f"{TESTDIR}/test_{uuid.uuid4()}_{stokes[0]}.fits"
            params.test_files.append(test_filename)
            hdul.writeto(test_filename, overwrite=True)
            hdul.close()

            filenames.append(test_filename)

        fname = filenames[0]
        fnames_quv = filenames[1:]
        myfits = FitsData(fname)

        myfits.add_axis("STOKES", idx=4, crval=1, cdelt=1, crpix=0, cunit="Jy")
        myfits.expand_along_axis_from_files("STOKES", fnames_quv)
