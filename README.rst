==========
fitstools
==========

Python libraries for interfacing with FITS images in an intuitive way

Features
--------

* Image and coordinate data indexing is consistent
* Simple API for adding, transposing and expanding axes (or dimensions)
* `Xarray support <https://docs.xarray.dev/en/stable/index.html>`_
* `Zarr support <https://zarr.readthedocs.io/en/stable/index.html>`_ 

Example Usage
-------------

.. code-block:: ipython

    In[1]: from fitstools.reader import FitsData

    In[2]: myfits = FitsData("example-image.fits")
           myfits.coord_names

    Out[2]: ['STOKES', 'FREQ', 'DEC', 'RA']

    In[3]: myfits.dshape
    Out[3]: (1, 504, 100, 100) # these dimensions match the labels above

    In[4]: myfits.coords

    Out[4]: 
    Coordinates:
    STOKES   (stokes) int32 4B dask.array<chunksize=(1,), meta=np.ndarray>
    FREQ     (spectral) float64 4kB 8.803e+08 8.804e+08 ... 9.328e+08 9.329e+08
    RA       (celestial.ra) float64 800B 53.16 53.16 53.16 ... 53.1 53.1 53.1
    DEC      (celesstial.dec) float64 800B -28.16 -28.16 ... -28.11 -28.11
