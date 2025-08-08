from astropy.io import fits
import dask.array as da
import xarray as xr
import numpy as np
from astropy.wcs import WCS
from scabha.basetypes import File
from typing import List, Dict, Any
from astropy.table import Table
from astropy.coordinates import SpectralCoord
from astropy import units

class FitsData:
    def __init__(self, fname: str, memmap: bool = True):
        self.fname = File(fname)
        if not self.fname.EXISTS:
            raise FileNotFoundError(f"Input FITS file '{fname}' does not exist")
            
        self.hdulist = fits.open(self.fname, memmap=memmap)
        self.phdu = self.hdulist[0]
        self.header = self.phdu.header
        self.wcs = WCS(self.header)
        self.dim_info = self.wcs.get_axis_types()[::-1]
        self.coord_names = self.wcs.axis_type_names[::-1]
        self.coords = xr.Coordinates()
        self.open_arrays = []
        self.data = da.asarray(self.phdu.data)
        self.data_units = self.header.get("BUNIT", "jy").lower()
        
        if self.dshape != self.wcs.array_shape:
            raise RuntimeError("Input FITS file WCS information does not match Image data")
    
        
        self.register_dimensions()

    def coord_index(self, name:str) -> int:
        """
        Returns index of given axis

        Args:
            name (str): Coordinate/axis name

        Returns:
        int : coordinate index
        """
        return self.coord_names.index(name)

    @property 
    def ndim(self):
        return self.data.ndim
    
    @property
    def dims(self):
        return [self.coords[name].dim for name in self.coord_names]
        

    @property
    def dshape(self):
        return self.data.shape

    @property
    def data(self):
        return self._data 
    
    @data.setter
    def data(self, value):
        self._data = value
    
    def set_coord_attrs(self, name:str, dim:str):
        """ 
        Add FITS pixel meta data to this instances coords attribute (xarray.Coordinates)

        Args:
            name (str): Name (or label) of coordinate to set
        """
        idx = self.coord_index(name)
        self.coords[name].attrs = {
            "name": name,
            "pixel_size": self.pixel_scales[name],
            "dim": dim,
            "ref_pixel": int(self.header[f"CRPIX{self.ndim - idx }"]) - 1, # FITS indexing is 1-based
            "units": self.wcs.world_axis_units[::-1][idx],
            "size": self.dshape[idx],
        }
    
    def register_dimensions(self):
        """
        Register (or set) FITS data coordinate and dimension data from the FITS header
        (including WCS information)
        
        Args:
            set_dims (List[str]): FITS coordinates (celestial, spectral, etc.) that must be set according to header WCS information. For example, if 'celestial' is part of this list, then the RA and DEC coordinate grid will be populated using the header WCS instead a dummy array. Avoid inlcuding dimesions for which you will need the coordinate grids for.
        """

        self.pixel_scales = {}
        pixel_scales = self.wcs.pixel_scale_matrix.diagonal()[::-1]
        
        names = self.coord_names 
        for idx,coord in enumerate(names):
            self.pixel_scales[coord] = pixel_scales[idx]
            diminfo = self.dim_info[idx]
            dim = diminfo["coordinate_type"]
            dtype = diminfo["group"]
            
            if dim == "celestial":
                continue
            elif dim == "spectral":    
                self.set_spectral_dimension(coord)
                continue
            elif dim == "stokes":
                self.set_stokes_dimensions(coord)
                continue
                
            dimsize = self.dshape[idx]
            try:
                dimgrid = getattr(self.wcs,
                                dim).array_index_to_world_values(da.arange(dimsize))
            except AttributeError:
                dimgrid = da.empty(dimsize, dtype=dtype)
                
            self.coords[coord] = (dim,), dimgrid
            self.set_coord_attrs(coord, dim)
            
        self.set_celestial_dimensions()
        
    def set_stokes_dimensions(self, coord_name:str):
        """
        Args:
            coord_name (str): _description_
        """
        
        idx = self.coord_index(coord_name)
        crval = int(self.header.get(f"CRVAL{self.ndim-idx}", 1))
        crpix = int(self.header.get(f"CRPIX{self.ndim-idx}", 1))
        cdelt = int(self.header.get(f"CDELT{self.ndim-idx}", 1))
        naxis = int(self.header.get(f"NAXIS{self.ndim-idx}"))
        size = cdelt*naxis
        grid = da.arange(0,size,cdelt) 
        shift = grid[crpix-1] + crval
        grid += shift
        
        self.coords[coord_name] = ("stokes",), grid
        self.set_coord_attrs(coord_name, "stokes")
        
        
    def set_spectral_dimension(self, coord_name, rest_freq=None):
        idx = self.coord_index(coord_name)
        dimsize = self.dshape[idx]
        
        dimgrid = self.wcs.spectral.array_index_to_world_values(da.arange(dimsize))                
        self.coords[coord_name] = ("spectral",), dimgrid
        self.set_coord_attrs(coord_name, "spectral")
        self.nchan = dimsize
        self.spectral_coord = coord_name
        self.spectral_refpix = self.coords[coord_name].ref_pixel
        self.spectral_units = self.wcs.spectral.world_axis_units[0]
        self.spectral_restfreq = rest_freq or self.header.get("RESTFRQ", None)
        
    def set_celestial_dimensions(self):
        """
        set celestial coordinates data using the FITS WCS infomation

        Args:
            empty (bool, optiona): Set the coordinate grids as an empty array. Defaults to True.
        """
        
        for idx, diminfo in enumerate(self.dim_info):
            dim = diminfo["coordinate_type"]
            dim_number = diminfo["number"]
            if dim == "celestial":
                if dim_number == 0:
                    ra_idx = idx
                elif dim_number == 1:
                    dec_idx = idx
                else:
                    raise ValueError(f"Unkown celestial dimension in WCS: {dim_number}")
                    
        ra_dim = self.wcs.axis_type_names[::-1 ][ra_idx]
        dec_dim = self.wcs.axis_type_names[::-1 ][dec_idx]
        ra_dimsize = self.dshape[ra_idx]
        dec_dimsize = self.dshape[dec_idx]
            
        grid = self.wcs.celestial.array_index_to_world_values(da.arange(ra_dimsize),
                                        da.arange(dec_dimsize))
        
        self.coords[ra_dim] = ("celestial.ra",), grid[0]
        self.set_coord_attrs(ra_dim, "celestial.ra")
        
        self.coords[dec_dim] = ("celestial.dec",), grid[1]
        self.set_coord_attrs(dec_dim, "celestial.dec")

        
    def get_freq_from_vrad(self, rest_freq_Hz=None):
        """
        Convert radio velocity coordinates to frequencies

        Returns:
            astropy.SpectralCoord: Astropy SpectralCoord instance
        """
        
        rest_freq_Hz = rest_freq_Hz or self.spectral_restfreq
        return SpectralCoord(self.coords["VRAD"], 
                unit=units.meter/units.second).to(units.Hz, 
                                            doppler_rest = rest_freq_Hz*units.Hz,
                                            doppler_convention="radio").value
        
    def get_freq_from_vopt(self, rest_freq_Hz=None):
        """
        Convert radio velocity coordinates to frequencies

        Returns:
            astropy.SpectralCoord: Astropy SpectralCoord instance
        """
        
        rest_freq_Hz = rest_freq_Hz or self.spectral_restfreq
        return SpectralCoord(self.coords["VOPT"],
                unit=units.meter/units.second).to(units.Hz, 
                                            doppler_rest = rest_freq_Hz*units.Hz,
                                            doppler_convention="optical").value

    
    def add_axis(self, name:str, idx:int, coord_type:str, axis_grid:np.ndarray, attrs:Dict):
        """ Add a new axis to FITS data

        Args:
            name (str): Name of new axis. e.g, STOKES, RA, DEC
            idx (int): index where new axis must be added
            coord_type (str): Axis type (e.g, stokes, spectral)
            axis_grid (np.ndarray): coordinate grid of new axis
            attrs (Dict): Pixel meta data. Example: dict(name='STOKES', pixel_size=1, ref_pixel=0, size=1, units=None, dim='stokes')

        Raises:
            RuntimeError: Dimensions not matching after axis was added
        """
        
        slc = [slice(None)] * (self.ndim + 1)
        slc[idx] = da.newaxis
        self.data = self.data[tuple(slc)]
        self.coords[name] = (coord_type,), axis_grid
        self.coords[name].attrs = attrs
        self.coord_names.insert(idx, name)
        if len(self.coord_names) != self.ndim:
            raise RuntimeError(f"New axis '{name}' could not added")
    
    def expand_along_axis(self, name:str , data:np.ndarray):
        """
        Expand data along the given dimension. 

        Args:
            name (str): Name of expansion coordinate
            data (np.ndarray): Data slice to add along the axis. Number of dimensions 
                must be equal or one less than the current data
        """
        
        idx = self.coord_index(name)
        if len(data.shape) == self.ndim - 1:
            slc = [slice(None)] * self.ndim
            slc[idx] = da.newaxis
            data = data[tuple(slc)]
        self.data = da.concatenate((self.data, data), axis=idx)
        self.coords[name].attrs.update({"size": self.dshape[idx]})
    
    def expand_along_axis_from_files(self, name, files:List[File]):
        idx = self.coord_index(name)
        
        for fname in files:
            with fits.open(fname, memmap=True) as hdul:
                slc = [slice(None)] * self.ndim
                if self.ndim != hdul[0].data.ndim:
                    slc[idx] = da.newaxis
                slc = tuple(slc)
                data = da.asarray(hdul[0].data[slc])
            self.data = da.concatenate((self.data, data), axis=idx)
            self.coords[name].attrs.update({"size": self.dshape[idx]})

    def register_beam_info(self):
        """
        Get FITS beam information and assign it to a self.beam_info attribute:
        self.beam_info = {\
            "bmaj": <array>,\
            "bmin": <array>,\
            "bpa": <array>, \
        }
        """
        header = self.header
        beam_table = None
        
        try:
            beam_table = Table.read(self.fname)
        except ValueError:
            pass
        
        beam_info = {}
        beam_info["bmaj"] = np.zeros(self.nchan)
        beam_info["bmin"] = np.zeros(self.nchan)
        beam_info["bpa"] = np.zeros(self.nchan)
    
        if beam_table:
            bunit = beam_table["BMAJ"].unit
            for chan in range(self.nchan):
                beam = beam_table[chan]
                beam_info["bmaj"][chan] = beam["BMAJ"]
                beam_info["bmin"][chan] = beam["BMIN"]
                beam_info["bpa"][chan] = beam["BPA"]
                
        elif header.get(f"BMAJ1", False):
            bunit = self.coords["RA"].attrs["units"]
            for chan in range(self.nchan):
                beam_info["bmaj"][chan] = header[f"BMAJ{chan+1}"]
                beam_info["bmin"][chan] = header[f"BMIN{chan+1}"]
                beam_info["bpa"][chan] = header[f"BPA{chan+1}"]
        elif header.get("BMAJ", False):
            bunit = self.coords["RA"].attrs["units"]
            if self.spectral_coord == "VRAD":
                freqs = self.get_freq_from_vrad()
            elif self.spectral_coord == "VOPT":
                freqs = self.get_freq_from_vopt()
            else:
                freqs = self.coords["FREQ"].data
            
            for chan in range(self.nchan):
                scale_factor = freqs[self.spectral_refpix]/freqs[chan]
                beam_info["bmaj"][chan] = header[f"BMAJ"] * scale_factor
                beam_info["bmin"][chan] = header[f"BMIN"] * scale_factor
                beam_info["bpa"][chan] = header[f"BPA"]
        else:
            #TODO(mika): Print warning
            return
        
        # convert to radians
        if isinstance(bunit, str):
            bunit = getattr(units, bunit)
        
        beam_info["bmaj"] = (beam_info["bmaj"]*bunit).to(units.rad).value
        beam_info["bmin"] = (beam_info["bmin"]*bunit).to(units.rad).value
        beam_info["bpa"] = (beam_info["bpa"]*bunit).to(units.rad).value
        
        self.beam_info = beam_info
        return 

    @property
    def data(self):
        return self._data 
    
    @data.setter
    def data(self, value):
        self._data = value
    
    def get_data(self, data_slice=None) -> np.ndarray:
        if data_slice:
            data = self.phdu.section[tuple(data_slice)]
        else:
            data = self.phdu.data
        self.open_arrays.append(data)
        
        return self.open_arrays[-1]
        
    def get_xds(self, data_slice=[],
                transpose=[], chunks=dict(RA=64,DEC=64), **kwargs):
        
        dim_chunks = {}
        # swap coords for dims in chunks
        for key,val in chunks.items():
            if key in self.coord_names:
                key = self.coords[key].attrs["dim"]
            dim_chunks[key] = val
            
        if len(transpose) == 0:
            transpose = self.dims
            
        # swap coords for dims in transpose
        dim_transpose = []
        for key in transpose:
            if key in self.coord_names:
                key = self.coords[key].attrs["dim"]
            dim_transpose.append(key)
        
        if len(data_slice) == 0:
            data_slice = [slice(None)]*self.ndim
            
        data = da.asarray(self.data[tuple(data_slice)])
        
        xds = xr.DataArray(data,
            coords = self.coords,
            attrs = {
                "header": dict(self.header.items()),
            },
            **kwargs,
        ).transpose(*dim_transpose)
        
        return xds.chunk(dim_chunks)
    
    def write_to_fits(self, fname:File, coord_names:list[str], data_slice:List[Any]=[],
                    chunks:Dict={}): 
        """ Write FitsData object into a FITS file

        Args:
            fname (File): Name of FITS image to write
            coord_names (list[str]): Coordinates to include in the FITS image. The ordering is the Python convention.
                    For example, to get a FITS image with RA -> NAXIS1, DEC -> NAXIS2, FREQ -> NAXIS3, STOKES -> NAXIS4,
                    you need to give coord_names=['STOKES', 'FREQ', 'DEC', 'RA'].
            data_slice (slice): Tuple of sclices
            chunks (Dict|Mapping): How to chunk data when writing to disk 
        """
        
        for i, coord in enumerate(coord_names):
            idx = self.ndim - i
            
            cdelt = self.coords[coord].pixel_size
            crpix = self.coords[coord].ref_pixel
            cunit = self.coords[coord].units
            crval = da.compute(self.coords[coord].data[crpix])
            
            header = fits.Header(self.header)
            opts = {
                f"CDELT{idx}": cdelt,
                f"CRPIX{idx}": crpix + 1,
                f"CUNIT{idx}": cunit,  
                f"CRVAL{idx}": crval,
            }
            header.update(opts)
        xds = self.get_xds(data_slice, coord_names, chunks)
        phdu = fits.PrimaryHDU(xds.data, 
                            header=header)
        
        phdu.writeto(fname, overwrite=True)
    
    def __close__(self):
        self.hdulist.close()
        for data in self.open_arrays:
            del data
    
    def  close(self):
        self.__close__()
            
    def __exit__(self):
        self.__close__()
