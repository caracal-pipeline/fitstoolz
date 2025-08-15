from scabha.basetypes import File
from astropy.io import fits
from astropy.table import Table
from astropy import units


def get_beam_table(fname:File):
    fname = File(fname)
    if not fname.EXISTS:
        raise FileNotFoundError(f"Input FITS file '{fname}' does not exist")

    beam_info = {
        "BMAJ": [],
        "BMIN": [],
        "BPA": [],
        "CHAN": [],
        "POL": [], # set to I for now
    }
    
    beam_table = None
    # accept the first beam table in hdulist
    with fits.open(fname) as hdulist:
        header = hdulist[0].header
        for hdu in hdulist:
            if isinstance(hdu, fits.BinTableHDU):
                tab = Table.read(hdu)
                if {"BMAJ", "BMIN", "BPA"}.issubset(tab.colnames):
                    beam_table = tab
                    break
                
    if isinstance(beam_table, Table):
        return beam_table
        
    bunit = getattr(units, header["CUNIT1"])
    chan = 1
    
    while header.get(f"BMAJ{chan}", None) is not None:
        beam_info["BMAJ"].append( header[f"BMAJ{chan}"] * bunit )
        beam_info["BMIN"].append( header[f"BMIN{chan}"] * bunit )
        beam_info["BPA"].append( header[f"BPA{chan}"] * bunit )
        beam_info["CHAN"].append( chan -1 )
        beam_info["POL"].append( 0 )
        chan +=1 
    
    if chan == 1 and header.get(f"BMAJ", None) is not None:
        beam_info["BMAJ"].append( header[f"BMAJ"] * bunit )
        beam_info["BMIN"].append( header[f"BMIN"] * bunit )
        beam_info["BPA"].append( header[f"BPA"] * bunit )
        beam_info["CHAN"].append( 0 )
        beam_info["POL"].append( 0 )
        
        chan = 2
    
    if chan > 1:
        return Table(beam_info)
    else:
        return False
