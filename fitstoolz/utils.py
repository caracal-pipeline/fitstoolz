from scabha.basetypes import File
from astropy.io import fits
from astropy.table import Table
from astropy import units


def get_beam_table(fname:File, hdu_index=1, freqs=None):
    fname = File(fname)
    if not fname.EXISTS:
        raise FileNotFoundError(f"Input FITS file '{fname}' does not exist")
    
    try:
        return Table.read(fname, hdu=hdu_index)
    except ValueError:
        pass
    
    with fits.open(fname) as hdulist:
        hdulist = fits.open(fname)
        header = hdulist[0].header
        
    beam_info = {
        "BMAJ": [],
        "BMIN": [],
        "BPA": [],
        "CHAN": [],
        "POL": [], # set to I for now
    }
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

    
        
        
        
        