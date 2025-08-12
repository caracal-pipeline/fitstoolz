from scabha.schema_utils import clickify_parameters, paramfile_loader
from scabha.basetypes import File
from . import thisdir, sources, cli
from fitstoolz.reader import FitsData
from omegaconf import OmegaConf
from fitstoolz import LOG

script = "stack"
parserfile = File(f"{thisdir}/{script}.yaml")
stack_config = paramfile_loader(parserfile, sources)[script]
@clickify_parameters(stack_config)
@cli.command("stack")
def stack(**kwargs):
    opts = OmegaConf.create(kwargs)
    
    fname0 = opts.fname[0]
    fnames = opts.fname[1:]
    myfits = FitsData(fname=fname0, memmap=opts.memmap)
    
    myfits.expand_along_axis_from_files(opts.axis, fnames)
    coord_names = myfits.coord_names
    spectral = myfits.spectral_coord
    chunks = {
        "RA": opts.ra_chunks or myfits.coords["RA"].size,
        "DEC": opts.dec_chunks or myfits.coords["DEC"].size,
        spectral: opts.spectral_chunks or myfits.nchan
    }
        
    myfits.write_to_fits(opts.stacked_fits, coord_names=coord_names, chunks=chunks)
    LOG.info(f"Wrote stacked file to: {opts.stacked_fits}")
