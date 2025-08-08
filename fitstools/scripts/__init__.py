import click
from scabha.schema_utils import clickify_parameters, paramfile_loader
import glob
import os

thisdir = os.path.dirname(__file__)

source_files = glob.glob(f"{thisdir}/library/*.yaml")
sources = [File(item) for item in source_files]
parserfile = File(f"{thisdir}/{command}.yaml")

stack_config = paramfile_loader(parserfile, sources)[command]


@click.group()
def fitstool():
    pass

@fitstool.command()
@clickify_parameters(config)
def stack(fname)