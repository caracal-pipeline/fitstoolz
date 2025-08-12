import click
from scabha.basetypes import File
import glob
import os


thisdir = os.path.dirname(__file__)

source_files = glob.glob(f"{thisdir}/cabs/*.yaml")
sources = [File(item) for item in source_files]

@click.group()
def cli():
    pass

def add_commands():
    # Importing the commands in a function to avoid a circular import error
    from .stack import stack
    
add_commands()

