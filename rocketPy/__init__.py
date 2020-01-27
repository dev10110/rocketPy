"""Top-level package for rocketPy."""

__author__ = """Devansh Ramgopal Agrawal"""
__email__ = 'devanshinspace@gmail.com'
__version__ = '0.1.3'

#print("Hello from rocketPy!")
# print("____________________")

import pint
ureg = pint.UnitRegistry(system='mks')
Q_ = ureg.Quantity

from .materials import *
from .util import *
from .components import *
from .rocket import *
