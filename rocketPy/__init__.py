"""Top-level package for rocketPy."""

__author__ = """Devansh Ramgopal Agrawal"""
__email__ = 'devanshinspace@gmail.com'
__version__ = '0.1.2'

print("Hello from rocketPy!")

import numpy as np
import pint
ureg = pint.UnitRegistry(system='mks')
Q_ = ureg.Quantity

from .materials import *
from .util import *
