"""Defines all the materials and their material properties.
Base class Material allows for users to define a material, while some specific commonly used materials are predefined to help speed up the design process.
Users should check that the right material properties are assumed for their parts."""

from . import Q_, ureg


class Material():
    """Base class for defining material properties

    Args:
        name (str): Name of material.

    Examples
        Examples should be written in doctest format, and
        should illustrate how to use the function/class.
        >>>

    Attributes:
        density (Pint.Quantity): Material Density.
        name (str): name

    """

    def __init__(self, name):

        self.name = name

        self.density = 0 * ureg.km / (ureg.m**3)

    def __repr__(self):

        return f'{self.name}: (Material))'

    def describe(self):
        for d in self.__dict__:
            print(f'{d:20s}: {str(self.__dict__[d]):20s}')


class Aluminium(Material):
    """Defines a basic aluminium.

    Args:
        name (str): Description of parameter `name`. Defaults to 'Al-6061-T6'.

    Examples
        Examples should be written in doctest format, and
        should illustrate how to use the function/class.
        >>>

    Attributes:
        density (Pint.Quantity): Description of parameter `density`.
        tensile_modulus (Pint.Quantity): Description of parameter `tensile_modulus`.
        tensile_strength (Pint.Quantity): Description of parameter `tensile_strength`.
        max_temp (Pint.Quantity): Description of parameter `max_temp`.

    """

    def __init__(self, name='Al-6061-T6'):

        super().__init__(name=name)

        self.density = 2.7 * ureg.g / (ureg.cm**3)
        self.tensile_modulus = 69 * ureg.GPa
        self.tensile_strength = 270 * ureg.MPa
        self.max_temp = 420 * ureg.degK


class PLA(Material):

    def __init__(self, name='PLA'):

        super().__init__(name=name)

        self.density = 1.05 * ureg.g / (ureg.cm**3)


class Phenolic(Material):

    def __init__(self, name='Phenolic'):

        super().__init__(name=name)

        self.density = 0.95 * ureg.g / (ureg.cm**3)


class Acrylic(Material):

    def __init__(self, name='Acrylic'):

        super().__init__(name=name)

        self.density = 1.19 * ureg.g / (ureg.cm**3)


class Plywood(Material):

    def __init__(self, name='Plywood'):

        super().__init__(name=name)

        self.density = 0.63 * ureg.g / (ureg.cm**3)


class Polycarbonate(Material):

    def __init__(self, name='Polycarbonate'):

        super().__init__(name=name)

        self.density = 1.2 * ureg.g / (ureg.cm**3)
