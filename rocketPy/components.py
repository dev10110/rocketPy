"""Creates the components of a rocket. Base classes and the specific useful components are created.
All components inherit from ``Component``
From here, it splits into
``InternalComponent`` or ``ExternalComponent``
where the difference is used to help the drag model, and plotting.

"""

import numpy as np
from prettytable import PrettyTable
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


from . import Q_, ureg
from . import mach_correction


class Component():
    """Base class for all components

    Args:
        name (str): Name of component. Defaults to 'Component'.
        mass (function or None or Quantity): Mass of component. Defaults to None.
            If ``None``: calls the ``self.estimate_mass()`` method to assign mass.
            If ``function``: calls the function and assigns mass
            If ``Pint.Quantity``: assigns mass directly.
        inertia (function or None or Quantity): Assigns inertia of the rocket, in (I_xx, I_yy, I_zz). Rest of the components are assumed to be 0. Defaults to None.
            If ``None``: calls the ``self.estimate_inertia()`` method to assign inertia.
            If ``function``: calls the function to assign the inertia. The function must return a tuple of (I_xx, I_yy, I_zz)
            If ``tuple of Pint.Quantity``: assigns the inertia direction. The tuple must be (I_xx, I_yy, I_zz) with each being a Pint.Quantity of the right units.

    Examples
        Examples should be written in doctest format, and
        should illustrate how to use the function/class.
        >>>

    Attributes:
        estimate_mass (type): Description of parameter `estimate_mass`.
        I_xx (type): Description of parameter `I_xx`.
        I_yy (type): Description of parameter `I_yy`.
        I_zz (type): Description of parameter `I_zz`.
        estimate_inertia (type): Description of parameter `estimate_inertia`.
        x_ref (type): Description of parameter `x_ref`.
        y_ref (type): Description of parameter `y_ref`.
        z_ref (type): Description of parameter `z_ref`.
        A_ref (type): Description of parameter `A_ref`.
        name
        mass

    """

    def __init__(self, name='Component', mass=None, inertia=None):

        self.name = name

        # mass properties
        if callable(mass):
            self.mass = mass()
        elif mass is None:
            self.mass = self.estimate_mass()
        else:
            self.mass = mass

        self.mass = self.mass.to(ureg.kg)

        # moment of inertias
        # assumed to be 0 initially, but should be replaced
        if callable(inertia):
            self.I_xx, self.I_yy, self.I_zz = inertia()
        elif mass is None:
            self.I_xx, self.I_yy, self.I_zz = self.estimate_inertia()
        else:
            self.I_xx, self.I_yy, self.I_zz = inertia

        # this is the leading reference coordinate of the component, x is along
        # rocket axis, increasing from nose backwards.
        self.x_ref = 0 * ureg.meter
        self.y_ref = 0 * ureg.meter
        self.z_ref = 0 * ureg.meter

        # reference area to be used for drag and lift force calculations.
        self.A_ref = None

    def estimate_mass(self):
        """Generic method to estimate the mass of the component - assume mass is 0.
        This method should be overridden for each component specified"""
        return 0.0 * ureg.kg

    def estimate_inertia(self):
        """Generic method to estimate the mass of the component - assume inertia is 0.
        This method should be overridden for each component specified"""

        return (0.0 * ureg.kg * (ureg.m**2), 0.0 * ureg.kg *
                (ureg.m**2), 0.0 * ureg.kg * (ureg.m**2))

    def set_position(
            self,
            start_of=None,
            end_of=None,
            middle_of=None,
            after=None,
            offset=0 *
            ureg.m):
        """Defines the x_ref of this component relative to other components.

        Args:
            start_of (Component): If not None, aligns self with other. Defaults to None.
            end_of (Component): If not None, aligns end of self with other. Defaults to None.
            middle_of (Component): If not None, aligns midpoints of self and other. Defaults to None.
            after (Component): If not None, aligns start of self to end of other. Defaults to None.
            offset (Pint.Quantity): Follows rules as above, but adds ``offset`` to self.x_ref. Defaults to 0*ureg.m.

        Examples
            Examples should be written in doctest format, and
            should illustrate how to use the function/class.
            >>>

        """

        if start_of is not None:

            self.x_ref = start_of.x_ref + offset

            return

        elif end_of is not None:

            self.x_ref = end_of.x_ref + end_of.length - self.length + offset

            return

        elif middle_of is not None:

            self.x_ref = middle_of.x_ref + middle_of.length / 2 - self.length / 2 + offset

            return

        elif after is not None:

            self.x_ref = after.x_ref + after.length + offset

            return

        else:
            raise ValueError('The location must be specified')

    def __repr__(self):
        return f'{self.name} (type: {self.__class__.__name__})'

    def plot(self, ax=None, rotation=0 * ureg.degree, unit=ureg.m):
        """Plots the component.

        Args:
            ax (type): Description of parameter `ax`. Defaults to None.
            rotation (type): Description of parameter `rotation`. Defaults to 0*ureg.degree.
            unit (type): Description of parameter `unit`. Defaults to ureg.m.

        Returns:
            type: Description of returned object.

        Examples
            Examples should be written in doctest format, and
            should illustrate how to use the function/class.
            >>>

        """

        if ax is None:
            ax = plt.gca()

        for coords in self.plot_coords(rotation=rotation):

            if coords is not None:

                # convert to specified units
                coords2 = np.array([[c[i].m_as(unit)
                                     for i in [0, 1]] for c in coords])

                poly = Polygon(coords2, facecolor='none', edgecolor='k')

                ax.add_patch(poly)

        plt.axis('equal')
        plt.grid(True)
        plt.xlabel(f'x [{str(unit)}]')
        plt.ylabel(f'y or z [{str(unit)}]')
        return ax

    def plot_coords(self):
        return None

    def xcp(self, *args):

        return 0 * ureg.m

    def xcg(self, *args):

        return 0 * ureg.m

    def ycg(self, *args):

        return 0 * ureg.m

    def zcg(self, *args):

        return 0 * ureg.m

    def describe(self):

        print(self)

        x = PrettyTable()

        x.field_names = ["Parameter", "Value (SI)", "Value"]

        for d in self.__dict__:
            try:
                x.add_row([d,
                           f'{self.__dict__[d].to_base_units():.3f~}',
                           f'{self.__dict__[d]:.3f~}'])
            except BaseException:
                x.add_row([d, self.__dict__[d], ""])
        print(x)


class InternalComponent(Component):

    def __init__(self, name='Internal Component', mass=None, inertia=None):

        super().__init__(name=name, mass=mass, inertia=inertia)

    def CN(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):

        return 0.0


class ExternalComponent(Component):

    def __init__(self, name='External Component', mass=None,
                 inertia=None, A_ref=np.pi * (3 * ureg.inch)**2):

        super().__init__(name=name, mass=mass, inertia=inertia)

        self.A_ref = A_ref

    def CNa(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):

        return 0.0 / ureg.rad

    def CN(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):

        return self.CNa(alpha, Re, Mach) * alpha


class Cylinder(InternalComponent):

    def __init__(
            self,
            name='Internal Cylinder',
            mass=None,
            inertia=None,
            diameter=6 *
            ureg.inch,
            length=6 *
            ureg.inch,
            density=None):
        super().__init__(name=name, mass=mass, inertia=inertia)

        self.diameter = diameter
        self.length = length
        self.density = density

    def estimate_mass(self):
        if not self.density.check(ureg.kg / (ureg.m**3)):
            raise ValueError('Density must have the right units')

        V = np.pi * (self.diameter / 2)**2 * self.length

        m = self.density * V

        return m

    def estimate_inertia(self):

        m = self.mass
        r = self.diameter / 2
        h = self.length

        I_xx = 0.5 * m * r**2
        I_yy = (1 / 12) * m * (3 * r**2 + h**2)

        #note, I_yy = I_zz
        return I_xx, I_yy, I_yy


class NoseCone(ExternalComponent):

    def __init__(
            self,
            name='Nose Cone',
            mass=None,
            inertia=None,
            shape='Conical',
            diameter=None,
            length=None,
            fineness=None,
            wall_thickness=2 *
            ureg.mm,
            material=None):
        # note, the plotting currently assumes its a conical nose.

        # use the base diameter for the reference area
        A_ref = np.pi * diameter**2 / 4

        self.shape = shape

        self.x_ref = 0 * ureg.m

        if diameter and length:
            self.diameter = diameter
            self.length = length
        elif diameter and fineness:
            self.diameter = diameter
            self.length = fineness * self.diameter
        elif length and diameter:
            self.length = length
            self.diameter = length / fineness

        self.wall_thickness = wall_thickness
        self.material = material

        super().__init__(name=name, mass=mass, inertia=inertia, A_ref=A_ref)

    def estimate_mass(self):
        """Method to estimate the mass of the nose cone"""

        L = self.length
        R = self.diameter / 2
        t = self.wall_thickness
        rho = self.material.density

        return rho * t * L * np.pi * R

    def estimate_inertia(self):

        L = self.length
        R = self.diameter / 2
        self.wall_thickness
        m = self.mass

        I_xx = m * (R**2 / 2)
        I_yy = m * (L**2 / 18)

        return I_xx, I_yy, I_yy

    def xcg(self):

        return (2 / 3) * self.length

    def CNa(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):
        # eq 25 of [1]

        CNa_incomp = 2 / (1 * ureg.rad)  # per radian

        CNa = CNa_incomp * mach_correction(Mach)

        return CNa

    def xcp(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):

        # at the moment we assume the Xcp doesnt change with Mach

        if self.shape == 'Conical':
            return (2 / 3) * self.length
        elif self.shape == 'Ogive':
            return (5 / 8) * self.length
        elif self.shape == 'Parabolic':
            return (3 / 5) * self.length
        else:
            raise RuntimeError(
                'Please set the nose cone shape to a supported string')

    def plot_coords(self, rotation=0 * ureg.degree):

        l = self.length
        r = self.diameter / 2

        coords = [[0 * ureg.m, 0 * ureg.m], [l, r], [l, -r]]

        coords_shift = [[self.x_ref + c[0], c[1]] for c in coords]

        yield coords_shift


class Transition(ExternalComponent):

    def __init__(
            self,
            name='Transition',
            mass=None,
            inertia=None,
            fore_dia=None,
            aft_dia=None,
            length=None,
            wall_thickness=2 *
            ureg.mm,
            material=None):

        # use the base diameter for the reference area
        self.A_ref = np.pi * fore_dia**2 / 4

        self.fore_dia = fore_dia
        self.aft_dia = aft_dia
        self.d_ref = 2 * (self.A_ref / np.pi)**0.5
        self.length = length
        self.wall_thickness = wall_thickness
        self.material = material

        super().__init__(name=name, mass=mass, inertia=inertia, A_ref=self.A_ref)

    def estimate_mass(self):

        r0 = self.fore_dia / 2
        r1 = self.aft_dia / 2
        L = self.length
        t = self.wall_thickness
        rho = self.material.density

        m = rho * L * (np.pi * (r0 + r1)) * t

        return m

    def estimate_inertia(self):

        L = self.length
        r0 = self.fore_dia / 2
        r1 = self.aft_dia / 2
        m = self.mass

        I_xx = m * (0.5 * (r0**2 + r1**2))

        I_yy = m * (L**2 * (r0**2 + 4 * r0 * r1 + r1**2)) / \
            (18 * (r0 + r1)**2)

        return I_xx, I_yy, I_yy

    def xcg(self):

        r0 = self.fore_dia / 2
        r1 = self.aft_dia / 2
        L = self.length

        return (L / 3) * (r0 + 2 * r1) / (r0 + r1)

    def CNa(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):

        CNa_incomp = 2 * ((self.aft_dia / self.d_ref)**2 -
                          (self.fore_dia / self.d_ref)**2) / (1 * ureg.rad)

        CNa = CNa_incomp * mach_correction(Mach)

        return CNa

    def xcp(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):

        d_fore = self.fore_dia
        d_aft = self.aft_dia

        # eqn 40 of ref[1]
        return (self.length / 3) * ((d_fore + 2 * d_aft) / (d_fore + d_aft))

    def plot_coords(self, rotation=0 * ureg.degree):

        l = self.length
        fore_r = self.fore_dia / 2
        aft_r = self.aft_dia / 2

        coords = [[0 * ureg.m, fore_r], [l, aft_r],
                  [l, -aft_r], [0 * ureg.m, -fore_r]]

        coords_shifted = [[self.x_ref + c[0], c[1]] for c in coords]

        yield coords_shifted


class BodyTube(ExternalComponent):

    def __init__(
            self,
            name='Body Tube',
            mass=None,
            inertia=None,
            diameter=None,
            length=None,
            wall_thickness=2 *
            ureg.mm,
            material=None):

        self.A_ref = np.pi * diameter**2 / 4
        self.d_ref = diameter
        self.diameter = diameter
        self.length = length
        self.wall_thickness = wall_thickness
        self.material = material

        super().__init__(name=name, mass=mass, inertia=inertia, A_ref=self.A_ref)

    def estimate_mass(self):

        r0 = self.diameter / 2
        L = self.length
        t = self.wall_thickness
        rho = self.material.density

        m = 2 * rho * L * np.pi * r0 * t

        return m

    def estimate_inertia(self):

        m = self.mass
        r0 = self.diameter / 2
        L = self.length

        I_xx = m * r0**2

        I_yy = m * (L**2 / 12)

        # note Iyy=Izz
        return I_xx, I_yy, I_yy

    def xcg(self):

        return self.length / 2

    def plot_coords(self, rotation=0 * ureg.degree):

        l = self.length
        r = self.diameter / 2

        coords = [[0 * ureg.m, r], [l, r], [l, -r], [0 * ureg.m, -r]]

        coords_shifted = [[self.x_ref + c[0], c[1]] for c in coords]

        yield coords_shifted

    def CNa(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3, K=1.1):

        planform_area = self.diameter * self.length

        CNa_incomp = (K * planform_area / self.A_ref) * \
            alpha / ((1 * ureg.rad)**2)

        CNa = CNa_incomp * mach_correction(Mach)

        return CNa

    def xcp(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):

        return self.length / 2


class FinSet(ExternalComponent):

    def __init__(
            self,
            name='Fins',
            mass=None,
            inertia=None,
            n=None,
            span=None,
            root_chord=None,
            tip_chord=None,
            mid_sweep=None,
            tube_dia=None,
            thickness=None,
            material=None):

        # midsweep is the sweep angle at the mid-chord locations

        self.A_ref = np.pi * tube_dia**2 / 4

        self.d_ref = 2 * (self.A_ref / np.pi)**0.5

        self.n = n
        self.span = span
        self.root_chord = root_chord
        self.tip_chord = tip_chord
        self.mid_sweep = mid_sweep
        self.mid_chord_span = self.span / np.cos(self.mid_sweep)

        self.tube_dia = tube_dia
        self.length = 0 * ureg.m  # used in calculating the overall length of the rocket
        self.thickness = thickness

        self.exposed_area = 0.5 * \
            (self.root_chord + self.tip_chord) * self.span  # per fin
        self.planform_area = self.exposed_area + 0.5 * \
            self.tube_dia * self.root_chord  # per fin

        self.material = material

        super().__init__(name=name, mass=mass, inertia=inertia, A_ref=self.A_ref)

    def estimate_mass(self):

        n = self.n
        A = self.exposed_area
        t = self.thickness

        rho = self.material.density

        m = n * A * t * rho

        return m

    def xcg(self):

        # using mathematica

        cr = self.root_chord
        ct = self.tip_chord
        s = self.span
        r0 = self.tube_dia / 2
        gamma = self.leading_sweep()

        xcg = (cr**2 + cr * ct + ct**2 + (3 * (cr + ct) * r0 +
                                          (cr + 2 * ct) * s) * np.tan(gamma)) / (3 * (cr + ct))

        return xcg

    def estimate_inertia(self):

        r0 = self.tube_dia / 2
        cr = self.root_chord
        ct = self.tube_dia
        s = self.span
        m = self.mass
        gamma = self.leading_sweep()

        #I_xx_m is I_xx/m
        # this is for a single fin, but since I'm scaling it by the number of
        # fins, its all good.
        I_xx_m = r0**2 + (2 * (cr + 2 * ct) * r0 * s) / (3 *
                                                         (cr + ct)) + ((cr + 3 * ct) * s**2) / (6 * (cr + ct))
        I_xx = m * I_xx_m

        I_yy_m = (1 / (18 * (cr + ct)**2)) * (cr**4 + 2 * cr**3 * ct + 2 * cr * ct**3 + ct **
                                              4 - (cr**2 + 4 * cr * ct + ct**2) * s * np.tan(gamma) * (cr - ct - s * np.tan(gamma)))
        I_yy = m * I_yy_m

        return I_xx, I_yy, I_yy

    def leading_sweep(self):
        """Return the leading edge sweep of the fins"""

        lr = self.root_chord
        lt = self.tip_chord
        ls = self.span
        sweep = self.mid_sweep

        tip_le = lr / 2 + ls * np.tan(sweep) - lt / 2  # tip leading edge

        leading_sweep = np.arctan(tip_le / ls)

        return leading_sweep

    def CNa(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):

        if self.n <= 4:
            body_influence = 1 + (self.tube_dia / 2) / \
                (self.span + self.tube_dia / 2)
        elif self.n > 4:
            body_influence = 1 + 0.5 * \
                (self.tube_dia / 2) / (self.span + self.tube_dia / 2)

        CNa_incomp = (body_influence * (4 * self.n * (self.span / self.d_ref)**2) / (1 + (1 + (
            self.mid_chord_span / (self.root_chord + self.tip_chord))**2)**0.5)) / (1 * ureg.rad)

        CNa = CNa_incomp * mach_correction(Mach)

        return CNa

    def xcp(self, alpha=0 * ureg.rad, Re=1e6, Mach=0.3):

        lm = self.mid_chord_span
        lr = self.root_chord
        lt = self.tip_chord

        a = lm * (lr + lt) / (3 * (lr + lt))
        b = (1 / 6) * (lr + lt - (lr * lt) / (lr + lt))

        return a + b

    def plot_coords(self, rotation=0 * ureg.rad):

        r = self.tube_dia / 2
        lr = self.root_chord
        lt = self.tip_chord
        self.mid_chord_span
        ls = self.span
        sweep = self.mid_sweep

        tip_le = lr / 2 + ls * np.tan(sweep) - lt / 2  # tip leading edge

        x_ref = self.x_ref

        coords_single = [[x_ref, r, 0 *
                          ureg.m], [x_ref +
                                    tip_le, r +
                                    ls, 0 *
                                    ureg.m], [x_ref +
                                              tip_le +
                                              lt, r +
                                              ls, 0 *
                                              ureg.m], [x_ref +
                                                        lr, r, 0 *
                                                        ureg.m]]

        coords_single_m = np.array(
            [[cc.m_as(ureg.m) for cc in c] for c in coords_single])

        th_set = np.linspace(0, 2 * np.pi, self.n, endpoint=False)

        for th in th_set:

            th2 = th + rotation.m_as(ureg.rad)

            # we only need the x and y components, thus a 2x3 rotation matrix
            # is created.
            R2 = np.array([[1, 0, 0], [0, np.cos(th2), -np.sin(th2)]])

            coords_rotated = (R2 @ coords_single_m.T).T

            coords_units = [[c[0] * ureg.m, c[1] * ureg.m]
                            for c in coords_rotated]

            yield coords_units

    # references:
    # [1]: Simon Box, 2009, Estimating the dynamic and aerodynamic paramters of passively controlled high power rockets for flight simulaton
