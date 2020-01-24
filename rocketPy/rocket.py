"""Module to describe the rocket"""

from collections.abc import Iterable
import numpy as np
from prettytable import PrettyTable
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from . import ureg, Q_
from .components import BodyTube
from .util import mach_correction, si

class Rocket():

    def __init__(self, name='Rocket'):

        self.name = name


        self.components = []

        self.nose_cone = None
        self.body_tube = None
        self.fins = None
        self.boat_tail = None

        self.A_ref = 0*ureg.inch**2

        return

    def __repr__(self):
        return f'{self.name}'

    def add(self, component):

        if isinstance(component, Iterable):
            component_unique = [comp for comp in component if comp not in self.components]
            self.components.extend(component_unique)
        elif component not in self.components:
            self.components.append(component)

        return None

    def set_boat_tail(self, component):

        """Set the boat tail component"""

        if component in self.components:
            # find the boat tail in the list of components
            boat_tail = [comp for comp in self.components if component == comp][0]
        else:
            boat_tail = component
            self.add(component)

        self.boat_tail = boat_tail

        return None

    def set_fins(self, component):

        """Set the fin set of the rocket"""

        if component in self.components:
            # find the fins in the list of components
            fins = [comp for comp in self.components if component == comp][0]
        else:
            fins = component
            self.add(fins)

        self.fins = fins
        return None

    def set_nose_cone(self, component):

        """Set the nose cone of the rocket"""

        if component in self.components:

            # find the nc in the list of components
            nc = [comp for comp in self.components if component == comp][0]
        else:
            nc = component
            self.add(nc)

        self.nose_cone = nc

        # define the A_ref here.
        self.A_ref = nc.A_ref

        return None

    def set_body_tube(self, component):

        """Set the body tube of the rocket"""

        if component in self.components:

            # find the body tube in the list of components
            bt = [comp for comp in self.components if component == comp][0]
        else:
            bt = component
            self.add(bt)

        self.body_tube = bt

        return None


    def plot(self, ax=None, unit=ureg.m, rotation=0*ureg.degree, plot_component_cp=True,plot_component_cg=True, alpha=0*ureg.degree, Re=1e6, Mach=0.3):

        if not ax:
            ax = plt.gca()

        for comp in self.components:
            comp.plot(ax, rotation=rotation, unit=unit)


        if plot_component_cp:
            self.plot_cp(ax, unit=unit, alpha=alpha, Re=Re, Mach=Mach)

        if plot_component_cg:
            self.plot_cg(ax, unit=unit)

        plt.axis('equal')
        plt.grid(True)
        plt.xlabel(f'x [{str(unit)}]')
        plt.ylabel(f'y or z [{str(unit)}]')

        return ax

    def plot_cp(self, ax=None, unit=ureg.m, alpha=0*ureg.degree, Re=1e6, Mach=0.3):

        if not ax:
            ax = plt.gca()

        # plot for components
        for comp in self.components:
            ax.plot((comp.x_ref+comp.xcp(alpha, Re, Mach)).m_as(unit), 0, 'kx')

        # plot for rocket
        ax.plot(self.xcp().m_as(unit), 0, 'rx')

        return ax

    def plot_cg(self, ax=None, unit=ureg.m):

        if not ax:
            ax = plt.gca()

        # plot for components
        for comp in self.components:
            ax.plot((comp.x_ref+comp.xcg()).m_as(unit), 0, 'ko')

        # plot for rocket
        ax.plot(self.xcg().m_as(unit), 0, 'ro')

        return ax


    def length(self):
        return sum(comp.length for comp in self.components)

    def CNa(self, alpha=0*ureg.rad, Re=1e6, Mach=0.3):

        CNa = sum(comp.CNa(alpha, Re, Mach) for comp in self.components)

        return CNa

    def CN(self, alpha=0*ureg.rad, Re=1e6, Mach=0.3):

        return self.CNa(alpha, Re, Mach) * alpha

    def xcp(self, alpha=0*ureg.rad, Re=1e6, Mach=0.3):

        xcp = sum(comp.CNa(alpha, Re, Mach) * (comp.A_ref/self.A_ref) * (comp.x_ref + comp.xcp(alpha, Re, Mach)) for comp in self.components) / self.CNa(alpha, Re, Mach)

        return xcp

    def mass(self):

        m = sum(comp.mass for comp in self.components)

        return m

    def xcg(self, mass=None):

        # TODO (high): update to accomodate changing mass of rocket.

        xcg = sum(comp.mass * (comp.x_ref + comp.xcg()) for comp in self.components) / self.mass()

        return xcg

    def inertia_matrix(self, mass=None):

        # TODO (high): change this to return the inertia at the given mass total

        I_xx, I_yy, I_zz = si(self.inertia())

        return np.diagflat([I_xx, I_yy, I_zz])

    def inertia(self):

        return self.inertia_xx(), self.inertia_yy(), self.inertia_zz()

    def inertia_xx(self):

        return sum(comp.I_xx + comp.mass*(comp.x_ref + comp.xcg())**2 for comp in self.components)

    def inertia_yy(self):

        return sum(comp.I_yy + comp.mass*(comp.y_ref + comp.ycg())**2 for comp in self.components)

    def inertia_zz(self):

        return sum(comp.I_zz + comp.mass*(comp.z_ref + comp.zcg())**2 for comp in self.components)

    def CD(self, alpha=0*ureg.rad, Re=1e6, Mach=0.3):
        """Calculate the drag force at some angle of attack, including compressibility"""

        CD0 = self.CD0(Re)
        CD_body_alpha = self.CD_body_alpha(alpha)
        if self.fins is not None:
            CD_fin_alpha = self.CD_fin_alpha(alpha)
        else:
            CD_fin_alpha = 0.0

        CD = CD0 + CD_body_alpha + CD_fin_alpha

        # perform a mach number correction
        CD = CD*mach_correction(Mach)

        return CD

    def CA(self, alpha=0*ureg.rad, Re=1e6, Mach=0.3):
        """Compute the axial drag force from the normal force and the axial force"""

        # get 0 mach CD:
        CD = self.CD(alpha, Re, Mach)
        CN = self.CN(alpha, Re, Mach)

        # eqn. 54 of Box 2009
        CA = (CD*np.cos(alpha) - 0.5*CN*np.sin(2*alpha))/(1-np.sin(alpha)**2)

        return CA


    def CD_body_alpha(self, alpha=0*ureg.rad):

        def delta(alpha):
            # in radians (by fitting to Fig. 4 of Box 2009)
            # collected raw data: (alpha (deg), delta) = {{4, 0.780021417}, {6, 0.857352918}, {8, 0.92048524}, {10, 0.940041875}, {12, 0.960026851}, {16, 0.975050746}, {18, 0.980015024}}
            return 1 - 0.518535*np.exp(-0.00378764*alpha)

        def eta(alpha):
            # in radians
            return 0.259348*alpha**0.153029

        l_TR = self.length()
        l_n = self.nose_cone.length
        alpha = alpha.m_as(ureg.rad) # ensure its in radians

        # maximum body diameter of rocket
        d_b = max([comp.diameter for comp in self.components if type(comp) is BodyTube])


        CD_body_alpha = 2*delta(alpha)*alpha**2 + (3.6 * eta(alpha) * ((1.36 * l_TR - 0.55* l_n))/(np.pi * d_b)) * alpha**3

        return CD_body_alpha

    def CD_fin_alpha(self, alpha):

        alpha = alpha.to(ureg.rad).magnitude #ensure its in radians

        A_fp = self.fins.planform_area
        A_fe = self.fins.exposed_area
        d_f = self.fins.tube_dia
        l_TS = d_f + 2*self.fins.span # total fin span
        n = self.fins.n

        Rs = l_TS/d_f

        k_fb = 0.8065*Rs**2 + 1.1553*Rs
        k_bf = 0.1935*Rs**2 + 0.8174*Rs + 1


        #TODO: check this equation (eqn. 50 in Box 2009) - the typography seems odd
        CD_falpha = alpha**2 * ((1.2 * n * A_fp) / (np.pi * d_f**2) + 3.12*(k_fb + k_bf - 1)*(n * A_fe)/(np.pi * d_f**2))

        return CD_falpha

    def CD0(self, Re=1e6):
        """Calcualte the zero angle-of-attack incompressible drag of the rocket.
        Generally uses DATCOM method (as specified by Box [1])
        Reynolds number refers to the reynolds number by the length of the rocket. """

        CD0_fb = self.CD0_fb(Re)
        CD0_b  = self.CD0_b(Re)

        CD0 = CD0_fb + CD0_b

        if self.fins is not None:
            CD0_f  = self.CD0_f(Re)
            CD0 = CD0 + CD0_f

        return CD0

    def CD0_fb(self, Re=1e6):
        """Calcualte the zero-angle of attack drag due to forebody of the rocket"""

        #total length of the rocket
        l_TR = self.length()

        if self.boat_tail is not None:
            # length of the boat tail
            l_c = self.boat_tail.length
            # diameter at boat tail
            d_d = self.boat_tail.aft_dia
        else:
            l_c = 0*ureg.m;
            d_d = 0*ureg.m;

        # maximum body diameter of rocket
        d_b = max([comp.diameter for comp in self.components if type(comp) is BodyTube])

        # length of body tube
        l_b = sum([comp.length for comp in self.components if type(comp) is BodyTube])

        # length of nose cone
        l_n = self.nose_cone.length

        # coefficient of friction of fore body
        Cf_fb = self.Cf(Re=Re)

        CD0_fb = (1+60/(l_TR/d_b)**3 + 0.0025*(l_b/d_b))*(2.7*(l_n/d_b) + 4*(l_b/d_b) + 2*(1 - d_d/d_b)*(l_c/d_b))*Cf_fb

        return CD0_fb

    def CD0_b(self, Re=1e6):
        """Calcualte the zero-angle of attack drag due to base drag"""
        # find max dia
        d_b = max([comp.diameter for comp in self.components if type(comp) is BodyTube])

        if self.boat_tail is None:
            d_d = self.body_tube.diameter
        else:
            # find boat dia aft dia
            d_d = self.boat_tail.aft_dia

        CD0_b = 0.029*(d_d/d_b)**3/(self.CD0_fb(Re))**0.5


        return CD0_b

    def CD0_f(self, Re=1e6):
        """Calcualte the zero-angle of attack drag due to the fins, including the effect of the interference"""

        if self.fins is None:
            raise RuntimeError("Please define the fins using rocket.set_fins(fins) first.")

        l_TR = self.length()

        l_m_fins = self.fins.mid_chord_span
        t_f = self.fins.thickness
        A_fp = self.fins.planform_area
        A_fe = self.fins.exposed_area
        d_f = self.fins.tube_dia
        n = self.fins.n

        Re_fins = Re*(l_m_fins/l_TR)
        Cf_f = self.Cf(Re_fins)


        CD0_f = 2 * Cf_f * (1 + 2*t_f/l_m_fins) * (4 * n * (2*A_fp - A_fe)) / (np.pi * d_f**2)

        return CD0_f


    def Cf(self, Re=1e6):
        """Return the viscous friction coefficient at a Reynolds number"""

        Re_c = 5e5; #critical reynolds number for transition

        if Re < Re_c:

            Cf = 1.328/np.sqrt(Re)

            return Cf

        else:
            B = Re_c*(0.074*Re**(-0.2) - 1.328*Re**(-0.5))

            Cf = 0.074*Re**(-0.2) - B/Re

            return Cf

    def describe(self, describe_components=False):


        print(f'Rocket: {self.name}')
        print()

        print('Rocket Details')

        x1 = PrettyTable()

        x1.field_names = ["Parameter", "Value", "Notes"]

        try:
            x1.add_row(["Total Mass", f'{self.mass():.4f~}', ""])
        except:
            x1.add_row(["Total Mass", 'ERROR', ""])

        try:
            x1.add_row(["Total Length", f'{self.length():.4f~}', ""])
        except:
            x1.add_row(["Total Length", 'ERROR', ""])

        try:
            x1.add_row(["X_CG", f'{self.xcg():4f~}', ""])
        except:
            x1.add_row(["X_CG", 'ERROR', ""])

        try:
            x1.add_row(["X_CP", f'{self.xcp():.4f~}', "At default values"])
        except:
            x1.add_row(["X_CP", 'ERROR', "At default values"])

        try:
            diameter = 2*(self.A_ref/np.pi)**0.5
            x1.add_row(["Static Margin (calibers)", f'{(self.xcp()-self.xcg())/self.diameter:.4f~}',"At default values"])
        except:
            x1.add_row(["Static Margin (calibers)", 'ERROR', "At default values"])


        try:
            x1.add_row(["CD", f'{self.CD():.4f~}',"At default values"])
        except:
            x1.add_row(["CD", 'ERROR',"At default values"])

        try:
            x1.add_row(["CNa", f'{self.CNa():.4f~}',"At default values"])
        except:
            x1.add_row(["CNa", 'ERROR',"At default values"])

        print(x1)

        print()

        print('Component Details')

        m = self.mass()

        x2 = PrettyTable()
        x2.field_names = ["Component", "Type", "Material", "Mass", "Mass Fraction %", "CNa"]
        for comp in self.components:
            x2.add_row([comp.name, comp.__class__.__name__, comp.material.name, f'{comp.mass:5.2f~}', f'{100*(comp.mass/m):.2f~}', f'{comp.CNa():.3f~}'])

        print(x2)

        print()
        print("Describing all components in full:")
        print()

        if describe_components:
            for comp in self.components:
                comp.describe()
