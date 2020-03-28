"""Utility functions for rocketPy"""

from . import Q_, ureg
import numpy as np
import numpy.linalg as la


def si(v):
    """Utility function to convert a Pint Quantity into a float in SI units.

    Args:
        v (Pint.Quantity or list/numpy.Array of Pint.Quantity): quantity or list of quantities to convert.

    Returns:
        float or list of floats: magnitudes in SI units.

    Examples
        Examples should be written in doctest format, and
        should illustrate how to use the function/class.
        >>> si(5.0*ureg.meter)
        5.0
        >> si(6.0*ureg.inch)
        0.1524
    """

    if isinstance(v, Quaternion):
        return v.q

    try:
        return v.to_base_units().magnitude
    except BaseException:
        try:
            return [si(vv) for vv in v]
        except BaseException:
            raise RuntimeError('Conversion to SI has failed')


def mach_correction(Ma=0.0, method='default'):
    """Performs the Prandtl-Glauert compressibility correction, extended for supersonic region.

    Args:
        Ma (dimensionless float): Mach Number. Defaults to 0.0.
        method (str): Choose the correction method ('default' or 'Cambridge'). Defaults to 'default'.

    Returns:
        dimensionless float: Mach correction multiplier (float)

    Examples
        Examples should be written in doctest format, and
        should illustrate how to use the function/class.
        >>> mach_correction(0.5)
        1.1547005383792517
        >> mach_correction(1.5)
        0.8944271909999159

    """

    if method == 'default':
        # my own correction:
        beta = max(np.sqrt(abs(1 - Ma**2)), np.sqrt(1 - 0.8**2))

        return 1 / beta

    if method == 'Cambridge':
        # Follows the formulation of eqn. 55-57 of Box:
        if Ma < 0.8:
            return 1 / np.sqrt(1 - Ma**2)
        elif Ma > 1.1:
            return 1 / np.sqrt(Ma**2 - 1)
        else:
            return 1 / np.sqrt(1 - 0.8**2)


def unit_vector(v):
    """Return a unit vector in the direction of v."""
    return v / la.norm(v)


def angle_between(va, vb):
    """ Return the angle between two column vectors using the dot product"""

    th = np.arccos(unit_vector(va) @ unit_vector(vb))

    return float(th)


class Quaternion():

    """Note, does not support Pint units """

    def __init__(self, s, vx, vy, vz):
        """ Creates a quaternion. Expects floats"""
        self.q = np.array([s, vx, vy, vz])

    @classmethod
    def from_angle(cls, theta, axis):
        """A rotation of angle `theta` about the axis `ax, ay, az`. Allows the axis to be not normalized"""

        axis = unit_vector(axis)

        s = np.cos(theta / 2)
        v = np.sin(theta / 2) * axis

        return cls(s, *v)

    def rot_matrix(self, qq=None):

        if qq is None:
            qq = self.q

        R = np.array([[1 -
                       2 *
                       qq[2]**2 -
                       2 *
                       qq[3]**2, 2 *
                       qq[1] *
                       qq[2] -
                       2 *
                       qq[0] *
                       qq[3], 2 *
                       qq[1] *
                       qq[3] +
                       2 *
                       qq[0] *
                       qq[2]], [2 *
                                qq[1] *
                                qq[2] +
                                2 *
                                qq[0] *
                                qq[3], 1 -
                                2 *
                                qq[1]**2 -
                                2 *
                                qq[3]**2, 2 *
                                qq[2] *
                                qq[3] -
                                2 *
                                qq[0] *
                                qq[1]], [2 *
                                         qq[1] *
                                         qq[3] -
                                         2 *
                                         qq[0] *
                                         qq[2], 2 *
                                         qq[2] *
                                         qq[3] +
                                         2 *
                                         qq[0] *
                                         qq[1], 1 -
                                         2 *
                                         qq[1]**2 -
                                         2 *
                                         qq[2]**2]])

        return R

    def rate_of_change(self, omega):
        """Return the rate of change of the quaternion based on an angular velocity"""
        # follows the formulation in Box
        s = self.q[0]
        v = self.q[1:4]

        sdot = 0.5 * omega @ v
        vdot = 0.5 * (s * omega + np.cross(omega, v))

        return np.hstack([[sdot], vdot])

    def __repr__(self):
        return str(self.q)
