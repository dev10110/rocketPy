import numpy as np

class Quaternion(np.ndarray):
    """Only works with unit quaternions (as needed to describe rotations)"""

    def __new__(cls, input_array=[0.,1.,0.,0.]):
        """By default works with 0 rotation about x axis"""

        n = np.linalg.norm(input_array)
        if n == 0.:
            raise ValueError("Quaternion axis must not be 0 (the norm of the [1,2,3] components of the quaternion cannot be 0.)")
        input_array = np.array(input_array)/np.linalg.norm(input_array)

        obj = np.asarray(input_array).view(cls)

        return obj

    def __array__finalize(self, obj):
        pass


    @classmethod
    def from_angle(cls, theta, axis):
        """A rotation of angle `theta` about the axis `ax, ay, az`. Allows the axis to be not normalized"""

        axis = np.array(axis)/np.linalg.norm(axis)

        s = np.cos(theta / 2)
        v = np.sin(theta / 2) * axis

        return cls([s, *v])

    def normalize(self):
        """Convert quaternion to unit vector"""

        self/np.linalg.norm(self)

    def rot_matrix(self):
        """Generate a rotation matrix"""

        R = np.array([[1 - 2 * self[2]**2 - 2 * self[3]**2,
                       2 * self[1] * self[2] - 2 * self[0] * self[3],
                       2 * self[1] * self[3] + 2 * self[0] * self[2]],
                [2 * self[1] * self[2] + 2 * self[0] * self[3],
                 1 - 2 * self[1]**2 - 2 * self[3]**2,
                 2 * self[2] * self[3] - 2 * self[0] * self[1]],
                [2 * self[1] * self[3] - 2 * self[0] * self[2],
                 2 * self[2] * self[3] + 2 * self[0] * self[1],
                 1 - 2 * self[1]**2 - 2 * self[2]**2]])

        return R

    def axis_angle(self):
        """Returns the axis of rotation and the angle of rotation (radians) as a tuple"""

        ax = np.array(self[1:])
        ax = ax/np.linalg.norm(ax)

        angle = 2*np.arccos(self[0])

        return ax, angle


    def rate_of_change(self, omega):
        """Return the rate of change of the quaternion based on an angular velocity"""

        # follows the formulation in Box
        # note that there are other methods, which may be more accurate

        omega = np.array(omega)

        s = self[0]
        v = np.array(self[1:])

        sdot = - 0.5 * (omega @ v) # mistake in Box eqn 7 (based on https://www.sciencedirect.com/topics/computer-science/quaternion-multiplication and http://web.cs.iastate.edu/~cs577/handouts/quaternion.pdf)
        vdot = 0.5 * (s * omega + np.cross(omega, v))

        # returns a numpy object, not a quaternion
        return np.array([sdot, *vdot])

    def __mul__(self, other):
        raise NotImplementedError
