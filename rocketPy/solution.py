import numpy as np
#from collections.abc import Iterable


class Solution():

    def __init__(self, sol, stage):
        """Wrapper for np.ODEresult to store extra information.
        Most importantly, allows multiple ODEresults to be chained together,
        with potential for storing extra information in them, and making it easy
        to access properties from each.

        Args:
            sol (ODEresult or List of ODEresult): The ODEresult objects to be stored
            stage (int or list of int): Corresponding stage counters

        Returns:
            Solution: Solution object
        """




        if type(sol) is list:
            # # if a list is passed, save it
            self.sols = sol
            self.stages = stage
        else:
            #if a single OdeResult, save just the one object as a list
            self.sols = [sol]
            self.stages = [stage]


    def DOF(self):
        """Returns the number of degrees of freedom in the state vector

        Returns:
            int: length of state vector
        """

        # access the first solution and get the shape
        return self.sols[0].y.shape[0]

    def sol(self, time, error='raise'):
        """Returns the state at requested times.
        Assumes that scipy.integrate was called with dense_output=True and thus uses scipy's interpolation method

        Args:
            time (float or np.array): Float time to request the solution for, or list/array of times that the state is requested for.
            error ('raise' or numeric): if the requested time is outside the interpolations capabilities, if error='raise', will raise a warning, else will return error in the state. Maintains the shape of the state vector.
        Returns:
            np.array: state vector. If time is float, returns simple array of the right length. if time is a list of length n and state has length m, returns a shape of size (n x m).

        """
        # convert to numpy object
        time  = np.asarray(time)

        # check if the request was for a single point of time
        if time.ndim == 0:
            return self._sol(time, error)

        # if not, iterate over the entire list of times
        # and collect the solutions into a list
        ys = [self._sol(t, error) for t in time]

        # transpose and return the state vectors
        return np.transpose(np.vstack(ys))

    def _sol(self, time, error='raise'):
        """Equivalent to Simulation.sol, except only works for single time point.

        Args:
            time (float): time
            error ('raise' or numeric): if the requested time is outside the interpolations capabilities, error='raise' will raise a warning, else will return [error] in the state.

        Returns:
            np.array: state vector

        """

        # check if the requested time is in the bounds:
        if time < self.t_min() or time > self.t_max():
            if error is 'raise':
                raise ValueError("Request time is outside computed times")
            else:
                # return a numpy array the size of the state
                return np.array([error,]*self.DOF())

        # else, it is in the computed times
        # so compute it
        # work backwards. check each stage, if the stage has the solution return it
        for sol in reversed(self.sols):
            if time >= sol.sol.t_min and time <= sol.sol.t_max:
                return sol.sol(time)

    def t_min(self):
        """Minimum time of all the solutions in this object.
        Computes a simple minimum, gaps are not accounted for.

        Returns:
            float: minimum time
        """
        return min(sol.sol.t_min for sol in self.sols)

    def t_max(self):
        """Maximim time of all the solutions in this object.
        Computes a simple maximum, gaps are not accounted for.

        Returns:
            float: maximum time
        """
        return max(sol.sol.t_max for sol in self.sols)

    def __add__(self, other):
        """Overloading the + operator to allow solutions to be 'added' together, essentially chaining them.

        Args:
            other (Solution)

        Returns:
            Solution: Solution object with both objects inside
        """


        if self.DOF() != other.DOF():
            raise RuntimeError("Adding Solution objects requires them to have the same state vector size")

        # add the two lists together
        stages = self.stages + other.stages
        sols = self.sols + other.sols

        # create a new solution
        newBranch = Solution(sols, stages)

        return newBranch

    def __repr__(self):
        """Display string"""
        return f'Solution:\n[\n Stages: {self.stages}\n ODEresults: {self.sols} \n ]\n'
