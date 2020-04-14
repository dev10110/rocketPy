# class to define simulations for dynamic objects
#

import numpy as np
import scipy.integrate as spint

from .solution import Solution

import logging



class Simulation():
    def __init__(self, object):

        # store the object that we will do dynamics on
        self.object = object

    def solve(self, t_span, y0, stage, user_events=[], **kwargs):
        """Thin wrapper for scipy.solve_ivp to handle stages and events more intuitively.
        Arguments generally follow scipy.solve_ivp.
        The dynamics are picked up from object.dynamics, and dense_output=True to allow for solutions to be queried later
        Solve method uses scipy default, but can be specified using kwargs

        Args:
            t_span ([float, float]): Simulation start and stop time.
            y0 (np.array): Initial state vector
            stage (float): stage number to use for the solve
            user_events (list of functions): each function is an event function. Defaults to [].

        Returns:
            rocketPy.Solution: solution object with the stage stored

        """

        # solve using solve_ivp
        sol = spint.solve_ivp(self.object.dynamics, t_span, y0, args=(stage,), events=user_events, dense_output=True, **kwargs)

        # return solution
        return Solution(sol, stage)


    def nominal_solve(self, t_span, y0, starting_stage, **kwargs):
        """Solves for a nominal flight

        Args:
            t_span ([float, float]): Timespan to simulate over
            y0 (np.array): Initial state
            starting_stage (int): which stage to start in


        Returns:
            rocketPy.Solution: Solution object

        """

        # simulate this stage

        logging.debug(f'\n Nominal Simulating: T: {t_span}, Y0: {y0}, Stage: {starting_stage}')

        sol = self.solve(t_span, y0, starting_stage, self.object.staging_functions, **kwargs)

        logging.debug(f'{sol.sols[0].t_events}')

        # for each of the staging events add the next stage
        # determine the earliest event, and use that to start the next stage

        min_i = np.inf
        min_j = np.inf
        min_t = np.inf

        # todo (low): refactor to simplify
        #
        # for each of the staging functions
        for i in range(len(self.object.staging_functions)):
            # if the event is relevant to the current stage
            if starting_stage in self.object.staging_functions[i].trigger_if_stage_in:
                # over all the times that this event was detected
                for j in range(len(sol.sols[0].t_events[i])):
                    # figure out if we need to update which is the minimum
                    if sol.sols[0].t_events[i][j] < min_t:
                        min_i = i
                        min_j = j
                        min_t = sol.sols[0].t_events[i][j]

        logging.debug(f'Determined min i={min_i},j={min_j}, t={min_t}')

        # if there were no t events, just return sol
        if min_t == np.inf:
            return sol

        # now figure out which is the next phase
        next_stage = self.object.staging_functions[min_i].nominal_next_stage

        # if there is no next stage despite the event, return sol
        if next_stage is None:
            return sol

        # figure out the new start time
        new_t_start = sol.sols[0].t_events[min_i][min_j] + self.object.staging_functions[min_i].t_offset
        new_t_span = [new_t_start, t_span[-1]];

        # determine if the state needs to be changed (allows for step changes in the state!)
        if self.object.staging_functions[min_i].modify_state is None:
            new_y0 = sol.sol(new_t_start) # unmodified for now
        else:
            # apply the change
            new_y0 = self.object.staging_functions[min_i].modify_state(self.object, sol.sol(new_t_start))

        # recursively call this function, and chain the stages
        return sol + self.nominal_solve(new_t_span, new_y0, next_stage)

    def full_solve(self, t_span, y0, starting_stage, **kwargs):
        """Solves every single possible outcome.

        #todo (low): add probability weighting

        Args:
            t_span ([float, float]): Timespan for the solve
            y0 (np.array): Initial state vector
            starting_stage (int): Initial stage to start from

        Returns:
            list of rocketPy.Solution: list of possible outcomes
        """

        # simulate this stage
        logging.debug(f'\n Full Simulating: T: {t_span}, Y0: {y0}, Stage: {starting_stage}')

        # solve the current stage
        sol = self.solve(t_span, y0, starting_stage, self.object.staging_functions, **kwargs)
        logging.debug('Got: ', sol.sols[0].t_events)

        # create a list with the solution
        # this allows for the thing to ignore the event detection
        new_sols = [sol]

        # ooph combinatorics might get bad on this problem
        # todo (low): maybe add a probability threshold

        # for each type of event function,
        for i in range(len(self.object.staging_functions)):
            # if this stage allows this trigger
            if starting_stage in self.object.staging_functions[i].trigger_if_stage_in:
                # determine the next stages
                next_stages = self.object.staging_functions[i].possible_next_stages
                # if there is a next stage
                if len(next_stages)>0:
                    # for each next stage
                    for next_stage in next_stages:
                        # for each time that the trigger was detected
                        for j in range(len(sol.sols[0].t_events[i])):
                            # as long as the event doesnt happen at the start of this simulation (as this is a false event trigger)
                            if sol.sols[0].t_events[i][j] > t_span[0]:

                                # create a branch
                                logging.debug(f'For tspan={t_span} and stage={starting_stage}, About to use i,j = i={i},j={j}, next_stage={next_stage} of {next_stages}')

                                # determine new t
                                new_t_start = sol.sols[0].t_events[i][j] + self.object.staging_functions[i].t_offset
                                new_t_span = [new_t_start, t_span[-1]];

                                # as long as the time offsets dont go outside the simulated range
                                if new_t_start > sol.t_min() and new_t_start < sol.t_max():

                                    # determine new y0
                                    if self.object.staging_functions[i].modify_state is None:
                                        new_y0 = sol.sol(new_t_start)
                                    else:
                                        # apply the modification
                                        new_y0 = self.object.staging_functions[i].modify_state(self.object, sol.sol(new_t_start))
                                    # append a new solution
                                    new_sols.extend([sol + s for s in self.full_solve(new_t_span, new_y0, next_stage, **kwargs)])
        # return
        return new_sols
