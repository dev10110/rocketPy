
"""This example file demonstrates the construction of a simple rocket.
The corresponding .ork file is the OpenRocket implementation of the same rocket for comparison.
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle

import rocketPy as rp
from rocketPy import ureg


## First, create a rocket.
r = rp.Rocket(name='Sporadic Impulse')

## create a nose cone
nc = rp.NoseCone(name='Nose Cone', diameter=6*ureg.inch, length=48*ureg.cm, material=rp.materials.PLA())
# assign to rocket
r.set_nose_cone(nc)

## create a BodyTube
bt = rp.BodyTube(name = 'Body Tube', diameter=6*ureg.inch, length=(354-48)*ureg.cm, wall_thickness=2*ureg.mm, material=rp.materials.Phenolic())
# define its location
bt.set_position(after=nc)
# assign to rocket
r.set_body_tube(bt)

# create a boat tail
boat_tail = rp.Transition(name='Boat Tail', fore_dia=6*ureg.inch, aft_dia=6*ureg.inch, length=0.01*ureg.inch, material=rp.materials.Phenolic())
# define its location
boat_tail.set_position(after=bt)
# assign to rocket
r.set_boat_tail(boat_tail)

## create the fins
fins = rp.FinSet(name='Fins', n=4, span=6*ureg.inch, root_chord=12*ureg.inch, tip_chord=6*ureg.inch, mid_sweep=10*ureg.degree, tube_dia=6*ureg.inch,  thickness=2*ureg.mm, material=rp.materials.Aluminium())
# define its location
fins.set_position(end_of=bt, offset=-fins.root_chord)
#  assign to rocket
r.set_fins(fins)

# plot the entire rocket
fig = plt.figure()
ax = plt.gca()
r.plot(ax, unit=ureg.inch)
plt.draw()

# describe the rocket
r.describe(describe_components=True)


pickle.dump(r,open('SporadicImpulse.p','wb'))
pickle.dump(ureg,open('SporadicImpulse_ureg.p','wb'))

plt.show()
