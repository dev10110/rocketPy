#!/usr/bin/env python

"""Tests for `rocketPy` package."""

import pytest

import rocketPy as rp
from rocketPy import ureg

import numpy.testing as npt


def test_materials():

    #create al:

    al = rp.materials.Aluminium()
    assert al.max_temp == 420*ureg.degK

    #create polycarb:
    pc = rp.materials.Polycarbonate()
    npt.assert_almost_equal(pc.density.to_base_units().m, (1200*ureg.kg/(ureg.m**3)).to_base_units().m)
