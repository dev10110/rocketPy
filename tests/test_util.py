#!/usr/bin/env python

"""Tests for `rocketPy` package."""

import pytest

import rocketPy as rp
from rocketPy import ureg


def test_si():

    assert rp.util.si(5.0*ureg.m) == 5.0


def test_materials():


    assert True
