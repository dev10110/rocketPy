#!/usr/bin/env python

"""Tests for `rocketPy.Quaternion`"""

import pytest

import rocketPy as rp
from rocketPy.quaternion import Quaternion

import numpy as np
import numpy.testing as nptest

# test if the creation is good:
##


def test_generate_array():

    nptest.assert_allclose(Quaternion(), Quaternion([0, 1, 0, 0]))

    # check the normalization
    nptest.assert_allclose(Quaternion(
        [3, -4, 0, 0]), Quaternion([0.6, -0.8, 0, 0]))
    nptest.assert_allclose(Quaternion(
        [3, 4, 0., 0]), Quaternion([0.6, 0.8, 0, 0]))
    nptest.assert_allclose(Quaternion(
        [3, -4, 0, 0]), np.array([0.6, -0.8, 0, 0]))

    with pytest.raises(ValueError):
        q1 = Quaternion([0, 0, 0, 0])


def test_from_angle():

    q1 = Quaternion.from_angle(np.pi / 2, [1, 0, 0])
    r2o2 = (1 / 2)**0.5  # root 2/2
    nptest.assert_allclose(q1, np.array([r2o2, r2o2, 0, 0]))

    q2 = Quaternion.from_angle(np.pi, [1, 1, 1])

    n = np.array([1, 1, 1])
    r1o3 = (1 / 3)**0.5

    nptest.assert_allclose(q2, np.array([0, r1o3, r1o3, r1o3]), 1e-7, 1e-7)


def test_rotation_matrix():

    q1 = Quaternion.from_angle(np.pi / 4, [1, 0, 0])
    q2 = Quaternion.from_angle(np.pi, [0, 1, 0])
    q3 = Quaternion.from_angle(30 * np.pi / 180, [0, 1, 1])
    q4 = Quaternion.from_angle(30 * np.pi / 180, [0, 2, 1])

    r1 = q1.rot_matrix()
    r2 = q2.rot_matrix()
    r3 = q3.rot_matrix()
    r4 = q4.rot_matrix()

    R1 = np.array([[1.0000000,  0.0000000,  0.0000000],
                   [0.0000000,  0.7071068, -0.7071068],
                   [0.0000000,  0.7071068,  0.7071068]])

    R2 = np.array([[-1.0000000,  0.0000000,  0.0000000],
                   [0.0000000,  1.0000000,  0.0000000],
                   [-0.0000000,  0.0000000, -1.0000000]])

    R3 = np.array([[0.8660254, -0.3535534,  0.3535534],
                   [0.3535534,  0.9330127,  0.0669873],
                   [-0.3535534,  0.0669873,  0.9330127]])

    R4 = np.array([[0.8660254, -0.2236068,  0.4472136],
                   [0.2236068,  0.9732051,  0.0535898],
                   [-0.4472136,  0.0535898,  0.8928203]])

    nptest.assert_allclose(r1, R1, 1e-7, 1e-7)
    nptest.assert_allclose(r2, R2, 1e-7, 1e-7)
    nptest.assert_allclose(r3, R3, 1e-7, 1e-7)
    nptest.assert_allclose(r4, R4, 1e-7, 1e-7)


def test_axis_angle():

    q1 = Quaternion.from_angle(np.pi / 4, [1, 0, 0])
    q2 = Quaternion.from_angle(np.pi, [0, 1, 0])
    q3 = Quaternion.from_angle(30 * np.pi / 180, [0, 1, 1])
    q4 = Quaternion.from_angle(30 * np.pi / 180, [0, 2, 1])

    nptest.assert_allclose(q1.axis_angle()[0], [1, 0, 0], 1e-7, 1e-7)
    nptest.assert_allclose(q1.axis_angle()[1], np.pi / 4, 1e-7, 1e-7)

    nptest.assert_allclose(q2.axis_angle()[0], [0, 1, 0], 1e-7, 1e-7)
    nptest.assert_allclose(q2.axis_angle()[1], np.pi, 1e-7, 1e-7)

    nptest.assert_allclose(q3.axis_angle()[0], [
                           0, 1 / 2**0.5, 1 / 2**0.5], 1e-7, 1e-7)
    nptest.assert_allclose(q3.axis_angle()[1], np.pi / 6, 1e-7, 1e-7)

    nptest.assert_allclose(q4.axis_angle()[0], [
                           0, 2 / 5**0.5, 1 / 5**0.5], 1e-7, 1e-7)
    nptest.assert_allclose(q4.axis_angle()[1], np.pi / 6, 1e-7, 1e-7)

    #nptest.assert_allclose(r2,R2, 1e-7,1e-7)
    #nptest.assert_allclose(r3,R3, 1e-7,1e-7)
    #nptest.assert_allclose(r4,R4, 1e-7,1e-7)


def test_derivative():
    pass


##
