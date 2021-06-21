""" This script get input from user to calculate the mean power production of
a satellite in orbit.
Author: Yannis Koveos (https://gitlab.com/ykoveos)
License: GNU GENERAL PUBLIC LICENSE Version 3
Requirements:
* SGP4: https://pypi.python.org/pypi/sgp4/
* PyKDL: https://pypi.org/project/PyKDL/
* numpy: version 1.16.2
* matplotlib : version 2.2.2
"""

from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from sgp4.ext import invjday
import PyKDL
import numpy as np
from math import ceil, pi
from tqdm import tqdm

axis_labels = ["X+", "X-", "Y+", "Y-", "Z+", "Z-"]


def sunpos_ECI(JD):
    """Calculate sun position vector.

    @package This algorithm is from "Fundamentals of Astrodynamics and
    Applications", David A. Vallado, 1997, McGraw-Hill, Alg. 18, pg.183.
    For quick debug use Ex. 3-8, or use an online sun position calculator
    @param Julian Date
    @return d,r Earth Sun distance (AU), Sun position vector in ECI (AU)
    """
    T_UT1 = (JD - 2451545.0) / 36525
    l_m = 280.4606184 + 36000.77005361 * T_UT1  # Mean longitude of the sun
    l_m = np.fmod(l_m, 360.0)
    if l_m < 0:
        l_m = l_m + 360.0
    T_TDB = T_UT1

    man_sun = 357.5277233 + 35999.05034 * T_TDB  # Mean anomaly of Sun
    man_sun = np.fmod(man_sun, 360.0)
    if man_sun < 0:
        man_sun = man_sun + 360.0
    man_sun = man_sun * np.pi / 180.0  # In RADIANS

    # Ecliptic longitude of Sun, in RADIANS
    l_ecl = (l_m + 1.914666471 * np.sin(man_sun) +
             0.019994643 * np.sin(2 * man_sun)) * np.pi / 180.0
    # Distance to Sun in AUs
    d = 1.000140612 - 0.016708617 * np.cos(man_sun)
    -0.000139589 * np.cos(2 * man_sun)
    eps = (23.439291 - 0.0130042 * T_TDB) * np.pi / 180.0  # In RADIANS
    r = (d * np.cos(l_ecl), d * np.cos(eps) * np.sin(l_ecl),
         d * np.sin(eps) * np.sin(l_ecl))
    return d, r


def get_orbit_period(TLE):
    """Determine orbit period in munutes from TLE."""
    satellite = twoline2rv(TLE[1], TLE[2], wgs72)
    return (2 * pi) / satellite.no_kozai


def power_calculation(TLE, JD_ini, orbits, ang_v, pose, pv_eff, pv_count,
                      pv_area, timestep=1, show_progress=False):
    """
    Calculate accumulated power for each side of the satellite.

            Parameters:
                    TLE (string list): 2 line orbit elements
                    JD_ini: Julian day for simulation start
                    orbits (int): Number of orbits to simulate
                    ang_v (list): Angular velocity for XYZ
                    pose (list): Initial pose for XYZ
                    pv_eff (list): PV Efficiency
                    pv_count (list): Number of PVs per side
                    pv_area (list): Area of each PV element
                    timestep (int): Simulation time step in seconds
                    show_progress (bool): Display progress bar

            Returns:
                    mean_power_coeff (list): Mean power coefficient for each axis
                    mean_power (list): Mean power for each axis
                    side_power_values (list): Power values for each axis
                    orbit_side_mean_power (list): Mean power per orbit for each axis
    """
    satellite = twoline2rv(TLE[1], TLE[2], wgs72)
    # Solar irradiance (kW/m^2) in a specific orbit
    si = 1.4
    orbit_period = get_orbit_period(TLE)*(60/timestep)
    total_time = ceil(orbits * orbit_period)  # in min

    # Mean power of each side in mW
    mean_power = np.array([0, 0, 0, 0, 0, 0])
    # Earth radius in km
    r_earth = 6371.0
    # Power coefficient of each side
    fxp = np.zeros(total_time)
    fxm = np.zeros(total_time)
    fyp = np.zeros(total_time)
    fym = np.zeros(total_time)
    fzp = np.zeros(total_time)
    fzm = np.zeros(total_time)
    # Orbit and daily means
    orbit_side_mean_coeff = []
    # Emulate the orbit
    if show_progress:
        pbar = tqdm(total=total_time,
                    desc=TLE[0],
                    bar_format="{l_bar}{bar} [ time left: {remaining} ]")
    eclipse = []
    eclipse_pc = []     # eclipse percentage per orbit
    for curr_time in range(0, total_time):
        # Time in Julian days, curr_time in minutes
        JD = JD_ini + curr_time / (1440.0*60/timestep)
        y, mon, d, h, mn, sec = invjday(JD)
        # d_sun Earth Sun distance (AU) , Sun position vector r_sun in ECI (AU)
        d_sun, r_sun = sunpos_ECI(JD)
        r_sun = PyKDL.Vector(r_sun[0], r_sun[1], r_sun[2])
        r_sun.Normalize()
        # Positon in km and vel in km/s from the center of the earth in ECI
        r, v = satellite.propagate(y, mon, d, h, mn, sec)
        # print(r,v, np.linalg.norm(r))
        # Orbital frame: Z is down, X is in velocity, projected to be perpendigular
        # to Z (since non circular orbits), in the same velocity - down plane
        # (which should be the orbital plane)
        Zof_eci = PyKDL.Vector(-r[0], -r[1], -r[2])
        Zof_eci.Normalize()
        Xof_eci = PyKDL.Vector(v[0], v[1], v[2])
        Xof_eci.Normalize()
        # Fix for non circular orbit, when Zof_eci and Xof_eci aren't perpedicular
        Xof_eci = Xof_eci - PyKDL.dot(Zof_eci, Xof_eci) * Zof_eci
        Xof_eci.Normalize()
        # Yof_eci is always perpedicular to the plane of Zof_eci and Xof_eci
        Yof_eci = Zof_eci * Xof_eci
        # Rotation matrix that transforms a vector from the orbit frame to ECI
        Rof = PyKDL.Rotation(Xof_eci, Yof_eci, Zof_eci)
        # Rotation matrix that transforms a vector from Body Frame to Orbital Frame
        Rbf = PyKDL.Rotation(PyKDL.Vector(1.0, 0.0, 0.0),
                             PyKDL.Vector(0.0, 1.0, 0.0),
                             PyKDL.Vector(0.0, 0.0, 1.0))
        # Add angular velocities in body frame for each axis, integration time 1min
        pose = (pose + ang_v * 1.0)
        Rbf.DoRotX(pose[0])
        Rbf.DoRotY(pose[1])
        Rbf.DoRotZ(pose[2])
        # Body vector in body frame
        Xbf = PyKDL.Vector(1.0, 0.0, 0.0)
        Ybf = PyKDL.Vector(0.0, 1.0, 0.0)
        Zbf = PyKDL.Vector(0.0, 0.0, 1.0)
        # Body vector in ECI
        Xbf_eci = Rof * (Rbf * Xbf)
        Ybf_eci = Rof * (Rbf * Ybf)
        Zbf_eci = Rof * (Rbf * Zbf)
        # Satellite vector in ECI to check if the sun is behind the earth
        r_sat = PyKDL.Vector(-r[0], -r[1], -r[2])
        # Check if the sun is behind the earth (Eclipse)
        tmp = PyKDL.dot(r_sat, r_sun)
        if ((r_earth > (r_sat - tmp * r_sun).Norm()) and (tmp > 0)):
            r_sun = 0.0 * r_sun
            eclipse.append(True)
        else:
            eclipse.append(False)
        # Calculate the projection of sun vector in each side of satellite,
        # power efficiency
        # X axis
        tmp = PyKDL.dot(Xbf_eci, r_sun)
        if tmp >= 0.0:
            fxp[curr_time] = tmp
        else:
            fxm[curr_time] = -tmp
        # Y axis
        tmp = PyKDL.dot(Ybf_eci, r_sun)
        if tmp >= 0.0:
            fyp[curr_time] = tmp
        else:
            fym[curr_time] = -tmp
        # Z axis
        tmp = PyKDL.dot(Zbf_eci, r_sun)
        if tmp >= 0.0:
            fzp[curr_time] = tmp
        else:
            fzm[curr_time] = -tmp

        # Per orbit calculations
        if (curr_time > 0) and (((curr_time + 1) % ceil(orbit_period)) == 0):
            # Store orbit mean pwr
            orbit_side_mean_coeff.append(
                calc_avg_power_coeff(
                    (fxp[(curr_time - ceil(orbit_period) + 1):curr_time],
                     fxm[(curr_time - ceil(orbit_period) + 1):curr_time],
                     fyp[(curr_time - ceil(orbit_period) + 1):curr_time],
                     fym[(curr_time - ceil(orbit_period) + 1):curr_time],
                     fzp[(curr_time - ceil(orbit_period) + 1):curr_time],
                     fzm[(curr_time - ceil(orbit_period) + 1):curr_time]),
                    ceil(orbit_period)))
            # Get eclipse time
            eclipse_pc.append(sum(eclipse)/orbit_period)
            eclipse = []

        if show_progress:
            pbar.update(1)
    p_eff_values = (fxp, fxm, fyp, fym, fzp, fzm)
    mean_power_coeff = calc_avg_power_coeff(p_eff_values, total_time)
    # PV mean power production
    mean_power = pv_area * pv_eff * pv_count * mean_power_coeff * si
    # Plot the power coefficient
    side_power_values = []
    orbit_side_mean_power = []
    for axis, side in enumerate(p_eff_values):
        side_power_values.append(side * si * pv_area[axis] * pv_eff[axis] * pv_count[axis])
    for axis, side in enumerate(np.transpose(orbit_side_mean_coeff)):
        orbit_side_mean_power.append(side * si * pv_area[axis] * pv_eff[axis] * pv_count[axis])
    return (mean_power_coeff, mean_power, side_power_values,
            orbit_side_mean_power, eclipse_pc)


def calc_avg_power_coeff(pwr_eff, total_time):
    """Calculate average power coefficient.

    Calculate the integral of power coefficient to get the energy coefficient
    and then the division with total time to get the average power coefficient
    """
    mean_power_coeff = []
    for axis_pwr in pwr_eff:
        mean_power_coeff.append(np.trapz(
            axis_pwr) / total_time)
    return mean_power_coeff
