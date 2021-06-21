# import solar_power_budget as spb
import satellite_solar_power_budget.solar_power_budget as spb
from sgp4.ext import jday
import numpy as np
import matplotlib.pyplot as plt
import yaml
import os.path

skip_plot_axis = []
axis_colors = [
    'red', 'darkred', 'limegreen', 'darkgreen', 'cornflowerblue', 'mediumblue'
]


def plot_stacked(plot,
                 title,
                 TLE,
                 values,
                 skip_plot_axis,
                 xlabel="Time (min)"):
    """Create stacked axis plot."""
    plot.figure(title, figsize=(8, 5))
    plot.xlabel(xlabel)
    plot.ylabel("Power (mW)")
    plot.title('\n'.join(TLE), loc='left')
    for axis, side in enumerate(values):
        if axis not in skip_plot_axis:
            plot.plot(np.arange(0, len(side)),
                      side,
                      axis_colors[axis],
                      label=spb.axis_labels[axis] + ': ' + f'{min(side):.1f}' +
                      f'/{max(side):.1f}' + f'/{mean_power[axis]:.1f} mW')
    plot.legend(title='Average power (min/max/mean)', loc='upper right')


def plot_split(plot, title, TLE, values, skip_plot_axis):
    """Create split axis plot."""
    fig, axs = plot.subplots(6 - len(skip_plot_axis),
                             sharex=True,
                             sharey=True,
                             num=title,
                             figsize=(8, 9))
    fig.suptitle('\n'.join(TLE))
    fig.text(0.5, 0.04, 'Time (min)', ha='center')
    fig.text(0.04, 0.5, 'Power (mW)', va='center', rotation='vertical')
    idx = 0
    for axis, side in enumerate(side_power_values):
        if axis not in skip_plot_axis:
            axs[idx].plot(np.arange(0, len(side)),
                          side,
                          axis_colors[axis],
                          label=spb.axis_labels[axis] + ': ' +
                          f'{min(side):.1f}' + f'/{max(side):.1f}' +
                          f'/{mean_power[axis]:.1f} mW')
            axs[idx].legend(title='Average power (mW)', loc='upper right')
            idx += 1



with open('../../../satellite-parameters.yaml', mode='r') as file:
    try:
        satellite_parameters = yaml.load(file, Loader=yaml.FullLoader)

        # PV efficient in each side [X+, X-, Y+, Y-, Z+, Z-]
        pv_eff = np.array([
            satellite_parameters['SolarPanels']['efficiency']['X+'],
            satellite_parameters['SolarPanels']['efficiency']['X-'],
            satellite_parameters['SolarPanels']['efficiency']['Y+'],
            satellite_parameters['SolarPanels']['efficiency']['Y-'],
            satellite_parameters['SolarPanels']['efficiency']['Z+'],
            satellite_parameters['SolarPanels']['efficiency']['Z-']
        ])
        # Number of PV in each side [X+, X-, Y+, Y-, Z+, Z-]
        pv_count = np.array([
            satellite_parameters['SolarPanels']['count']['X+'],
            satellite_parameters['SolarPanels']['count']['X-'],
            satellite_parameters['SolarPanels']['count']['Y+'],
            satellite_parameters['SolarPanels']['count']['Y-'],
            satellite_parameters['SolarPanels']['count']['Z+'],
            satellite_parameters['SolarPanels']['count']['Z-']
        ])
        # Active area of each PV in mm^2, [X+, X-, Y+, Y-, Z+, Z-]
        pv_area = np.array([
            satellite_parameters['SolarPanels']['area']['X+'],
            satellite_parameters['SolarPanels']['area']['X-'],
            satellite_parameters['SolarPanels']['area']['Y+'],
            satellite_parameters['SolarPanels']['area']['Y-'],
            satellite_parameters['SolarPanels']['area']['Z+'],
            satellite_parameters['SolarPanels']['area']['Z-']
        ])
    except Exception as e:
        print('Error reading satellite parameters', e)
        exit(1)

with open('simulation-parameters.yaml', mode='r') as file:
    try:
        simulation_parameters = yaml.load(file, Loader=yaml.FullLoader)
        simulation_parameters['TLE'].insert(0, 'TLE')
        TLE = tuple(simulation_parameters['TLE'])

        # Initialize angular velocities in rad/min [X, Y, Z]
        ang_v = np.array([
            simulation_parameters['Spacecraft']['ang_v']['X'],
            simulation_parameters['Spacecraft']['ang_v']['Y'],
            simulation_parameters['Spacecraft']['ang_v']['Z']
        ])
        # Initial angle conditions in deg [X, Y, Z]
        pose0 = np.deg2rad([
            simulation_parameters['Spacecraft']['pose0']['X'],
            simulation_parameters['Spacecraft']['pose0']['Y'],
            simulation_parameters['Spacecraft']['pose0']['Z']
        ])
        # Set the start day for simulation, Julian day
        JD_ini = jday(simulation_parameters['Sim']['JD_start']['year'],
                      simulation_parameters['Sim']['JD_start']['month'],
                      simulation_parameters['Sim']['JD_start']['day'],
                      simulation_parameters['Sim']['JD_start']['hour'],
                      simulation_parameters['Sim']['JD_start']['min'],
                      simulation_parameters['Sim']['JD_start']['second'])
        orbit_count = simulation_parameters['Sim']['orbits']
        for axis, skip in enumerate(simulation_parameters['Sim']['skip_plot_axis']):
            if (skip):
                skip_plot_axis.append(axis)
        try:
            output_path = os.path.expanduser(simulation_parameters['Sim']['output_path'])
        except Exception:
            output_path = False

    except Exception as e:
        print('Error reading simulation parameters', e)
        exit(1)

# Calculation
try:
    orbit_count = simulation_parameters['Sim']['duration'] * 24 * 60 / spb.get_orbit_period(TLE)
    orbit_duration_txt = str(simulation_parameters['Sim']['duration'])+'d'
except Exception:
    orbit_duration_txt = f'{orbit_count:.1f}'
    pass
print('------------------', TLE[0], '------------------')
print('Orbit period:', spb.get_orbit_period(TLE), 'min')
print("Orbit count: ", orbit_count)
mean_p_coeff, mean_power, side_power_values, orbit_side_mean_power, eclipse_pc = \
    spb.power_calculation(TLE, JD_ini,
                            orbit_count,
                            ang_v, pose0, pv_eff,
                            pv_count, pv_area, timestep=simulation_parameters['Sim']['timestep'],
                            show_progress=True)
print("Eclipse percentage (min/max)",
        f'{min(eclipse_pc)*100:.1f}',
        f'{max(eclipse_pc)*100:.1f}')
# Mean power coefficient
print("Mean Power coefficient of each side")
for v, k in tuple(zip(spb.axis_labels, mean_p_coeff)):
    print(v, k)
print("Total Mean Power coefficient", sum(mean_p_coeff))
print()
# Mean power
print("Mean Power of each side")
for v, k in tuple(zip(spb.axis_labels, mean_power)):
    print(v, k)
print("Total Mean Power in mW", mean_power.sum(), 'mW')
# Plot power
if simulation_parameters['Sim']['stack_axis'] or len(skip_plot_axis) == 5:
    plot_stacked(plt, TLE[0] + " " + orbit_duration_txt, TLE,
                    side_power_values, skip_plot_axis)
else:
    plot_split(plt, TLE[0] + " " + orbit_duration_txt, TLE,
                side_power_values, skip_plot_axis)
if output_path:
    plt.savefig(os.path.join(output_path,
                TLE[0] + " " + orbit_duration_txt + ".svg"))
    plt.savefig(os.path.join(output_path,
                TLE[0] + " " + orbit_duration_txt + ".png"))
# Plot power per orbit
plot_stacked(plt,
                TLE[0] + " Orbit Mean " + orbit_duration_txt,
                TLE,
                orbit_side_mean_power,
                skip_plot_axis,
                xlabel='Orbit')
if output_path:
    plt.savefig(os.path.join(output_path,
                TLE[0] + " Orbit Mean " + orbit_duration_txt + ".svg"))
    plt.savefig(os.path.join(output_path,
                TLE[0] + " Orbit Mean " + orbit_duration_txt + ".png"))
plt.show()
