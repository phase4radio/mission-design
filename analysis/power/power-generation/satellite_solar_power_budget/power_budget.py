import solar_power_budget as spb
from sgp4.ext import jday
import numpy as np
import matplotlib.pyplot as plt
import yaml
import os.path

# Parameters
TLEs = []
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


with open('parameters.yaml', mode='r') as file:
    try:
        parameters = yaml.load(file, Loader=yaml.FullLoader)
        for name, spacecraft in parameters['TLEs'].items():
            # print(name, spacecraft)
            spacecraft['TLE'].insert(0, name)
            TLEs.append(tuple(spacecraft['TLE']))
        # Initialize angular velocities in rad/min [X, Y, Z]
        ang_v = np.array([
            parameters['Spacecraft']['ang_v']['X'],
            parameters['Spacecraft']['ang_v']['Y'],
            parameters['Spacecraft']['ang_v']['Z']
        ])
        # Initial angle conditions in deg [X, Y, Z]
        pose0 = np.deg2rad([
            parameters['Spacecraft']['pose0']['X'],
            parameters['Spacecraft']['pose0']['Y'],
            parameters['Spacecraft']['pose0']['Z']
        ])
        # PV efficient in each side [X+, X-, Y+, Y-, Z+, Z-]
        pv_eff = np.array([
            parameters['SolarPanels']['efficiency']['X+'],
            parameters['SolarPanels']['efficiency']['X-'],
            parameters['SolarPanels']['efficiency']['Y+'],
            parameters['SolarPanels']['efficiency']['Y-'],
            parameters['SolarPanels']['efficiency']['Z+'],
            parameters['SolarPanels']['efficiency']['Z-']
        ])
        # Number of PV in each side [X+, X-, Y+, Y-, Z+, Z-]
        pv_count = np.array([
            parameters['SolarPanels']['count']['X+'],
            parameters['SolarPanels']['count']['X-'],
            parameters['SolarPanels']['count']['Y+'],
            parameters['SolarPanels']['count']['Y-'],
            parameters['SolarPanels']['count']['Z+'],
            parameters['SolarPanels']['count']['Z-']
        ])
        # Active area of each PV in mm^2, [X+, X-, Y+, Y-, Z+, Z-]
        pv_area = np.array([
            parameters['SolarPanels']['area']['X+'],
            parameters['SolarPanels']['area']['X-'],
            parameters['SolarPanels']['area']['Y+'],
            parameters['SolarPanels']['area']['Y-'],
            parameters['SolarPanels']['area']['Z+'],
            parameters['SolarPanels']['area']['Z-']
        ])
        # Set the start day for simulation, Julian day
        JD_ini = jday(parameters['Sim']['JD_start']['year'],
                      parameters['Sim']['JD_start']['month'],
                      parameters['Sim']['JD_start']['day'],
                      parameters['Sim']['JD_start']['hour'],
                      parameters['Sim']['JD_start']['min'],
                      parameters['Sim']['JD_start']['second'])
        orbit_count = parameters['Sim']['orbits']
        for axis, skip in enumerate(parameters['Sim']['skip_plot_axis']):
            if (skip):
                skip_plot_axis.append(axis)
        try:
            output_path = os.path.expanduser(parameters['Sim']['output_path'])
        except Exception:
            output_path = False

    except Exception as e:
        print('Error reading parameters', e)
        exit(1)

# Calculation
for TLE in TLEs:
    try:
        orbit_count = parameters['Sim']['duration'] * 24 * 60 / spb.get_orbit_period(TLE)
        orbit_duration_txt = str(parameters['Sim']['duration'])+'d'
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
                              pv_count, pv_area, timestep=parameters['Sim']['timestep'],
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
    if parameters['Sim']['stack_axis'] or len(skip_plot_axis) == 5:
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
