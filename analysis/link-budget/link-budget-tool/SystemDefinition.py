import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import yaml
import scipy.constants
import random
import itur

import pylink

from utils import Loader

class SystemDefinition:
    
    def __init__(self, satellite_configuration, groundstation_configuration, channel_configuration, monte_carlo=False, warn=True):
        
        # store the parameters internally
        self.satellite_configuration = satellite_configuration
        self.groundstation_configuration = groundstation_configuration
        self.channel_configuration = channel_configuration
        self.monte_carlo = monte_carlo
        self.warn = warn
        
        # if we're doing Monte-Carlo we need to choose an elevation to use
        if self.monte_carlo:
            
            distribution = channel_configuration['statistical']['geometry']['elevation']
            
            assert distribution[0] == 'Uniform'
            
            minimum_elevation =  float(distribution[1][0])
            maximum_elevation =  float(distribution[1][1])
            
            self.elevation = (maximum_elevation - minimum_elevation)*random.random() + minimum_elevation
                    
        else:
            self.elevation = channel_configuration['nominal']['geometry']['elevation']

        # set the components
        self.satellite_init()
        self.groundstation_init()
        self.channel_init()

        # defaults to DVB-S2X
        self.modulation = pylink.Modulation()

        # form the downlink signal chain
        self.downlink = pylink.DAGModel([
               self.geometry,
               self.gs_rx_antenna,
               self.sat_transmitter,
               self.sat_tx_antenna,
               self.gs_receiver,
               self.downlink_channel,
               self.rx_interconnect,
               self.tx_interconnect,
               self.modulation,
               pylink.LinkBudget(
                   name        = 'Downlink',
                   is_downlink = True)])

        # form the downlink signal chain
        self.uplink = pylink.DAGModel([
               self.geometry,
               self.sat_rx_antenna,
               self.sat_receiver,
               self.gs_transmitter,
               self.gs_tx_antenna,
               self.uplink_channel,
               pylink.LinkBudget(
                   name        = 'Uplink',
                   is_downlink = False)])
        
        
        
    def satellite_init(self):

        # build the receive RF chain
        self.sat_rf_chain = [
            pylink.Element(
               name            = 'Cables',
               gain_db         = self.sampler(self.satellite_configuration, ['rx', 'cable', 'gain']),
               noise_figure_db = self.sampler(self.satellite_configuration, ['rx', 'cable', 'noise_figure'])),
            pylink.Element(
               name            = 'LNA',
               gain_db         = self.sampler(self.satellite_configuration, ['rx', 'lna', 'gain']),
               noise_figure_db = self.sampler(self.satellite_configuration, ['rx', 'lna', 'noise_figure'])),
            pylink.Element(
               name            = 'Filter',
               gain_db         = self.sampler(self.satellite_configuration, ['rx', 'filter', 'gain']),
               noise_figure_db = self.sampler(self.satellite_configuration, ['rx', 'filter', 'noise_figure'])),
            pylink.Element(
               name            = 'Demodulator',
               gain_db         = self.sampler(self.satellite_configuration, ['rx', 'demodulator', 'gain']),
               noise_figure_db = self.sampler(self.satellite_configuration, ['rx', 'demodulator', 'noise_figure'])),
            ]

        # if no gain given, calculate from the antenna aperture
        if 'gain' in self.satellite_configuration['nominal']['rx']['antenna']:
            sat_rx_antenna_gain = self.sampler(self.satellite_configuration, ['rx', 'antenna', 'gain'])

            if 'aperture' in self.satellite_configuration['nominal']['rx']['antenna'] and self.warn:
                print("Satellite Rx Antenna: Both gain and aperature provided, will use gain")

        else:
            sat_rx_antenna_gain = 10*np.log10(
                                    self.sampler(self.satellite_configuration, ['rx', 'antenna', 'aperture_efficiency'])
                                    * self.sampler(self.satellite_configuration, ['rx', 'antenna', 'aperture'])
                                    * 4 * np.pi
                                    * (1e6 * self.sampler(self.channel_configuration, ['uplink', 'center_freq']))**2
                                    / scipy.constants.c**2)

        # specify the receive antenna
        self.sat_rx_antenna = pylink.Antenna(
               gain             = sat_rx_antenna_gain,
               polarization     = self.satellite_configuration['nominal']['rx']['antenna']['polarization'],
               pattern          = self.sampler(self.satellite_configuration, ['rx', 'antenna', 'pattern']),
               is_rx            = True,
               tracking         = self.satellite_configuration['nominal']['rx']['antenna']['tracking'],
               pointing_loss_db = self.sampler(self.satellite_configuration, ['rx', 'antenna', 'pointing_loss']),
               rx_noise_temp_k  = self.sampler(self.satellite_configuration, ['rx', 'antenna', 'noise_temp']))
        
        
        # if no gain given, calculate from the antenna aperture
        if 'gain' in self.satellite_configuration['nominal']['tx']['antenna']:
            sat_tx_antenna_gain = self.sampler(self.satellite_configuration, ['tx', 'antenna', 'gain'])

            if 'aperture' in self.satellite_configuration['nominal']['tx']['antenna'] and self.warn:
                print("Satellite Tx Antenna: Both gain and aperature provided, will use gain")

        else:
            sat_tx_antenna_gain = 10*np.log10(
                                    self.sampler(self.satellite_configuration, ['tx', 'antenna', 'aperture_efficiency'])
                                    * self.sampler(self.satellite_configuration, ['tx', 'antenna', 'aperture'])
                                    * 4 * np.pi
                                    * (1e6 * self.sampler(self.channel_configuration, ['downlink', 'center_freq']))**2
                                    / scipy.constants.c**2)

        self.sat_tx_antenna = pylink.Antenna(
               gain             = sat_tx_antenna_gain,
               polarization     = self.satellite_configuration['nominal']['tx']['antenna']['polarization'],
               pattern          = self.sampler(self.satellite_configuration, ['tx', 'antenna', 'pattern']),
               is_rx            = False,
               tracking         = self.satellite_configuration['nominal']['tx']['antenna']['tracking'],
               pointing_loss_db = self.sampler(self.satellite_configuration, ['rx', 'antenna', 'pointing_loss']))

        self.sat_receiver = pylink.Receiver(
               rf_chain        = self.sat_rf_chain,
               implementation_loss_db = self.sampler(self.satellite_configuration, ['rx', 'receiver', 'implementation_loss']),
               name = 'Satellite ' + self.channel_configuration['nominal']['uplink']['name'] + ' Receiver')

        self.sat_transmitter = pylink.Transmitter(
               tx_power_at_pa_dbw = self.sampler(self.satellite_configuration, ['tx', 'pa', 'power']),
               name = 'Satellite ' + self.channel_configuration['nominal']['downlink']['name'] + ' Transmitter')

        
    def groundstation_init(self):
        
        ## define the ground station
        self.gs_rf_chain = [
            pylink.Element(
               name            = 'Cables',
               gain_db         = self.sampler(self.groundstation_configuration, ['rx', 'cable', 'gain']),
               noise_figure_db = self.sampler(self.groundstation_configuration, ['rx', 'cable', 'noise_figure'])),
            pylink.Element(
               name            = 'LNA',
               gain_db         = self.sampler(self.groundstation_configuration, ['rx', 'lna', 'gain']),
               noise_figure_db = self.sampler(self.groundstation_configuration, ['rx', 'lna', 'noise_figure'])),
            pylink.Element(
               name            = 'Filter',
               gain_db         = self.sampler(self.groundstation_configuration, ['rx', 'filter', 'gain']),
               noise_figure_db = self.sampler(self.groundstation_configuration, ['rx', 'filter', 'noise_figure'])),
            pylink.Element(
               name            = 'Demodulator',
               gain_db         = self.sampler(self.groundstation_configuration, ['rx', 'demodulator', 'gain']),
               noise_figure_db = self.sampler(self.groundstation_configuration, ['rx', 'demodulator', 'noise_figure'])),
            ]
        
        # if no gain given, calculate from the antenna aperture
        if 'gain' in self.groundstation_configuration['nominal']['rx']['antenna']:
            self.gs_rx_antenna_gain = self.sampler(self.groundstation_configuration, ['rx', 'antenna', 'gain'])

            if 'aperture' in self.groundstation_configuration['nominal']['rx']['antenna'] and self.warn:
                print("Ground Station Rx Antenna: Both gain and aperature provided, will use gain")

        else:
            self.gs_rx_antenna_gain = 10*np.log10(
                                    self.sampler(self.groundstation_configuration, ['rx', 'antenna', 'aperture_efficiency'])
                                    * self.sampler(self.groundstation_configuration, ['rx', 'antenna', 'aperture'])
                                    * 4 * np.pi
                                    * (1e6 * self.sampler(self.channel_configuration, ['downlink', 'center_freq']))**2
                                    / scipy.constants.c**2)

        self.gs_rx_antenna = pylink.Antenna(
               gain             = self.gs_rx_antenna_gain,
               pattern          = self.sampler(self.groundstation_configuration, ['rx', 'antenna', 'pattern']),
               polarization     = self.groundstation_configuration['nominal']['rx']['antenna']['polarization'],
               is_rx            = True,
               tracking         = self.groundstation_configuration['nominal']['rx']['antenna']['tracking'],
               pointing_loss_db = self.sampler(self.groundstation_configuration, ['rx', 'antenna', 'pointing_loss']),
               rx_noise_temp_k  = self.sampler(self.groundstation_configuration, ['rx', 'antenna', 'noise_temp']))
        
        # if no gain given, calculate from the antenna aperture
        if 'gain' in self.groundstation_configuration['nominal']['tx']['antenna']:
            self.gs_tx_antenna_gain = self.sampler(self.groundstation_configuration, ['tx', 'antenna', 'gain'])

            if 'aperture' in self.groundstation_configuration['nominal']['tx']['antenna'] and self.warn:
                print("Ground Station Tx Antenna: Both gain and aperature provided, will use gain")
                
        else:
            self.gs_tx_antenna_gain = 10*np.log10(
                                    self.sampler(self.groundstation_configuration, ['tx', 'antenna', 'aperture_efficiency'])
                                    * self.sampler(self.groundstation_configuration, ['tx', 'antenna', 'aperture'])
                                    * 4 * np.pi
                                    * (1e6 * self.sampler(self.channel_configuration, ['uplink', 'center_freq']))**2
                                    / scipy.constants.c**2)

        self.gs_tx_antenna = pylink.Antenna(
               gain             = self.gs_tx_antenna_gain,
               pattern          = self.sampler(self.groundstation_configuration, ['tx', 'antenna', 'pattern']),
               polarization     = self.groundstation_configuration['nominal']['tx']['antenna']['polarization'],
               is_rx            = False,
               tracking         = self.groundstation_configuration['nominal']['tx']['antenna']['tracking'],
               pointing_loss_db = self.sampler(self.groundstation_configuration, ['tx', 'antenna', 'pointing_loss']))


        self.gs_receiver = pylink.Receiver(
               rf_chain        = self.gs_rf_chain,
               name            = 'Ground ' + self.channel_configuration['nominal']['downlink']['name'] + ' Receiver')

        self.gs_transmitter = pylink.Transmitter(
               tx_power_at_pa_dbw = self.sampler(self.groundstation_configuration, ['tx', 'pa', 'power']),
               name               = 'Ground ' + self.channel_configuration['nominal']['uplink']['name'] + ' Transmitter')
        
        
        
    def channel_init(self):
              
        self.rx_interconnect = pylink.Interconnect(is_rx=True)
        self.tx_interconnect = pylink.Interconnect(is_rx=False)

        altitude = self.sampler(self.channel_configuration, ['geometry', 'altitude'])
        self.geometry = pylink.Geometry(
               apoapsis_altitude_km   = altitude,
               periapsis_altitude_km  = altitude,
               min_elevation_deg      = self.elevation)

        self.downlink_channel = pylink.Channel(
               bitrate_hz                    = self.sampler(self.channel_configuration, ['downlink', 'bitrate']),
               allocation_hz                 = self.sampler(self.channel_configuration, ['downlink', 'allocation']),
               center_freq_mhz               = self.sampler(self.channel_configuration, ['downlink', 'center_freq']),
               atmospheric_loss_db           = self.sampler(self.channel_configuration, ['downlink', 'atmospheric_loss']),
               ionospheric_loss_db           = self.sampler(self.channel_configuration, ['downlink', 'ionospheric_loss']),
               rain_loss_db                  = self.sampler(self.channel_configuration, ['downlink', 'rain_loss']),
               multipath_fading_db           = self.sampler(self.channel_configuration, ['downlink', 'multipath_fading']),
               polarization_mismatch_loss_db = self.sampler(self.channel_configuration, ['downlink', 'polarization_mismatch_loss']))
        
        self.uplink_channel = pylink.Channel(
               bitrate_hz                    = self.sampler(self.channel_configuration, ['uplink', 'bitrate']),
               allocation_hz                 = self.sampler(self.channel_configuration, ['uplink', 'allocation']),
               center_freq_mhz               = self.sampler(self.channel_configuration, ['uplink', 'center_freq']),
               atmospheric_loss_db           = self.sampler(self.channel_configuration, ['uplink', 'atmospheric_loss']),
               ionospheric_loss_db           = self.sampler(self.channel_configuration, ['uplink', 'ionospheric_loss']),
               rain_loss_db                  = self.sampler(self.channel_configuration, ['uplink', 'rain_loss']),
               multipath_fading_db           = self.sampler(self.channel_configuration, ['uplink', 'multipath_fading']),
               polarization_mismatch_loss_db = self.sampler(self.channel_configuration, ['uplink', 'polarization_mismatch_loss']))
        
        
    def sampler(self, configuration, item_list):
        
        # run down the heirarchy of each dictonary item defined in item_list
        level_nominal = configuration['nominal']
        for item in item_list:
            level_nominal = level_nominal[item]
            
        # we've reached the bottom, turn to float
        try:
            if type(level_nominal) == list:
                mean = [float(_) for _ in level_nominal]
            elif level_nominal == 'ITU':  # should deal with this better
                mean = level_nominal
            else:
                mean = float(level_nominal)
        
            # no sampling performed, just return the mean
            if not self.monte_carlo:

                # check if we need to calculate the values
                if mean == 'ITU':

                    # pull out the info
                    frequency = float(configuration['nominal'][item_list[0]]['center_freq']) * 1e6
                    link = item_list[0]

                    # create a random number between 0 and 1
                    percentage = random.random()

                    # find the losses
                    rain_loss, atmospheric_loss = self.ITU_loss(frequency=frequency, link=link)

                    # split out the results
                    if item_list[-1] == 'atmospheric_loss':
                        return atmospheric_loss
                    elif item_list[-1] == 'rain_loss':
                        return rain_loss
                else:
                    return mean
            else:

                # run down the heirarchy of each dictonary item defined in item_list
                level_statistical = configuration['statistical']
                for item in item_list:
                    level_statistical = level_statistical[item]
                distribution = level_statistical

                # find out what distribution to use
                if distribution[0] == 'Uniform':
                    minimum_elevation =  float(distribution[1][0])
                    maximum_elevation =  float(distribution[1][1])
                    temp = (maximum_elevation - minimum_elevation)*random.random() + minimum_elevation                
                    return temp
                
                elif distribution[0] == 'Gaussian':
                    if type(mean) == list:
                        return [random.gauss(_, distribution[1]) for _ in mean]
                    else:
                        temp = random.gauss(mean, distribution[1])
                        return temp

                elif distribution[0] == 'ITU':

                    # pull out the info
                    location = distribution[1]
                    frequency = float(configuration['nominal'][item_list[0]]['center_freq']) * 1e6
                    link = item_list[0]

                    # create a random number between 0 and 1 - bound it on the low end
                    percentage = max(random.random(), 0.01)
                    
                    # find the losses
                    rain_loss, atmospheric_loss = self.ITU_loss(frequency=frequency, link=link, location=location, 
                                                                percentage=percentage)

                    # split out the results
                    if item_list[-1] == 'atmospheric_loss':
                        return atmospheric_loss
                    elif item_list[-1] == 'rain_loss':
                        return rain_loss

                elif distribution[0] == 'NotValid':
                    return mean
                elif distribution[0] == 'Fixed':
                    return mean
                else:
                    print("Statistical distribution '%s' not supported yet" % distibution[0])
                    return NaN
        except:
            return level_nominal
            
            
    def update_gs_rx_antenna_gain(self, gain):
        

        # create new antenna
        self.gs_rx_antenna = pylink.Antenna(
               gain             = gain,
               pattern          = self.sampler(self.groundstation_configuration, ['rx', 'antenna', 'pattern']),
               polarization     = self.groundstation_configuration['nominal']['rx']['antenna']['polarization'],
               is_rx            = True,
               tracking         = self.groundstation_configuration['nominal']['rx']['antenna']['tracking'],
               pointing_loss_db = self.sampler(self.groundstation_configuration, ['rx', 'antenna', 'pointing_loss']),
               rx_noise_temp_k  = self.sampler(self.groundstation_configuration, ['rx', 'antenna', 'noise_temp']))


        # form the downlink signal chain
        self.downlink = pylink.DAGModel([
               self.geometry,
               self.gs_rx_antenna,
               self.sat_transmitter,
               self.sat_tx_antenna,
               self.gs_receiver,
               self.downlink_channel,
               self.rx_interconnect,
               self.tx_interconnect,
               self.modulation,
               pylink.LinkBudget(
                   name        = 'Downlink',
                   is_downlink = True)])

        
    def set_elevation(self, elevation):
        
        self.downlink.override(self.downlink.enum.min_elevation_deg, elevation)
        self.elevation = elevation
        
        
                                          
    def ITU_loss(self, frequency, link, location=[41.39, -71.05], percentage=0.1):
            
        # calculate the matching dish diameter -  assume efficiency of 65%
        if link == 'downlink':
            effective_diameter = (scipy.constants.c/(np.pi*frequency)) * np.sqrt(10**(self.gs_rx_antenna_gain/10)/0.65) * itur.u.m 
        else:
            effective_diameter = (scipy.constants.c/(np.pi*frequency)) * np.sqrt(10**(self.gs_tx_antenna_gain/10)/0.65) * itur.u.m
        
        elevation = self.elevation
            
        # convert frequency to GHz
        frequency_ghz = (frequency / 1e9) * itur.u.GHz

        # Compute atmospheric parameters
        hs = itur.topographic_altitude(41.39, -71.05)
        T = itur.surface_mean_temperature(41.39, -71.05)
        P = itur.models.itu835.pressure(41.39, hs)
        rho_p = itur.surface_water_vapour_density(location[0], location[1], percentage, hs)
        rho_sa = itur.models.itu835.water_vapour_density(location[0], hs)
        T_sa = itur.models.itu835.temperature(location[0], hs)
        V = itur.models.itu836.total_water_vapour_content(location[0], location[1], percentage, hs)

        # Compute rain and cloud-related parameters
        R_prob = itur.models.itu618.rain_attenuation_probability(location[0], location[1], elevation, hs)
        R_pct_prob = itur.models.itu837.rainfall_probability(location[0], location[1])
        R001 = itur.models.itu837.rainfall_rate(location[0], location[1], percentage)
        h_0 = itur.models.itu839.isoterm_0(location[0], location[1])
        h_rain = itur.models.itu839.rain_height(location[0], location[1])
        L_red = itur.models.itu840.columnar_content_reduced_liquid(location[0], location[1], percentage)
        A_w = itur.models.itu676.zenit_water_vapour_attenuation(location[0], location[1], percentage, frequency_ghz, h=hs)

        # Compute attenuation values
        A_g = itur.gaseous_attenuation_slant_path(frequency_ghz, elevation, rho_p, P, T)
        A_r = itur.rain_attenuation(location[0], location[1], frequency_ghz, elevation, hs=hs, p=percentage)
        A_c = itur.cloud_attenuation(location[0], location[1], elevation, frequency_ghz, percentage)
        A_s = itur.scintillation_attenuation(location[0], location[1], frequency_ghz, elevation, percentage, effective_diameter)
        A_t = itur.atmospheric_attenuation_slant_path(location[0], location[1], frequency_ghz, elevation, percentage, effective_diameter)

        # rain attenuation plus all other attenuations
        rain = A_r
        atmosphere = A_g + A_c + A_s
        return float(rain.value), float(atmosphere.value)
    
    def print_link_info(self, link):
        
        # select the link to use
        if link == 'downlink':
            link = self.downlink
        else:
            link = self.uplink
            
        
        # print out the info
        print('Elevation (degress):            %3g'  % link.min_elevation_deg)
        print('Slant Range (km):               %3g'  % link.slant_range_km)
        print('Antenna Angle (deg):            %3g'  % link.satellite_antenna_angle_deg)
        print('Total Rx Noise Temperature (K): %3g'  % link.rx_noise_temp_k)
        print('Receiver Noise Temperature (K): %3g'  % link.rx_system_noise_temp_k)
        print('Noise Factor:                   %3g'  % link.rx_system_noise_factor)
        print('Noise Figure (dB):              %3g'  % link.rx_system_noise_figure)
        print('Transmit Power (dBW):           %3g'  % link.tx_power_at_antenna_dbw)
        print('Transmit EIRP (dBW):            %3g'  % link.tx_eirp_dbw)
        print('UGPL (dB):                      %3g'  % link.unity_gain_propagation_loss_db)
        print('Total Channel Loss (dB):        %3g'  % link.total_channel_loss_db)
        print('Link Margin (dB):               %3g'  % link.link_margin_db)
        print('Noise BW Loss (dB):             %-3g' % link.excess_noise_bandwidth_loss_db)
        print('C / N0 (dB):                    %-3g' % link.cn0_db)
        print('Rx Power (dBW):                 %-3g' % link.rx_power_dbw)
        print('Rx N0 (dBW/Hz):                 %-3g' % link.rx_n0_dbw_per_hz)
        print('Eb / N0 (dB):                   %-3g' % link.rx_ebn0_db)
        print('Eb (dB):                        %-3g' % link.rx_eb)
        print('Bitrate Hz (dB):                %-3g' % link.bitrate_dbhz)

        code = link.best_modulation_code
        print("=== Modulation: %s ===" % link.modulation_name)
        print("  Best Code:                    %s" % code.name)
        print("  Transmit Spectral Efficiency: %f" % code.tx_eff)
        print("  Receive Spectral Efficiency:  %f" % code.rx_eff)
        print("  Required Es/N0:               %f" % code.esn0_db)
        print("  Required Eb/N0:               %f" % code.ebn0_db)
        
        if link.link_margin_db > 0:
            bandwidth = link.allocation_end_hz - link.allocation_start_hz
            print("  Achieved Datarate (Mb/s):     %f" % (bandwidth*code.rx_eff/1e6))
        else:
            print("  Achieved Datarate (Mb/s):     LINK DOES NOT CLOSE!")
            
    def find_confidence_intervals(self, data, intervals=[0.25, 0.5, 0.75]):
        
        # sort the data from lowest to highest
        data.sort()
        length = len(data)
        
        # find the median
        confidence = [data[int(0.5*length)]]
        
        # loop through the other intervals
        for interval in intervals:
            offset = interval/2
            confidence.append(data[int((0.5-offset)*length)])
            confidence.append(data[int((0.5+offset)*length)])
                      
        return confidence


    def plot_confidence_intervals(self, data_x, data_y, info, step=False, negative_shade=False):
                      
        # create new figure
        plt.figure()

        # plot either step or normal
        if step:
            plt.step(data_x, data_y[:, 0], 'b', where='post')
            for i in range(int(len(data_y[0])/2)):      
                plt.fill_between(data_x, data_y[:, 2*i+1], data_y[:, 2*i+2], color='b', alpha=.05, step='post')
                  
        else:
            plt.plot(data_x, data_y[:, 0], 'b')
            for i in range(int(len(data_y[0])/2)):
                plt.fill_between(data_x, data_y[:, 2*i+1], data_y[:, 2*i+2], color='b', alpha=.05)
        
        # add red shading to highlight negative margin
        if negative_shade:
            min_value = min(data_y.flatten())
            max_value = max(data_y.flatten())
            plt.fill_between([-1.2*data_x[0], 1.2*data_x[-1]], [0,0], 2*min_value, facecolor="red", alpha=0.1)
            plt.xlim([data_x[0], data_x[-1]])
            plt.ylim([min_value, max_value*1.2])
        
        # add the plot info
        plt.suptitle(info[0])
        plt.xlabel(info[1])
        plt.ylabel(info[2])
        plt.title('Monte-Carlo Analysis', size=10)
        plt.grid()
        plt.show()
                      
                      
    def link_metrics(self):
                      
        rx_power = self.downlink.rx_power_dbw
        link_margin = self.downlink.link_margin_db

        if self.downlink.link_margin_db > 0:
            bandwidth = self.downlink.allocation_end_hz - self.downlink.allocation_start_hz
            datarate = self.downlink.best_modulation_code.rx_eff * bandwidth
        else:
            datarate = 0
        
        return rx_power, link_margin, datarate