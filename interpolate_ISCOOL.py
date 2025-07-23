import os
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from typing import Tuple
from numpy.typing import ArrayLike

def filter_ISCOOL(path: str, mass_range: list, filename: str, ISCOOL_voltage_multiplier: float, plot: bool = False, 
				  iscool_columns: list = ['timestamp', 'offset', 'iscool_voltage'], method: str = 'spline_iscool', **args) -> interp1d:
	'''Main function to call for filtering ISCOOL. Two methods: linear spline or Savitzky-Golay filter

        Parameters
        ----------
        path: str
            path to the raw data
        mass_range: list
            range of masses that it should filter the ISCOOL values for
        filename: str
        	name of the iscool file, eg. iscool or iscool2
        ISCOOL_voltage_multiplier: float
        	factor to account for the voltage divider
        plot: bool
        	True if plot of ISCOOL vs timestamp in hours
        iscool_columns: list, default: ['timestamp', 'offset', 'iscool_voltage']
        	column names for the iscool .csv file
        method: str
        	method to use for the ISCOOL filter

        **args
        Only needed when method == 'savgol'
        See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html
            window_length: int
            polyorder: int
            deriv: int, optional, default: 0
            delta: float, optional, default: 1.0
            axis: int, optional, default: -1
            mode: str, optional, default: 'interp'
            cval: scalar, optional, default: 0.0
        Returns
        -------
        scipy.interpolate.interp1d
        '''
	if method.lower() == 'spline_iscool':
		return interpolateISCOOL(path = path, mass_range = mass_range, filename = filename, ISCOOL_voltage_multiplier = ISCOOL_voltage_multiplier, plot = plot)
	elif method.lower() == 'savgol':
		args.setdefault('deriv',0) 
		args.setdefault('delta',1) 
		args.setdefault('axis',-1) 
		args.setdefault('mode','interp') 
		args.setdefault('cval', 0.0)
		iscool,timestamp = get_array_iscool(path, mass_range, filename)
		iscool = iscool * ISCOOL_voltage_multiplier
		filtered_iscool = savgol_filter(iscool, window_length = args['window_length'], polyorder = args['polyorder'], deriv = args['deriv'], delta = args['delta'], 
			axis = args['axis'], mode = args['mode'], cval = args['cval'])
		if plot:
			fig,ax = plt.subplots(figsize = (14,9))
			ax.ticklabel_format(useOffset=False)
			ax.tick_params(axis='both', which='major', labelsize=20)
			plt.plot(timestamp/3600000,iscool,'-', label = 'Measured ISCOOL voltages')
			plt.plot(timestamp/3600000, filtered_iscool, label = 'Filtered ISCOOL voltage')
			plt.legend(fontsize = 20)
			plt.xlabel('timestamp [hours]', fontsize = 25)
			plt.ylabel('ISCOOL voltage', fontsize = 25)
			plt.show()
		filtered_ISCOOL_spline = interp1d(timestamp, filtered_iscool, bounds_error=False, fill_value='extrapolate')
		return filtered_ISCOOL_spline
	else:
		raise RuntimeError('Not yet implemented, implement the filter u want yourself, u lazy fk')

def interpolateISCOOL(path: str, mass_range: list, filename: str, ISCOOL_voltage_multiplier: float, plot: bool = False, iscool_columns: list = ['timestamp', 'offset', 'iscool_voltage']) -> interp1d:
	'''Creates a linear spline for ISCOOL interpolation. Adds one value at the beginning and end to not raise the input x is outside of interpolation range

        Parameters
        ----------
        path: str
            path to the raw data
        mass_range: list
            range of masses that it should filter the ISCOOL values for
        filename: str
        	name of the iscool file, eg. iscool or iscool2
        ISCOOL_voltage_multiplier: float
        	factor to account for the voltage divider
        plot: bool
        	True if plot of ISCOOL vs timestamp in hours
        iscool_columns: list, default: ['timestamp', 'offset', 'iscool_voltage']
        	column names for the iscool .csv file
        Returns
        -------
        scipy.interpolate.interp1d
        '''
	iscool,timestamp = get_array_iscool(path, mass_range, filename)
	iscool = iscool * ISCOOL_voltage_multiplier 
	iscool = np.append(iscool,iscool[-1]) # in case your first or last scan has no iscool, i just take the last recorded value and extrapolate it. Its not nice but it will have to do i guess
	timestamp = np.append(timestamp,timestamp[-1]+400000000) 
	ISCOOL_spline = interp1d(timestamp, iscool, bounds_error=False, fill_value='extrapolate')# take a spline, if you want to use the filter is sent. you will have to change this
	if plot: # for plotting
		fig,ax = plt.subplots(figsize = (14,9))
		test_x = np.arange(min(timestamp), max(timestamp),1000) 
		# since i take steps of 1000 s, the plot you see here will not be correct but if you take steps of 1 (or 5 s) you will get the correct plot but then you cannot see the iscool values anymore 
		test_y = ISCOOL_spline(test_x)
		ax.ticklabel_format(useOffset=False)
		ax.tick_params(axis='both', which='major', labelsize=20)
		plt.plot(timestamp/3600000,iscool,'-', label = 'Measured ISCOOL voltages')
		plt.plot(test_x/3600000, test_y, label = 'Interpolation ISCOOL voltage')
		plt.legend(fontsize = 20)
		plt.xlabel('timestamp [hours]', fontsize = 25)
		plt.ylabel('ISCOOL voltage', fontsize = 25)
		plt.show()
	return ISCOOL_spline

def get_array_iscool(path: str, mass_range: list, filename: str, iscool_columns: list = ['timestamp', 'offset', 'iscool_voltage']) -> Tuple[ArrayLike,ArrayLike]:
	'''Reads the timestamp and ISCOOL readout for a given mass range and returns in an array

        Parameters
        ----------
        path: str
            path to the raw data
        mass_range: list
            range of masses that it should filter the ISCOOL values for
        filename: str
        	name of the iscool file, eg. iscool or iscool2
        iscool_columns: list, default: ['timestamp', 'offset', 'iscool_voltage']
        	column names for the iscool .csv file
        Returns
        -------
        tuple of raw iscool readout and timestamp
        '''
	iscool = np.array([])
	timestamp = np.array([])
	for mass in os.listdir(path): #iterate over all masses
		try:
			mass = int(mass)
			# you can leave out the first one in the if statement because it was specifically for silver since in silver data there was one mass folder with mass 1094 
			if (mass >= min(mass_range)) and (mass <= max(mass_range)):
				mass_path = path + '\\' + str(mass) + '\\' # define path to the mass folder
				for scan in os.listdir(path + '\\' + str(mass)): #iterate over all scans
					iscool_path = f'{mass_path}{scan}\\{filename}_ds.csv' # define the path to the iscool file
					try:
						tempiscool = pd.read_csv(iscool_path, sep=';', header=None, names=iscool_columns) # read iscool
						iscool = np.concatenate((iscool, np.array(tempiscool['iscool_voltage']))) # concatenate it to an array with the rest of the iscool values
						timestamp = np.concatenate((timestamp,np.array(tempiscool['timestamp']))) # concatenate it to an array with the rest of the timestamp values
					except:
						pass
		except:
			continue
		else:
			continue
	array = np.array([timestamp,iscool])
	array = array[:,array[0].argsort()]
	iscool = array[1]
	timestamp = array[0]
	return iscool,timestamp