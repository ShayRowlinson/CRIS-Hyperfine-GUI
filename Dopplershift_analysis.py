import numpy as np
import os
import scipy.constants as csts
from uncertainties import ufloat as uf
import uncertainties.unumpy as unumpy
from scipy.stats import chi2
from scipy.signal import find_peaks
from typing import Union
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
import pandas as pd 
from interpolate_ISCOOL import filter_ISCOOL

# change the some of the class variables to what you need. This is just for silver
class Dopplershift_data:

    _signal_ch = 1
    _noise_ch = 3
    _ref_channel = 'wavenumber_2' 
    _reference = 'diode' # 'diode' or 'hene'
    diode_wavenumber = 12578.84915 # This needs to be checked! #12816.470353 # /cm
    hene_wavenumber = 15798         # /cm
    amu = csts.physical_constants['atomic mass unit-electron volt relationship'][0] #eV/c^2
    tof_lower = 0       # microseconds
    tof_upper = 1000    # microseconds
    tof_binsize = 0.1   # microseconds
    _devices = {"tagger": ['timestamp', 'offset', 'bunch_no', 'events_per_bunch', 'channel', 'delta_t'], 
                "wavemeter": ['timestamp', 'offset', 'wavenumber_1', 'wavenumber_2', 'wavenumber_3', 'wavenumber_4', 'wavenumber_mixed'], # for silver there was an extra mixed channel, leave it out for other analysis
                "wavemeter_pdl": ['timestamp', 'offset', 'wavenumber_pdl_1'],
                "iscool2": ['timestamp', 'offset', 'iscool_voltage'],
                "diodes": ['timestamp', 'offset', 'photodiode_1', 'photodiode_2', 'photodiode_3', 'photodiode_4'],
                "powermeter_1": ['timestamp', 'offset', 'powermeter_1'],
                "powermeter_3": ['timestamp', 'offset', 'powermeter_3'],
                "cec_voltage": ['timestamp', 'offset', 'cec_measured_voltage', 'cec_set_voltage']}
    # dict that shows which data has been loaded already
    _loaded_devices = {}
    # you can change all of these in this file or after you initialise the class in main.py by using eg. Dopplershift_data.bin_size_voltage = 4
    bin_size_voltage = 2
    bin_size_MHz = 30
    _bin_size_cm = bin_size_MHz*0.0000334
    min_TOFgate_width = 10 #microseconds (in case the Q switch of the non res in TOF is higher than the actual bunch, the Q switch gives a very narrow TOF so this is the minimum width the TOF gate should be)
    _transition_wavenumber = 45945.340 # in cm^-1
    # mass range over which you want to interpolate ISCOOL over (tbh you can choose 1, 100000 if you want, its more in the case if you want to omit a certain mass)
    mass_range = [1,99999999] 
    _masses = {}

    def __init__(self, mass: int, scan: int, voltage_scanning: bool, wn_channel: str,
             wn_bounds: list = [0, 1e15], exact_mass: float = None,
             data_folder: str = None, save_data_path: str = None, save_fig_path: str = None) -> None:
        """
        Initialises the class to Doppler-shift CRIS data

        Parameters
        ----------
        mass : int
            Mass number of the isotope (e.g., 117 for Sb-117)
        scan : int
            Scan number to process
        voltage_scanning : bool
            Whether the data was taken using voltage scanning
        wn_channel : str
            Spectroscopy wavenumber channel to use (e.g., 'wavenumber_1')
        wn_bounds : list, optional
            Bounds to gate the wavenumber on, by default [0, 1e15]
        exact_mass : float, optional
            Exact mass of the isotope (in amu), required if not using internal mass table
        data_folder : str
            Path to raw data folder, should contain mass/scan_xxxx/
        save_data_path : str
            Directory to save binned `.csv` files
        save_fig_path : str
            Directory to save `.png` figures
        """
        self._mass = mass
        self._exact_mass = exact_mass if exact_mass is not None else self._masses[mass]
        self._scan = scan
        self._voltage_scanning = voltage_scanning
        self._wn_channel = wn_channel
        self._wn_bounds = wn_bounds
        print(exact_mass)

        if not all([data_folder, save_data_path, save_fig_path]):
            raise ValueError("All path arguments (data_folder, save_data_path, save_fig_path) must be supplied.")

        # Define internal paths for this scan
        self._PATH = os.path.join(data_folder, str(mass), f"scan_{scan}") + os.sep
        self._SAVEPATH = os.path.join(save_data_path, str(mass)) + os.sep
        self._SAVEPATH_FIG = os.path.join(save_fig_path, str(mass)) + os.sep

        # Create save directories if they don't exist
        os.makedirs(self._SAVEPATH, exist_ok=True)
        os.makedirs(self._SAVEPATH_FIG, exist_ok=True)
    
    def extract_raw_data(self, devices_to_read: dict, path: str) -> pd.DataFrame:
        '''Reads and saves all data from the given devices to a Pandas DataFrame indexed and sorted by timestamp

        Parameters
        ----------
        devices_to_read: dict
            dictionary containing each device that needs to be read with as value the column names
        path: str
            path to the data

        Returns
        -------
        pd.DataFrame
        '''
        devices_to_concat = []
        for device in devices_to_read.keys():
            device_path = path + device + '_ds.csv' 
            try:
                temp_data = pd.read_csv(device_path, sep = ';', header = None, names = devices_to_read[device])
                temp_data.drop(labels = ['offset'], axis = 1, inplace = True)
                devices_to_concat.append(temp_data)
                self._loaded_devices[device] = True
            except:
                self._loaded_devices[device] = False
                print(device + ' file is missing')
        try:
            data = pd.concat(devices_to_concat)
        except:
            raise RuntimeError('Either the mass/scan combination is wrong or the path is wrong')
        data['voltage'] = np.zeros(len(data['timestamp']))
        data['timestamp_copy'] = data['timestamp']
        data = data.set_index(['timestamp'])
        data = data.sort_index()

        # fill blank values in independent parameters and drop lines with NaN values
        # If still using the old ISCOOL device (not iscool2):
        #   do not use the ISCOOL from this since it is read out every 5 sec and this code just forwardfills the previous ISCOOL value, so you get a step function every 5 sec
        # just use the interpolation, filter or average
        for device_values in devices_to_read.values():
            for column in device_values:
                if column not in ['bunch_no','events_per_bunch','delta_t']: # do not fill in these values since then you make up data
                    try:
                        data[column] = data[column].ffill()
                        data[column] = data[column].bfill()
                    except:
                        pass
        data = data.dropna()
        if self._loaded_devices['tagger']:
            data['delta_t'] = data['delta_t'] / 2000 # convert delta_t to microseconds
        else:
            raise print('No tagger loaded')
        return data

    def AdvCutNoise(self, data: pd.DataFrame, threshold: float) -> pd.DataFrame:
        '''Advanced method to cut noise. BEWARE this can take quite some time
        The code evaluates every bunch that has a count in the magneTOF signal channel
        In each of these bunches it flags each row of signal if its delta_t is within range of a interval [dT-threshold, dT+threshold] with dT any delta_t of a count in the noise channel in that bunch
        If the bunch has some rows (but not all) that are flagged as noise, then they will be discarded
        Else all rows that are flagged as noise will be converted to empty bunches
        Then all data not recorded on the signal or empty bunches channel is discarded

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        threshold: float
            the delta_t threshold that determines the maximum time difference to differentiate two counts in the signal and noise channel
        Returns
        -------
        pd.DataFrame
        '''
        data['noise'] = False
        signal_data = data[(data['channel'] == self._signal_ch)]
        noise_data = data[(data['channel'] == self._noise_ch)]
        for timestamp in signal_data['timestamp_copy']:
            if timestamp in noise_data['timestamp_copy']:
                try:
                    diff = abs(noise_data['delta_t'][timestamp] - data['delta_t'][timestamp]).values.reshape(len(noise_data['delta_t'][timestamp]), len(data['delta_t'][timestamp]))
                    indices = np.argwhere(diff<threshold)
                    indices = indices.T[1]
                except:
                    diff = abs(noise_data['delta_t'][timestamp] - data['delta_t'][timestamp]).values
                    indices = np.argwhere(diff<threshold)
                    indices = indices.T[0]
                data['noise'][timestamp].iloc[indices] = True
                if not np.all(data['noise'][timestamp]) and np.any(data['noise'][timestamp]):
                    noise_bool = (data['noise'][timestamp].values)
                    indices_noise = np.argwhere(noise_bool).T[0]
                    ix = data.index == timestamp
                    rows = data.loc[ix].iloc[indices_noise].index
                    data.loc[rows, 'channel'] = 5
                    data.loc[rows, 'noise'] = False
        data.loc[(data['channel'] == self._noise_ch),'noise'] = True
        data.loc[(data['noise'] == True),'channel'] = -1
        data.loc[(data['noise'] == True),'delta_t'] = -0.005
        s = np.logical_or.reduce([
                                data['channel'] == -1, 
                                data['channel'] == self._signal_ch,
                                ])
        data = data[s]
        return data

    def SimpleCutNoise(self, data: pd.DataFrame) -> pd.DataFrame:
        '''Simple method to cut noise. The code cuts all data that is not in the signal or empty bunches channel

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        Returns
        -------
        pd.DataFrame
        '''
        s = np.logical_or.reduce([
                                data['channel'] == -1, 
                                data['channel'] == self._signal_ch,
                                ])
        data = data[s]
        return data

    def gate_tof(self, data: pd.DataFrame, manual: Union[list,bool] = False) -> pd.DataFrame:
        '''TOF gate function to automatically detect a TOF peak and gate on this peak, or gate on a manual TOF interval.

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        manual: list or bool, default: False
            interval for manual TOF gate, or False if you dont want to manually specify the gate
        Returns
        -------
        pd.DataFrame
        '''
        tof_bins = np.arange(self.tof_lower, self.tof_upper+self.tof_binsize, self.tof_binsize)
        tof, tof_edges = np.histogram(a = data['delta_t'], bins = tof_bins)
        tof_centers = tof_edges[:-1] + np.diff(tof_edges)/2
        TOF_spectrum = np.array([tof_centers,tof])
        if manual:
            tfig, tax = plt.subplots()#figsize = (14,9))
            plt.plot(tof_centers, tof, label = 'Counts per bunch', color = 'k')
            plt.axvspan(manual[0], manual[1], color = '#0066FF', alpha = 0.5, label = 'Time gate')
            plt.xlim([0,manual[1] + 50])            
            plt.xlabel(r'Time since bunch release ($\mu$s)', fontsize = 14)
            tax.tick_params(axis='both', which='major', labelsize=14)
            plt.ylabel('Counts', fontsize = 14)
            plt.legend(fontsize = 14)
            plt.title(str(self._mass)+"Sb: Scan "+str(self._scan))
            #plt.show()
            data.loc[((data['delta_t'] < manual[0]) | (data['delta_t'] > manual[1])), 'delta_t'] = -0.005
            data.loc[((data['delta_t'] < manual[0]) | (data['delta_t'] > manual[1])), 'channel'] = -1
            return data.loc[(data['delta_t'] < 0) | ( (data['delta_t'] > manual[0]) & (data['delta_t'] < manual[1]))]
        indexmax = int(np.where(TOF_spectrum[1] == max(tof))[0][0])
        indexleftbound = np.array(np.where(TOF_spectrum[1][0:indexmax] < (1-0.95) * max(tof)))[0][-1]
        indexrightbound = np.array(np.where(TOF_spectrum[1][indexmax:] < (1-0.95) * max(tof)))[0][0] + indexmax
        if indexrightbound - indexleftbound < self.min_TOFgate_width: 
            raise RuntimeError('TOF gate too narrow')
        tfig, tax = plt.subplots()#figsize = (14,9))
        plt.plot(tof_centers, tof, label = 'Counts per bunch', color = 'k')
        plt.axvspan(TOF_spectrum[0][indexleftbound], TOF_spectrum[0][indexrightbound], color = '#0066FF', alpha = 0.5, label = 'Time gate')
        plt.xlim([0,TOF_spectrum[0][indexrightbound] + 50])
        plt.xlabel(r'Time since bunch release ($\mu$s)', fontsize = 14)
        tax.tick_params(axis='both', which='major', labelsize=14)
        plt.ylabel('Counts', fontsize = 14)
        plt.legend(fontsize = 14)
        plt.title(str(self._mass)+"Sb: Scan "+str(self._scan))
        plt.savefig(f'{self._SAVEPATH_FIG}{self._scan}_tof.png', dpi = 400, bbox_inches='tight')
        #plt.show()
        data.loc[((data['delta_t'] < TOF_spectrum[0][indexleftbound]) | (data['delta_t'] > TOF_spectrum[0][indexrightbound])), 'delta_t'] = -0.005
        data.loc[((data['delta_t'] < TOF_spectrum[0][indexleftbound]) | (data['delta_t'] > TOF_spectrum[0][indexrightbound])), 'channel'] = -1
        return data.loc[(data['delta_t'] < 0) | ( (data['delta_t'] > TOF_spectrum[0][indexleftbound]) & (data['delta_t'] < TOF_spectrum[0][indexrightbound]))]

    def filter_scatter(self, data: pd.DataFrame, filename: str, method: str, ISCOOL_voltage_multiplier: float, **args) -> pd.DataFrame:
        '''Filtering of ISCOOL readout. Three methods are possible: average, linear spline or Savitzky-Golay filter

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        filename: str
            name of the ISCOOL file to read
        method: str
            name of method to use
        ISCOOL_voltage_multiplier: float
            factor to multiply the readout voltage with to account for a voltage divider
        
        **args
        Only needed when method == 'savgol'
        See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html
            window_length: int
            polyorder: int
            deriv: int, optional
            delta: float, optional
            axis: int, optional
            mode: str, optional
            cval: scalar, optional 
        Returns
        -------
        pd.DataFrame
        '''
        if method.lower() == 'avg':
             data.loc[data['iscool_voltage'] > 0, 'iscool_voltage'] = np.average(data['iscool_voltage'])*ISCOOL_voltage_multiplier
        elif method.lower() == 'spline_iscool':
            data.loc[data['iscool_voltage'] > 0, 'iscool_voltage'] = filter_ISCOOL(self._DATA_FOLDER[:-1], self.mass_range, filename = filename, ISCOOL_voltage_multiplier = ISCOOL_voltage_multiplier)(data['timestamp_copy'])
        elif method.lower() == 'savgol': 
            data.loc[data['iscool_voltage'] > 0, 'iscool_voltage'] = filter_ISCOOL(self._DATA_FOLDER[:-1], self.mass_range, filename = filename, ISCOOL_voltage_multiplier = ISCOOL_voltage_multiplier, method = 'savgol', **args)(data['timestamp_copy'])
        else:
            raise RuntimeError('Not yet implemented, implement the filter u want urself, u lazy fk')
        return data

    def calibrate_CRIS_voltagemeter(self, data: pd.DataFrame, calibration_factor: float = 1) -> pd.DataFrame:
        '''Calibrate the voltages measured with a calibration factor applied to the CEC voltage and save total voltage applied to ion bunch in column 'voltage'

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        calibration_factor: float, default: 1
            Factor to multiply the CEC readout
        Returns
        -------
        pd.DataFrame
        '''
        if self._voltage_scanning:
            try:
                data.loc[data['cec_measured_voltage'] > -999999, 'cec_measured_voltage'] = data['cec_measured_voltage']*1000*calibration_factor #because voltage divider
                data.loc[data['voltage'] > -999999, 'voltage'] = data['iscool_voltage'] - data['cec_measured_voltage']
            except:
                raise RuntimeError('No cec_voltage file loaded, change voltage_scanning to False')
        else:
            data.loc[data['voltage'] > -999999, 'voltage'] = data['iscool_voltage']
        return data

    def apply_wavenumber_correction(self, data: pd.DataFrame) -> pd.DataFrame:
        '''Apply wavenumber correction on each wavenumber channel with HeNe reference or diode reference.
        Prints warning message if correction is larger than 100 1/cm or if the reference channel is not valid

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        Returns
        -------
        pd.DataFrame
        '''
        if self._reference.lower() == 'diode':
            correction = data[self._ref_channel] - self.diode_wavenumber
        elif self._reference.lower() == 'hene':
            correction = data[self._ref_channel] - self.hene_wavenumber
        else:
            correction = 0
            print('you did not choose a valid reference option')
        if np.all(correction > 100):
            print('DIODE CORRECTION BIGGER THAN 100 /CM. Reference diode wavenumber might be off')
        for name in self._devices['wavemeter'][2:]: # [2:] because we do not use timestamp or offset here
            if name != self._ref_channel:
                data.loc[data[name] > -999999, name] = data[name] - correction
        return data

    def gate_wavenumber(self, data: pd.DataFrame, wavenumber: str) -> pd.DataFrame:
        '''Apply gate on supplied wavenumber
        If the supplied channel is mixed, then it takes the sum of wavenumber_3 and wavenumber_4 since the mixed readout can have a ~100-ish MHz offset

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        wavenumber: str
            name of spectroscopy wavenumber channel
        Returns
        -------
        pd.DataFrame
        '''
        if wavenumber.lower() == 'wavenumber_mixed':
            data[wavenumber] = data['wavenumber_4'] + data['wavenumber_3']
        try:
            return data[(data[wavenumber] > self._wn_bounds[0]) & (data[wavenumber] < self._wn_bounds[1])]
        except:
            raise print('No valid wavenumber channel was given')

    def dopplershift_wn(self, data: pd.DataFrame, freq_multiplier: int) -> ArrayLike:
        '''Frequency multiplies and dopplershifts the spectroscopy channel

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        freq_multiplier: int
            frequency multiplier to account for doubling or tripling of laser frequency
        Returns
        -------
        ArrayLike
        '''
        try:
            wn = data[self._wn_channel].to_numpy() * freq_multiplier
        except:
            raise print('Wavenumber is not loaded or not valid')
        voltage = data['voltage'].to_numpy()
        beta_voltage = np.sqrt(1 - ( ((self._exact_mass * self.amu)**2) / ( (voltage+(self._exact_mass*self.amu))*(voltage+(self._exact_mass*self.amu)) ) ))
        return wn * (1 - beta_voltage) / np.sqrt(1 - (beta_voltage*beta_voltage)) 

    def bin_data(self, data: pd.DataFrame, freq_multiplier: int) -> pd.DataFrame:
        '''Bin data with binning function for voltage or laser scanning

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        freq_multiplier: int
            frequency multiplier to account for doubling or tripling of laser frequency
        Returns
        -------
        pd.DataFrame
        '''
        if not self._loaded_devices['cec_voltage'] or not self._voltage_scanning:
            return self.bin_wm_data(data,freq_multiplier)
        else:
            return self.bin_voltage_data(data, freq_multiplier)

    def bin_voltage_data(self, data: pd.DataFrame, freq_multiplier: int) -> pd.DataFrame:
        '''Bin data with binning function for voltage scanning
        First the bins are initialised with each bin centered around data in bin
        Then the data is grouped in each bin by voltage.
        Then the data is dopplershifted
        In each bin the x and xerr is converted to MHz and transformed to relative frequency wrt the transition frequency 
        the counts per bin are calculated by counting the amount of rows minus the amount of rows that are empty bunches
        yerr is the sqrt of this
        The bunches are calculated by counting the amount of unique bunch numbers in each bin
        The binned data is sorted wrt x
        for completely empty bins the yerr is approximated by 1/sqrt(amount of bunches in that bin)

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        freq_multiplier: int
            frequency multiplier to account for doubling or tripling of laser frequency
        Returns
        -------
        pd.DataFrame
        '''
        voltage = (data['voltage']).to_numpy()
        voltage_bins = np.arange(voltage.min()-self.bin_size_voltage/2, voltage.max()+self.bin_size_voltage/2, self.bin_size_voltage)
        data['digit_index'] = np.digitize(voltage, voltage_bins) 
        data[self._wn_channel] = self.dopplershift_wn(data, freq_multiplier)
        groups = data.groupby('digit_index')
        df = pd.DataFrame()
        df[['x','xerr']] = groups[self._wn_channel].agg(['mean','std'])
        df['xerr'] = df['xerr'].fillna(0)
        df['x'] = (df['x'] - self._transition_wavenumber) * 29979.2458 
        df['xerr'] = df['xerr'] * 29979.2458
        df['y'] = (groups.count()['timestamp_copy'] - groups['delta_t'].agg(lambda x: x.le(0).sum())) # minus the empty bunches, doesnt matter if you take timestamp or smt else
        df['yerr'] = np.sqrt(df['y'])
        df['bunches'] = ([len(array) for array in list(groups['timestamp_copy'].unique())])
        df = df.sort_values(by = ['x']) # put it in dataframe and sort by x for convenience but it is not needed at all
        yerr_h, yerr_l = (self.poisson_interval_high(df['y'])-df['y']),(df['y'] - self.poisson_interval_low(df['y']))
        df['yerr'] = yerr_h
        df['median_scan_time'] = np.median(data['timestamp_copy'])
        return df

    def bin_wm_data(self, data: pd.DataFrame, freq_multiplier: int) -> pd.DataFrame:
        '''Bin data with binning function for laser scanning
        first the data is dopplershifted
        then the bins are initialised 
        the data is grouped in each bin by laser frequency in the atom frame.
        In each bin the x and xerr is converted to MHz and transformed to relative frequency wrt the transition frequency 
        the counts per bin are calculated by counting the amount of rows minus the amount of rows that are empty bunches
        yerr is the sqrt of this
        The bunches are calculated by counting the amount of unique bunch numbers in each bin
        The binned data is sorted wrt x
        for completely empty bins the yerr is approximated by 1/sqrt(amount of bunches in that bin)

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        freq_multiplier: int
            frequency multiplier to account for doubling or tripling of laser frequency
        Returns
        -------
        pd.DataFrame
        '''
        data[self._wn_channel] = self.dopplershift_wn(data, freq_multiplier)
        spectrum_bins = np.arange(data[self._wn_channel].min(), data[self._wn_channel].max() + self.bin_size_MHz*0.0000334, self.bin_size_MHz*0.0000334)
        data['digit_index'] = np.digitize(data[self._wn_channel], spectrum_bins) 
        groups = data.groupby('digit_index')
        df = pd.DataFrame()
        df[['x','xerr']] = groups[self._wn_channel].agg(['mean','std'])
        df['xerr'] = df['xerr'].fillna(0)
        df['x'] = (df['x'] - self._transition_wavenumber) * 29979.2458
        df['xerr'] = df['xerr'] * 29979.2458
        df['y'] = (groups.count()['timestamp_copy'] - groups['delta_t'].agg(lambda x: x.le(0).sum())) # minus the empty bunches, doesnt matter if you take timestamp or smt else
        df['yerr'] = np.sqrt(df['y'])
        df['bunches'] = ([len(array) for array in list(groups['bunch_no'].unique())])
        df = df.sort_values(by = ['x']) # put it in dataframe and sort by x for convenience but it is not needed at all
        yerr_h, yerr_l = (self.poisson_interval_high(df['y'])-df['y']),(df['y'] - self.poisson_interval_low(df['y']))
        df['yerr'] = yerr_h
        df['median_scan_time'] = np.median(data['timestamp_copy'])
        return df

    def bin_time_data(self, data: pd.DataFrame, binsize: Union[float,int]) -> pd.DataFrame:
        '''Bin time data for decay spectroscopy
        first the bins are initialised with the given binsize
        the data is grouped in each bin by time
        In each bin the x and xerr is calculated
        the counts per bin are calculated by counting the amount of rows minus the amount of rows that are empty bunches
        yerr is the sqrt of this
        The bunches are calculated by counting the amount of unique bunch numbers in each bin
        The binned data is sorted wrt x
        for completely empty bins the yerr is approximated by the biggest error given by the chi2 distribution

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        binsize: float, int
            binsize in seconds
        Returns
        -------
        pd.DataFrame
        '''
        data['timestamp_copy'] = data['timestamp_copy']
        time_bins = np.arange(np.min(data['timestamp_copy']), np.max(data['timestamp_copy']),binsize)
        data['digit_index'] = np.digitize(data['timestamp_copy'], time_bins) 
        groups = data.groupby('digit_index')
        df = pd.DataFrame()
        df[['x', 'xerr']] = groups['timestamp_copy'].agg(['mean','std'])
        df['xerr'] = df['xerr'].fillna(0)
        df['x'] = df['x'] - np.min(df['x'])
        df['y'] = (groups.count()['timestamp_copy'] - groups['delta_t'].agg(lambda x: x.le(0).sum())) 
        df['yerr'] = np.sqrt(df['y'])
        df['bunches'] = ([len(array) for array in list(groups['bunch_no'].unique())])
        df = df.sort_values(by = ['x']) # put it in dataframe and sort by x for convenience but it is not needed at all
        yerr_h, yerr_l = (self.poisson_interval_high(df['y'])-df['y']),(df['y'] - self.poisson_interval_low(df['y']))
        df['yerr'] = yerr_h
        df['median_scan_time'] = np.median(data['timestamp_copy'])
        return df

    def cut_data(self, data: pd.DataFrame, borders: ArrayLike):
        '''Specialised function to cut implantation away from decay spectroscopy data

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        borders: ArrayLike
            Array of all the borders for the intervals
        Returns
        -------
        pd.DataFrame
        '''
        concat_data = []
        for i,border in enumerate(borders):
            if i+1 < len(borders) and not i % 2:
                temp_data = data.loc[(data['x'] >= borders[i]) & (data['x'] < borders[i+1])]
                temp_data['x'] = temp_data['x'] - np.min(temp_data['x'])
                concat_data.append(temp_data)
        df = pd.concat(concat_data)
        return df

    def rebin(self, data: pd.DataFrame, binsize: Union[float,int]) -> pd.DataFrame:
        '''Specialised function to rebin decay spectroscopy data

        Parameters
        ----------
        data: pd.Dataframe
            all the data for this scan
        binsize: float, int
            binsize in seconds
        Returns
        -------
        pd.DataFrame
        '''
        time_bins = np.arange(np.min(data['x']), np.max(data['x']),binsize)
        data['digit'] = np.digitize(data['x'], time_bins) 
        groups = data.groupby('digit')
        df = pd.DataFrame()
        df[['x', 'xerr']] = groups['x'].agg([np.mean,np.std])
        df['xerr'] = df['xerr'].fillna(0)
        df['y'] = groups.sum()['y']
        df['yerr'] = np.sqrt(df['y'])
        df['bunches'] = ([len(array) for array in list(groups['y'])])
        df = df.sort_values(by = ['x']) # put it in dataframe and sort by x for convenience but it is not needed at all
        yerr_h, yerr_l = (self.poisson_interval_high(df['y'])-df['y']),(df['y'] - self.poisson_interval_low(df['y']))
        df['yerr'] = yerr_h
        df['median_scan_time'] = np.median(data['timestamp_copy'])
        return df

    def plot_time(self, data: pd.DataFrame, fig: plt.figure, ax: Union[plt.Axes, ArrayLike], save: bool = False, save_format: str = 'png', **kwargs) -> plt.figure:
        '''plot and save figure in given format for time data

        Parameters
        ----------
        data: pd.Dataframe
            Analysed data for this scanm
        fig: plt.figure
            figure object to plot in
        ax: plt.Axes
            axes object for each subplot
        save: bool, default: False
            True if the figure needs to be saved
        save_format: str, default: 'png'
            format for the saved figure, eg. png, svg
        **kwargs
            plt.errorbar arguments
            See: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.errorbar.html
        Returns
        -------
        matplotlib.pyplot.figure
        '''
        # ax.errorbar(data['x'], data['y'], xerr = data['xerr'], yerr = data['yerr'], fmt = 'r.', fillstyle = 'none', markersize = 10, capsize = 2, ecolor = 'k', label = 'Dopplershifted data')
        plt.errorbar(data['x'], data['y']/data['bunches'], xerr = data['xerr'], yerr = data['yerr']/data['bunches'], capsize = 2, ecolor = 'k', **kwargs)
        plt.legend(fontsize = 14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel('Timestamp [s]', fontsize = 14)
        plt.ylabel('Countrate per ion bunch', fontsize = 14)
        plt.tick_params(axis='both', which='major', labelsize=14)
        if save:
            try:
                fig.savefig(f'{self._SAVEPATH_FIG}{self._scan}_wm_{self._wn_channel[11:]}.{save_format}', dpi = 400, bbox_inches='tight')
            except:
                os.mkdir(self._SAVEPATH_FIG)
                fig.savefig(f'{self._SAVEPATH_FIG}{self._scan}_wm_{self._wn_channel[11:]}.{save_format}', dpi = 400, bbox_inches='tight')
        return fig

    def save_data(self, data: pd.DataFrame) -> None:
        '''Save data in .csv format

        Parameters
        ----------
        data: pd.Dataframe
            Analysed data for this scan
        '''
        try:
            data.to_csv(self._SAVEPATH + str(self._scan) + '.csv', sep = ';')
        except:
            os.mkdir(self._SAVEPATH)
            data.to_csv(self._SAVEPATH + str(self._scan) + '.csv', sep = ';')

    def plot(self, data: pd.DataFrame, fig: plt.figure, ax: Union[plt.Axes, ArrayLike], save: bool = False, save_format: str = 'png', **kwargs) -> plt.figure:
        '''Save data in .csv format

        Parameters
        ----------
        data: pd.Dataframe
            Analysed data for this scanm
        fig: plt.figure
            figure object to plot in
        ax: plt.Axes
            axes object for each subplot
        save: bool, default: False
            True if the figure needs to be saved
        save_format: str, default: 'png'
            format for the saved figure, eg. png, svg
        **kwargs
            plt.errorbar arguments
            See: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.errorbar.html
        Returns
        -------
        matplotlib.pyplot.figure
        '''
        # ax.errorbar(data['x'], data['y'], xerr = data['xerr'], yerr = data['yerr'], fmt = 'r.', fillstyle = 'none', markersize = 10, capsize = 2, ecolor = 'k', label = 'Dopplershifted data')
        plt.errorbar(data['x'], data['y']/data['bunches'], xerr = data['xerr'], yerr = data['yerr']/data['bunches'], capsize = 2, ecolor = 'k', color = '#0066FF', **kwargs)
        plt.legend(fontsize = 14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel('Observed frequency - transition frequency (MHz)', fontsize = 14)
        #plt.xlabel('Relative frequency (MHz)', fontsize = 14)
        plt.ylabel('Countrate per ion bunch', fontsize = 14)
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.title(str(self._mass)+"Sb: Scan "+str(self._scan))
        if save:
            try:
                fig.savefig(f'{self._SAVEPATH_FIG}{self._scan}_wm_{self._wn_channel[11:]}.{save_format}', dpi = 400, bbox_inches='tight')
            except:
                os.mkdir(self._SAVEPATH_FIG)
                fig.savefig(f'{self._SAVEPATH_FIG}{self._scan}_wm_{self._wn_channel[11:]}.{save_format}', dpi = 400, bbox_inches='tight')
        return fig

    def poisson_interval_high(self, data: ArrayLike, alpha: Union[int,float] = 0.32) -> ArrayLike:
        '''Top error (1-sigma) for low count statistics

        Parameters
        ----------
        data: ArrayLike
            total amount of counts
        alpha: int, float, default: 0.32
            alpha value for chi2 distribution
        Returns
        -------
        ArrayLike
        '''
        high = chi2.ppf(1 - alpha / 2, 2 * data + 2) / 2
        return high

    def poisson_interval_low(self, data: ArrayLike, alpha: Union[int,float] = 0.32) -> ArrayLike:
        '''bottom error (1-sigma) for low count statistics

        Parameters
        ----------
        data: ArrayLike
            total amount of counts
        alpha: int, float, default: 0.32
            alpha value for chi2 distribution
        Returns
        -------
        ArrayLike
        '''
        low = chi2.ppf(alpha / 2, 2 * data) / 2
        low = np.nan_to_num(low)
        return low

    def plot_asymmetric(self, data: pd.DataFrame, fig: plt.figure, ax: Union[plt.Axes, ArrayLike], save: bool = False, save_format: str = 'png', **kwargs) -> plt.figure:
        '''plot and save figure in given format with asymmetric errorbars

        Parameters
        ----------
        data: pd.Dataframe
            Analysed data for this scanm
        fig: plt.figure
            figure object to plot in
        ax: plt.Axes
            axes object for each subplot
        save: bool, default: False
            True if the figure needs to be saved
        save_format: str, default: 'png'
            format for the saved figure, eg. png, svg
        **kwargs
            plt.errorbar arguments
            See: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.errorbar.html
        Returns
        -------
        matplotlib.pyplot.figure
        '''
        yerr_h, yerr_l = (self.poisson_interval_high(data['y'])-data['y'])/data['bunches'],(data['y'] - self.poisson_interval_low(data['y']))/data['bunches']
        plt.errorbar(data['x'], data['y']/data['bunches'], xerr = data['xerr'], yerr = [yerr_l,yerr_h], capsize = 2, ecolor = 'k', color = '#0066FF', **kwargs)
        plt.legend(fontsize = 14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel('Observed frequency - transition frequency (MHz)', fontsize = 14)
        #plt.xlabel('Relative frequency (MHz)', fontsize = 14)
        plt.ylabel('Countrate per ion bunch', fontsize = 14)
        plt.tick_params(axis='both', which='major', labelsize=14)
        if save:
            try:
                fig.savefig(f'{self._SAVEPATH_FIG}{self._scan}_wm_{self._wn_channel[11:]}.{save_format}', dpi = 400, bbox_inches='tight')
            except:
                os.mkdir(self._SAVEPATH_FIG)
                fig.savefig(f'{self._SAVEPATH_FIG}{self._scan}_wm_{self._wn_channel[11:]}.{save_format}', dpi = 400, bbox_inches='tight')
        return fig
