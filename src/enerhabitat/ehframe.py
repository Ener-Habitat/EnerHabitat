import pandas as pd
import pvlib
import warnings

from datetime import datetime, timedelta, timezone

from .ehtools import *
from .config import config

# Cap on the day-to-day convergence iterations of the 1D solver (the daily
# tolerance is 5e-4 °C, same as config2d.tol_day). Mirrors config2d.max_days.
MAX_DAYS = 60


class Location:
    """
    Location class to handle climate data from an EPW file.

    Attributes:
        file (str): Path to the EPW file containing climate data.
        city (str): City of the location.
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.
        altitude (float): Altitude of the location in meters.
        timezone (datetime.timezone): Fixed-offset timezone built from the EPW
            header's decimal UTC offset (local standard time, no DST rules;
            fractional offsets such as +5.5 or +5.75 are preserved).
             

    Methods:
        meanDay(day, month, year): Calculates the ambient temperature per second for the average day
        info(): Prints Location's attributes information.
        flag(): Returns a copy of the internal flag dictionary.
    """
    
    def __init__(self, epw_file:str):
        self.file = epw_file
        
        self.__meanday_dataframe = None
        self.__flag = {"recalculate": True,
                       "date": None,
                       "day": "15",
                       "month": "current_month",
                       "year": "current_year"}
        
    def info(self):
        """
        Prints Location information.
        """
        print("<class 'enerhabitat.Location'>")
        print(f'City: {self.city}')
        print(f'Timezone: {self.timezone}')
        print(f'Latitude: {self.latitude}°')
        print(f'Longitude: {self.longitude}°')
        print(f'Altitude: {self.altitude} m')
        print(f'File: {self.file}')
    
    def meanDay(self,
        day = "15",
        month = "current",
        year = "current",
        ) -> pd.DataFrame:
        """
        Calculates the ambient temperature per second for the average day based on Location data.
        """
        if month == "current": month = datetime.now().month
        if year == "current": year = datetime.now().year

        day = str(day)
        month = str(month)
        year = str(year)
        
        if day != self.__flag['day']:
            self.__flag['day'] = day
            self.__invalidate_cache()
        if month != self.__flag['month']:
            self.__flag['month'] = month
            self.__invalidate_cache()
        if year != self.__flag['year']:
            self.__flag['year'] = year
            self.__invalidate_cache()

        recalculate = (self.__meanday_dataframe is None or 
                       self.__updated)
        
        self.__flag['recalculate'] = recalculate
        self.__flag['date'] = self.__flag['day'] + '-' + self.__flag['month'] + '-' + self.__flag['year']
        
        if recalculate:
            self.__meanday_dataframe = self.__calc_meanday()
            self.__updated = False

        return self.__meanday_dataframe

    def flag(self):
        return self.__flag.copy()
    
    def __calc_meanday(self) -> pd.DataFrame:
        """
        Calculates the ambient temperature per second for the average day based on Location data.

        Args:
            day (str, optional): Day of interest. Defaults to 15.
            month (str, optional): Month of interest. Defaults to current month.
            year (str, optional): Year of interest. Defaults to current year.

        Returns:
            DataFrame: Predicted ambient temperature ( Ta ), global ( Ig ), beam ( Ib ) 
            and diffuse irradiance ( Id ) per second for the average day of the specified month and year.
        """
        
        # print("Calculating mean day...")
        day = self.__flag['day']
        month = self.__flag['month']
        year = self.__flag['year']

        f1 = f'{year}-{month}-{day} 00:00'
        f2 = f'{year}-{month}-{day} 23:59'


        epw_data = self.__epw_format_data(year=year)

        dia_promedio = pd.date_range(start=f1, end=f2, freq='1s',tz=self.timezone)
        # tz="UTC" is inert here: get_solarposition uses the tz-aware index
        # above (pvlib's Location.tz cannot hold fractional UTC offsets).
        location = pvlib.location.Location(latitude = self.latitude,
                                           longitude= self.longitude,
                                           altitude= self.altitude,
                                           tz="UTC")

        dia_promedio = location.get_solarposition(dia_promedio)
        del dia_promedio['apparent_zenith']
        del dia_promedio['apparent_elevation']

        sunrise,_ = get_sunrise_sunset_times(dia_promedio)
        tTmax,Tmin,Tmax = calculate_tTmaxTminTmax(month, epw_data)

        # Calculate ambient temperature y add to the DataFrame
        dia_promedio = add_temperature_model(dia_promedio, Tmin, Tmax, sunrise, tTmax)

        # Add Ig, Ib, Id y Tn a dia_promedio 
        dia_promedio = add_IgIbId_Tn(dia_promedio, epw_data, month, f1, f2, self.timezone)

        # Add DeltaTn
        DeltaTa= dia_promedio.Ta.max() - dia_promedio.Ta.min()
        dia_promedio['DeltaTn'] = calculate_DtaTn(DeltaTa)

        return dia_promedio

    def __epw_format_data(self, year = None, warns = False, alias = True) -> pd.DataFrame:
        """
        Reads Location's EPW file and returns a formatted DataFrame.
            year : None default to leave intact the year or change if desired. It raises a warning.
            alias : True default, to change columns to To, Ig, Ib, Ws, RH, ...
            warns : False default, True to enable warnings.
        """
        
        names = ['Year',
                 'Month',
                 'Day',
                 'Hour',
                 'Minute',
                 'Data Source and Uncertainty Flags',
                 'Dry Bulb Temperature',
                 'Dew Point Temperature',
                 'Relative Humidity',
                 'Atmospheric Station Pressure',
                 'Extraterrestrial Horizontal Radiation',
                 'Extraterrestrial Direct Normal Radiation',
                 'Horizontal Infrared Radiation Intensity',
                 'Global Horizontal Radiation',
                 'Direct Normal Radiation',
                 'Diffuse Horizontal Radiation',
                 'Global Horizontal Illuminance',
                 'Direct Normal Illuminance',
                 'Diffuse Horizontal Illuminance',
                 'Zenith Luminance',
                 'Wind Direction',
                 'Wind Speed',
                 'Total Sky Cover',
                 'Opaque Sky Cover',
                 'Visibility',
                 'Ceiling Height',
                 'Present Weather Observation',
                 'Present Weather Codes',
                 'Precipitable Water',
                 'Aerosol Optical Depth',
                 'Snow Depth',
                 'Days Since Last Snowfall',
                 'Albedo',
                 'Liquid Precipitation Depth',
                 'Liquid Precipitation Quantity']

        rename = {'Dry Bulb Temperature'       :'To',
                 'Relative Humidity'           :'RH',
                 'Atmospheric Station Pressure':'P' ,
                 'Global Horizontal Radiation' :'Ig',
                 'Direct Normal Radiation'     :'Ib',
                 'Diffuse Horizontal Radiation':'Id',
                 'Wind Direction'              :'Wd',
                 'Wind Speed'                  :'Ws'}

        data = pd.read_csv(self.file, skiprows=8, header=None, names=names, usecols=range(35))
        data.Hour = data.Hour -1
        if year != None:
            data.Year = year
            if warns == True:
                warnings.warn("Year has been changed, be carefull")
        try:
            data['tiempo'] = data.Year.astype('str') + '-' + data.Month.astype('str')  + '-' + data.Day.astype('str') + ' ' + data.Hour.astype('str') + ':' + data.Minute.astype('str') 
            data.tiempo = pd.to_datetime(data.tiempo,format='%Y-%m-%d %H:%M')
        except:
            data.Minute = 0
            data['tiempo'] = data.Year.astype('str') + '-' + data.Month.astype('str')  + '-' + data.Day.astype('str') + ' ' + data.Hour.astype('str') + ':' + data.Minute.astype('str') 
            data.tiempo = pd.to_datetime(data.tiempo,format='%Y-%m-%d %H:%M')

        data.set_index('tiempo',inplace=True)
        del data['Year']
        del data['Month']
        del data['Day']
        del data['Hour']
        del data['Minute']
        if alias:
            data.rename(columns=rename,inplace=True)
        
        return data

    def __invalidate_cache(self):
        self.__updated = True
        
    @property
    def file(self):
        return self.__epw_path
    @file.setter
    def file(self, file):
        """
        EPW file containing climate data. Attributes timezone, longitude, latitude, altitude are taken from this file.
        """
        datos=[]
        
        with open(file,'r') as epw:
            datos=epw.readline().split(',')
            
        self.__epw_path = file
        self.__city = str(datos[1]) + ", " + str(datos[2])
        self.__latitude = float(datos[6])
        self.__longitude = float(datos[7])
        self.__altitude = float(datos[9])
        
        # EPW field 8 is the decimal UTC offset (e.g. -6.0, +5.5, +5.75).
        # A stdlib fixed offset preserves fractional zones exactly;
        # 'Etc/GMT±N' only exists for integer offsets (and used to truncate
        # +5.5 → +5).
        tmz = float(datos[8])
        self.__timezone = timezone(timedelta(minutes=round(tmz * 60)))
        
        self.__invalidate_cache()
    
    @property
    def city(self):
        return self.__city
    @city.setter
    def city(self, value):
        pass
    
    @property
    def timezone(self):
        return self.__timezone
    @timezone.setter
    def timezone(self, value):
        pass
    
    @property
    def latitude(self):
        return self.__latitude
    @latitude.setter
    def latitude(self, value):
        pass
    
    @property
    def longitude(self):
        return self.__longitude
    @longitude.setter
    def longitude(self, value):
        pass
    
    @property
    def altitude(self):
        return self.__altitude
    @altitude.setter
    def altitude(self, value):
        pass    
    
class System():
    """
    System class to model a constructive system and calculate its interior temperature
    based on the sun-air temperature experienced by the surface.
    Attributes:
        location (Location): Location object containing climate data.
        tilt (float): Tilt angle of the surface in degrees.
        azimuth (float): Azimuth angle of the surface in degrees.
        absortance (float): Surface absortance of the system's external material.
        layers (list): List of tuples from outside to inside with material and width.
        energy_transfer (float): Total energy transfered to the constuctive system computed by solve()
        cooling_energy (float): Total cooling energy needed to mantain a constant Ti computed by solveAC()
        heating_energy (float): Total heating energy needed to mantain a constant Ti computed by solveAC()
        
    Methods:
        Tsa(): Calculates the sun-air temperature at the dt grid for the average day experienced by a surface.
        solve(energy): Solves the constructive system's inside temperature.
        solveAC(): Solves the constructive system's required cooling and heating energy to maintain the inside temperature.
        info(): Prints System information.
        add_layer(material, width): Adds a layer to the constructive system.
        remove_layer(index): Removes a layer from the constructive system by index.
    """
    
    def __init__(self, location:Location , tilt = 90, azimuth = 0, absortance = 0.8, layers = []):
        self.tilt = tilt
        self.azimuth = azimuth
        self.absortance = absortance
        self.location = location
        self.layers= layers
        
        # Caches governed by signatures (same pattern as System2D): the current
        # inputs are hashed into a tuple and compared on every call, so ANY
        # change (including mutating `layers` in place) triggers recomputation,
        # and solve()/solveAC() keep separate caches per mode.
        self.__tsa_dataframe = None
        self.__tsa_sig = None
        self.__free_df = None
        self.__free_sig = None
        self.__free_state = None
        self.__ac_df = None
        self.__ac_sig = None
        self.__ac_state = None
        self.__Tso = None
        self.__Tsi = None
        # AC setpoint: if None, solveAC() holds Ti at Tn.mean() (like the 2D)
        self.setpoint = None

        self.__last_solve = None
        self.__energy_transfer = None
        self.__cooling_energy = None
        self.__heating_energy = None
        self.__days = None
        self.__day_error = None
        self.__converged = None
        self.__energy_imbalance = None
        
    def Tsa(self,
            # solar_absortance:float=None,
            # surface_tilt:float=None,
            # surface_azimuth:float=None
            ) -> pd.DataFrame: 
        """
        Sun-air temperature for the average day experienced by a surface based
        on a meanDay dataframe from System's Location (Ta, Ig, Ib and Id).

        Sampled at the solver time step (config.dt), the same grid as solve(),
        so Ti, Tsa, Is and Ig share identical timestamps and a plain concat
        aligns without NaN.

        Returns:
            DataFrame: Predicted sun-air temperature ( Tsa ) and solar irradiance ( Is )
            for the average day, at the dt grid.
        """
        
        """
        if solar_absortance is not None:
            self.absortance = solar_absortance
        if surface_tilt is not None:
            self.tilt = surface_tilt
        if surface_azimuth is not None:
            self.azimuth = surface_azimuth
        """
        
        if self.location.flag()['date'] is None:
            self.location.meanDay()      # never selected: default current month
        sig = self._tsa_signature()
        if self.__tsa_dataframe is None or self.__tsa_sig != sig:
            self.__tsa_dataframe = self.__calc_tsa()
            self.__tsa_sig = sig
        return self.__tsa_dataframe

    def _tsa_signature(self):
        """Inputs that determine ``Tsa()``: location (identity and selected
        mean day), orientation, absorptance and the global config."""
        return (id(self.location), self.location.flag().get("date"),
                self.tilt, self.azimuth, self.absortance, config.version)

    def _solve_signature(self):
        """Inputs that determine a solve: the Tsa inputs plus the layers
        (read at call time, so in-place mutations are detected)."""
        return self._tsa_signature() + (tuple(tuple(l) for l in self.layers),)

    def __results_snapshot(self):
        return (self.__Tso, self.__Tsi, self.__energy_transfer,
                self.__cooling_energy, self.__heating_energy, self.__days,
                self.__day_error, self.__converged, self.__energy_imbalance,
                self.__last_solve)

    def __results_restore(self, state):
        (self.__Tso, self.__Tsi, self.__energy_transfer,
         self.__cooling_energy, self.__heating_energy, self.__days,
         self.__day_error, self.__converged, self.__energy_imbalance,
         self.__last_solve) = state
    
    def solve(self) -> pd.Series:
        """
        Solves the constructive system's inside temperature with the Tsa simulation dataframe.

        Returns:
            Ti (pandas.Series named "Ti"): Interior temperature for the
            constructive system, on the ``config.dt`` time grid. Use
            ``.to_frame("Ti")`` when a DataFrame is required. Also stores
            ``energy_transfer``, ``days``, ``day_error``, ``converged`` and
            ``energy_imbalance`` on the instance.
        """
        if len(self.layers) == 0:
            raise ValueError("Constructive system layers are not defined.")

        sig = self._solve_signature()
        if self.__free_df is not None and self.__free_sig == sig:
            self.__results_restore(self.__free_state)
            return self.__free_df

        self.__free_df = self.__calc_solve(AC=False)
        self.__free_sig = sig
        self.__free_state = self.__results_snapshot()
        return self.__free_df
    
    def solveAC(self) -> pd.Series:
        """
        Solves the constructive system's required cooling and heating energy to
        maintain the interior temperature with the Tsa simulation dataframe.
        The indoor temperature is held at ``setpoint`` if set, otherwise at the
        neutrality temperature ``Tn.mean()``.

        Returns:
            Ti (pandas.Series named "Ti"): Interior temperature for the
            constructive system. The energies are stored on the instance as
            ``cooling_energy`` and ``heating_energy`` (J/(m²·day)), along with
            ``days``, ``day_error`` and ``converged``.
        """
        if len(self.layers) == 0:
            raise ValueError("Constructive system layers are not defined.")

        sig = self._solve_signature() + (self.setpoint,)
        if self.__ac_df is not None and self.__ac_sig == sig:
            self.__results_restore(self.__ac_state)
            return self.__ac_df

        self.__ac_df = self.__calc_solve(AC=True)
        self.__ac_sig = sig
        self.__ac_state = self.__results_snapshot()
        return self.__ac_df

    def info(self):
        """
        Prints System information.
        """
        print("<class 'enerhabitat.System'>")
        print(f"Location: {self.location.city}")
        print(f"meanDay date: {self.location.flag()['date']}")
        print(f"Tilt: {self.tilt}°")
        print(f"Azimuth: {self.azimuth}°")
        print(f"Absortance: {self.absortance}")
        print(f"Energy transfer: {self.energy_transfer}")
        print(f"Cooling energy: {self.cooling_energy}")
        print(f"Heating energy: {self.heating_energy}")
        if len(self.layers) != 0:
            print("Layers:")
            for i, (material, width) in enumerate(self.layers):
                print(f"\t{i+1}: {material}, {width} m")
        else:
            print("Layers: No layers defined")
    
    def add_layer(self, material:str, width:float):
        """
        Adds a layer to the constructive system.

        Args:
            material (str): Material name.
            width (float): Width of the material in meters.
        """
        self.__layers.append((material, width))
        return self.layers
    
    def remove_layer(self, index:int):
        """
        Removes a layer from the constructive system by index.

        Args:
            index (int): Positive index of the layer to remove.
        """
        if index < 0 or index >= len(self.__layers):
            raise IndexError("Layer index out of range.")
        del self.__layers[index]
        return self.layers
    
    def __calc_tsa(self) -> pd.DataFrame:

        # Fetch the location's CURRENT mean-day selection (re-passing its own
        # flag values does not change the selection, so the cached frame is
        # reused); the System follows the Location, as the 2D does. Tsa()
        # guarantees a selection exists before calling here.
        f = self.location.flag()
        tsa_dataframe = self.location.meanDay(day=f['day'], month=f['month'],
                                              year=f['year']).copy()
            
        absortance = self.absortance
        tilt = self.tilt
        azimuth = self.azimuth
        
        outside_convection_heat_transfer = config.ho

        # Long-wave radiation factor RF: 3.9 °C for a horizontal roof facing
        # the cold sky, decreasing linearly with tilt to 0 at a vertical wall
        # (the 2016 online tool's rule), 0 beyond 90°. Identical to the
        # previous binary rule at tilt = 0 and tilt = 90.
        LWR = 3.9 * max(0.0, 1.0 - tilt / 90.0)

        total_irradiance = pvlib.irradiance.get_total_irradiance(
            surface_tilt=tilt,
            surface_azimuth=azimuth,
            dni=tsa_dataframe['Ib'],
            ghi=tsa_dataframe['Ig'],
            dhi=tsa_dataframe['Id'],
            solar_zenith=tsa_dataframe['zenith'],
            solar_azimuth=tsa_dataframe['azimuth']
        )

        # Add Is
        tsa_dataframe['Is'] = total_irradiance.poa_global

        # Add Tsa
        tsa_dataframe['Tsa'] = tsa_dataframe.Ta + tsa_dataframe.Is*absortance/outside_convection_heat_transfer - LWR

        # Tsa and solve() share the same time grid: it is computed on the mean
        # day per second (required by pvlib) and returned already subsampled to
        # the solver time step dt, so that Ti, Tsa, Is and Ig live at the same
        # instants and a concat aligns without NaN.
        tsa_dataframe = tsa_dataframe.iloc[::config.dt]

        return tsa_dataframe

    def __calc_solve(self, AC=False) -> pd.DataFrame:
        """
        Solves the constructive system's inside temperature with the Tsa simulation dataframe.

        Args:
            constructive_system (list): List of tuples from outside to inside with material and width.
            Tsa_dataframe (DataFrame): Predicted sun-air temperature ( Tsa ) at the dt grid for the average day DataFrame.

        Returns:
            Ti (pandas.Series named "Ti"): Interior temperature for the
            constructive system. The energies are stored on the instance:
            ``energy_transfer`` (free-running) or ``cooling_energy`` /
            ``heating_energy`` (AC), plus the convergence diagnostics.
        """
        
        La = config.La # Length of the dummy frame
        Nx = config.Nx # Number of elements to discretize
        ho = config.ho # Outside convection heat transfer
        hi = config.hi # Inside convection heat transfer
        dt = config.dt # Time step
        AIR_DENSITY = config.AIR_DENSITY
        AIR_HEAT_CAPACITY = config.AIR_HEAT_CAPACITY

        SC_dataframe = self.Tsa().copy()
        constructive_system = self.layers
        
        materials = config.materials

        # Roof-like surfaces: hi depends on the heat-flow direction at the
        # indoor surface, evaluated every time step (NOM-020): hi_down when
        # Tsi > Ti (stable stratification, downward flow), hi_up when Tsi <= Ti
        # (unstable, upward flow). The wall treatment (fixed hi) applies for
        # tilt >= 60 deg — the same boundary EnergyPlus uses to classify a
        # surface as roof (tilt < 60) vs wall (60-120), and equivalent to the
        # ISO 6946 rule (heat flow within +-30 deg of the horizontal is
        # "horizontal"). config.hi_flow = False forces the fixed hi everywhere.
        hi_flow = bool(config.hi_flow) and float(self.tilt) < 60.0
        hi_up = config.hi_up if hi_flow else hi
        hi_down = config.hi_down if hi_flow else hi

        cs = set_construction(materials, constructive_system)
        k, rhoc, dx, Gf = set_k_rhoc(cs, Nx)
        # hi only enters the last-node diagonal: build the "up" set and derive
        # the "down" one from it.
        mass_coeff, a_up, b_static, c_static = prepare_static_coefficients(k, rhoc, dx, dt, ho, hi_up, Gf)
        a_down = a_up.copy()
        a_down[Nx - 1] += hi_down - hi_up

        d = np.empty(Nx)
        P = np.empty(Nx)
        Q = np.empty(Nx)
        Tn_aux = np.empty(Nx)
        cap_up = hi_up * dt / (AIR_DENSITY * AIR_HEAT_CAPACITY * La)
        cap_down = hi_down * dt / (AIR_DENSITY * AIR_HEAT_CAPACITY * La)

        T = np.full(Nx, SC_dataframe.Tn.mean())
        # AC: the setpoint (System.setpoint) overrides the neutrality
        # temperature, as in the 2D; free-running ignores it and Tn.mean()
        # is only the initial condition of the indoor-air node.
        if AC and self.setpoint is not None:
            SC_dataframe['Ti'] = float(self.setpoint)
        else:
            SC_dataframe['Ti'] = SC_dataframe.Tn.mean()

        # Tsa() already comes at step dt (see __calc_tsa); no resampling here.
        Tsa_vals = SC_dataframe['Tsa'].to_numpy()
        Ti_vals = SC_dataframe['Ti'].to_numpy(copy=True)
        n_steps = Tsa_vals.shape[0]
        Tso_vals = np.empty(n_steps)
        Tsi_vals = np.empty(n_steps)

        tol_day = 5e-4          # °C, same as config2d.tol_day
        days = 0
        day_error = 1.0e9

        if AC:  # AC = True
            # Periodic convergence over ALL persisted states: the solid field
            # AND the indoor series Ti_vals (it feeds the next day's steps).
            Ti_prev = Ti_vals.copy()
            while day_error > tol_day and days < MAX_DAYS:
                Told = T.copy()
                Qcool = Qheat = 0.
                for idx in range(n_steps):
                    if T[Nx - 1] > Ti_vals[idx]:
                        hi_s, a_s = hi_down, a_down
                    else:
                        hi_s, a_s = hi_up, a_up
                    calculate_coefficients(mass_coeff, T, Tsa_vals[idx], ho, Ti_vals[idx], hi_s, d)
                    # Llamado de funcion para Acc
                    T, Ti = solve_PQ_AC(a_s, b_static, c_static, d, T, Nx, Ti_vals[idx], P, Q, Tn_aux)
                    Tso_vals[idx] = T[0]
                    Tsi_vals[idx] = T[Nx - 1]
                    if (T[Nx-1] > Ti):
                        Qcool += hi_s*dt*(T[Nx-1]-Ti)
                    if (T[Nx-1] < Ti):
                        Qheat += hi_s*dt*(Ti-T[Nx-1])
                    Ti_vals[idx] = Ti
                C = np.abs(Told - T).mean()
                dTi = np.abs(Ti_vals - Ti_prev).max()
                Ti_prev = Ti_vals.copy()
                day_error = max(C, dTi)
                days += 1

            SC_dataframe['Ti'] = Ti_vals
            self.__Tso = pd.Series(Tso_vals, index=SC_dataframe.index, name="Tso")
            self.__Tsi = pd.Series(Tsi_vals, index=SC_dataframe.index, name="Tsi")

            self.__last_solve = 'ac'
            self.__energy_transfer = None
            self.__cooling_energy = Qcool
            self.__heating_energy = Qheat
            self.__store_convergence(days, day_error, tol_day, None, "solveAC")

            return SC_dataframe['Ti']

        else:
            # Free-running: the indoor air is a lumped-capacitance node
            # integrated IN TIME. 'tint' is a scalar that advances step by step
            # and persists between iterations (like the wall T, and like the
            # original C solver). In the oscillatory steady state the cycle
            # closes (Qin == Qout) and the transferred energy is either one; it
            # is measured at the indoor surface [Nx-1] with hi and the air at
            # the same instant (tinn, before the update), integrating over dt.
            tint = float(Ti_vals[0])   # = Tn.mean(), initialized only once
            Qin = Qout = 0.0
            # Periodic convergence over ALL persisted states: the solid field
            # AND the indoor-air node tint (it carries over between days).
            while day_error > tol_day and days < MAX_DAYS:
                Told = T.copy()
                tint_prev = tint
                Qin = Qout = 0.0
                for idx in range(n_steps):
                    tinn = tint
                    if T[Nx - 1] > tinn:
                        hi_s, a_s, cap_s = hi_down, a_down, cap_down
                    else:
                        hi_s, a_s, cap_s = hi_up, a_up, cap_up
                    calculate_coefficients(mass_coeff, T, Tsa_vals[idx], ho, tinn, hi_s, d)
                    T, tint = solve_PQ(a_s, b_static, c_static, d, T, Nx, tinn, cap_s, P, Q, Tn_aux)
                    Ti_vals[idx] = tint
                    Tso_vals[idx] = T[0]
                    Tsi_vals[idx] = T[Nx - 1]
                    flux = hi_s * (T[Nx - 1] - tinn) * dt
                    if flux > 0:
                        Qin += flux
                    else:
                        Qout -= flux
                C = np.abs(Told - T).mean()
                dTi = abs(tint - tint_prev)
                day_error = max(C, dTi)
                days += 1

            SC_dataframe['Ti'] = Ti_vals
            self.__Tso = pd.Series(Tso_vals, index=SC_dataframe.index, name="Tso")
            self.__Tsi = pd.Series(Tsi_vals, index=SC_dataframe.index, name="Tsi")

            self.__last_solve = 'temp'
            self.__energy_transfer = Qin   # oscillatory steady state: Qin == Qout
            self.__cooling_energy = None
            self.__heating_energy = None
            # Energy-closure diagnostic: in the periodic regime Qin == Qout.
            qmax = max(Qin, Qout)
            imbalance = abs(Qin - Qout) / qmax if qmax > 0.0 else 0.0
            self.__store_convergence(days, day_error, tol_day, imbalance, "solve")

            return SC_dataframe['Ti']

    def __store_convergence(self, days, day_error, tol_day, imbalance, who):
        """Store the daily-convergence diagnostics and warn if not converged."""
        self.__days = days
        self.__day_error = float(day_error)
        self.__converged = day_error <= tol_day
        self.__energy_imbalance = imbalance
        if not self.__converged:
            warnings.warn(
                f"System.{who}: not converged after {days} days "
                f"(day_error={day_error:.3e} °C > tol {tol_day}); results may "
                f"not be periodic. See enerhabitat.ehframe.MAX_DAYS.",
                RuntimeWarning)

    # Los setters no invalidan nada explícitamente: las firmas de caché leen
    # los valores actuales en cada llamada (ver _tsa_signature/_solve_signature).
    @property
    def layers(self):
        return self.__layers
    @layers.setter
    def layers(self, layers:list):
        """
        List of tuples from outside to inside with material and width.
        Example: [('Brick',0.1), ('Insulation',0.05), ('Adobe',0.02)]
        """
        self.__layers = layers

    @property
    def location(self):
        return self.__instance_location
    @location.setter
    def location(self, loc:Location):
        """
        Location object containing climate data.
        """
        self.__instance_location = loc

    @property
    def tilt(self):
        return self.__tilt
    @tilt.setter
    def tilt(self, angle:float):
        """
        Tilt angle of the surface in degrees.
        """
        self.__tilt = angle

    @property
    def azimuth(self):
        return self.__azimuth
    @azimuth.setter
    def azimuth(self, angle:float):
        """
        Azimuth angle of the surface in degrees.
        """
        self.__azimuth = angle

    @property
    def absortance(self):
        return self.__absortance
    @absortance.setter
    def absortance(self, value:float):
        """
        Surface absortance of the system's external material.

        Raises:
            ValueError: if ``value`` is outside ``[0, 1]``.
        """
        if not (0.0 <= value <= 1.0):
            raise ValueError(f"absortance must be in [0, 1], got {value!r}")
        self.__absortance = value

    # Solo lectura
    @property
    def energy_transfer(self):
        return self.__energy_transfer
    @energy_transfer.setter
    def energy_transfer(self, value):
        pass
    
    @property
    def heating_energy(self):
        return self.__heating_energy
    @heating_energy.setter
    def heating_energy(self, value):
        pass
    
    @property
    def cooling_energy(self):
        return self.__cooling_energy
    @cooling_energy.setter
    def cooling_energy(self, value):
        pass

    @property
    def days(self):
        """Day-to-day iterations used by the last solve."""
        return self.__days
    @days.setter
    def days(self, value):
        pass

    @property
    def day_error(self):
        """Final day-to-day error (°C) of the last solve: max of the solid-field
        mean change and the indoor-state change |ΔT_i|."""
        return self.__day_error
    @day_error.setter
    def day_error(self, value):
        pass

    @property
    def converged(self):
        """True if the last solve reached the periodic regime (day_error ≤ 5e-4 °C)
        within MAX_DAYS; a RuntimeWarning is emitted otherwise."""
        return self.__converged
    @converged.setter
    def converged(self, value):
        pass

    @property
    def energy_imbalance(self):
        """Relative energy-closure diagnostic of the last free-running solve:
        |Qin − Qout| / max(Qin, Qout) (0 = perfect periodic closure).
        None after solveAC()."""
        return self.__energy_imbalance
    @energy_imbalance.setter
    def energy_imbalance(self, value):
        pass

    @property
    def Tso(self):
        """Outdoor-surface temperature Series (named "Tso") of the last solve,
        on the Tsa() time grid; None before solving. Concatenates directly:
        pd.concat([wall.solve(), wall.Tso, wall.Tsi, wall.Tsa()], axis=1)."""
        return self.__Tso
    @Tso.setter
    def Tso(self, value):
        pass

    @property
    def Tsi(self):
        """Indoor-surface temperature Series (named "Tsi") of the last solve,
        on the Tsa() time grid; None before solving."""
        return self.__Tsi
    @Tsi.setter
    def Tsi(self, value):
        pass
