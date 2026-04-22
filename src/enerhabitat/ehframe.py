import pandas as pd
import pvlib
import pytz
import warnings

from datetime import datetime

from .ehtools import *
from .config import config

class Location:
    """
    Location class to handle climate data from an EPW file.

    Attributes:
        file (str): Path to the EPW file containing climate data.
        city (str): City of the location.
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.
        altitude (float): Altitude of the location in meters.
        timezone (pytz.timezone): Timezone of the location.
             

    Methods:
        meanDay(day, month, year): Calculates the ambient temperature per second for the average day
        info(): Prints Location's attributes information.
        copy(): Returns a copy of the Location instance.
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

    def copy(self):
        """
        Returns a copy of the Location instance.
        """
        return Location(self.file)
    
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
        location = pvlib.location.Location(latitude = self.latitude, 
                                           longitude= self.longitude, 
                                           altitude= self.altitude,
                                           tz=self.timezone)

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
        
        tmz = int(datos[8].split('.')[0])
        self.__timezone = pytz.timezone('Etc/GMT'+f'{(-tmz):+}')
        
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
        Tsa(): Calculates the sun-air temperature per second for the average day experienced by a surface.
        solve(energy): Solves the constructive system's inside temperature.
        solveAC(): Solves the constructive system's required cooling and heating energy to maintain the inside temperature.
        info(): Prints System information.
        copy(): Returns a copy of the System instance.
        flag(): Returns a copy of the internal flag dictionary.
        add_layer(material, width): Adds a layer to the constructive system.
        remove_layer(index): Removes a layer from the constructive system by index.
    """
    
    def __init__(self, location:Location , tilt = 90, azimuth = 0, absortance = 0.8, layers = []):
        self.tilt = tilt
        self.azimuth = azimuth
        self.absortance = absortance
        self.location = location
        self.layers= layers
        
        self.__flag = {"recalculate": True,
                       "tsa_date": None,
                       "solve_date": None,
                       "config": config.to_dict()
                       }
        # self.__update_flag_config()
        self.__tsa_dataframe = None
        self.__solve_dataframe = None
        
        self.__tsa_solver_version = -1
        self.__solve_solver_version = -1
        self.__last_solve = None
        
        self.__energy_transfer = None
        self.__cooling_energy = None
        self.__heating_energy = None
        self.__invalidate_cache()
        
    def Tsa(self,
            # solar_absortance:float=None,
            # surface_tilt:float=None,
            # surface_azimuth:float=None
            ) -> pd.DataFrame: 
        """
        Sun-air temperature per second for the average day experienced
        by a surface based on a meanDay dataframe from System's Location
        (Ta, Ig, Ib and Id).

        Returns:
            DataFrame: Predicted sun-air temperature ( Tsa ) and solar irradiance ( Is )
            per second for the average day.
        """
        
        """
        if solar_absortance is not None:
            self.absortance = solar_absortance
        if surface_tilt is not None:
            self.tilt = surface_tilt
        if surface_azimuth is not None:
            self.azimuth = surface_azimuth
        """
        
        mean_date = self.location.flag()["date"]  # Asegura que el DataFrame del día medio esté actualizado
        
        if  self.__flag["tsa_date"] != mean_date:
            self.__flag["tsa_date"] = mean_date
            self.__invalidate_cache()
        
        self.__update_flag_config()
        
        recalculate = (self.__tsa_dataframe is None or
                       self.__updated or 
                       self.__tsa_solver_version != config.version
                       )
        
        if recalculate:
            self.__tsa_dataframe = self.__calc_tsa()  # el método que calcula Tsa
            self.__updated = False
        
        self.__flag['recalculate'] = recalculate
        return self.__tsa_dataframe
    
    def solve(self) -> pd.DataFrame:
        """
        Solves the constructive system's inside temperature with the Tsa simulation dataframe.

        Args:
            energy (bool): If True, returns also the energy transfer ET.
        
        Returns:
            Ti (DataFrame): Interior temperature for the constructive system.
        """
        constructive_system = self.layers
        if len(constructive_system) == 0:
            raise ValueError("Constructive system layers are not defined.")
        
        if self.__flag["tsa_date"] != self.__flag['solve_date']:
            self.__flag["solve_date"] = self.__flag['tsa_date']
            self.__invalidate_cache()
        
        self.__update_flag_config()
        
        recalculate = (self.__updated or 
                        self.__solve_dataframe is None or 
                        self.__solve_solver_version != config.version or
                        self.__last_solve != 'temp'
                        )
        
        self.__flag['recalculate'] = recalculate
        
        if recalculate:    
            self.__solve_dataframe = self.__calc_solve(AC=False)
            self.__updated = False
            
        return self.__solve_dataframe
    
    def solveAC(self) -> pd.DataFrame:
        """
        Solves the constructive system's required cooling and heating energy to 
        maintain the interior temperature with the Tsa simulation dataframe.

        Returns:
            Ti (DataFrame): Interior temperature for the constructive system.
            Qcool, Qheat (float): Cooling energy and heating energy values.
        """
        constructive_system = self.layers
        if len(constructive_system) == 0:
            raise ValueError("Constructive system layers are not defined.")
        
        if self.__flag["tsa_date"] != self.__flag['solve_date']:
            self.__flag["solve_date"] = self.__flag['tsa_date']
            self.__invalidate_cache()
        
        self.__update_flag_config()
        
        recalculate = (self.__updated or 
                        self.__solve_dataframe is None or 
                        self.__solve_solver_version != config.version or
                        self.__last_solve != 'ac'
                       )

        self.__flag['recalculate'] = recalculate
        
        if recalculate:
            self.__solve_dataframe = self.__calc_solve(AC=True)
            self.__updated = False
        
        return self.__solve_dataframe

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
    
    def copy(self):
        """
        Returns a copy of the System instance.
        """
        return System(
            location=self.location.copy(),
            tilt=self.tilt,
            azimuth=self.azimuth,
            absortance=self.absortance,
            layers=self.layers.copy()
        )
    
    def flag(self):
       return self.__flag.copy()
    
    def add_layer(self, material:str, width:float):
        """
        Adds a layer to the constructive system.

        Args:
            material (str): Material name.
            width (float): Width of the material in meters.
        """
        self.__capas.append((material, width))
        self.__invalidate_cache()
        return self.layers
    
    def remove_layer(self, index:int):
        """
        Removes a layer from the constructive system by index.

        Args:
            index (int): Positive index of the layer to remove.
        """
        if index < 0 or index >= len(self.__capas):
            raise IndexError("Layer index out of range.")
        del self.__capas[index]
        self.__invalidate_cache()
        return self.layers
    
    def __update_flag_config(self):
        """
        Checks for changes and updates the flag["config"] in the System instance.
        """
        if self.__flag["config"] != config.to_dict():
            self.__flag["config"] = config.to_dict()
            self.__invalidate_cache()

    def __calc_tsa(self) -> pd.DataFrame:
        
        if self.__flag['tsa_date'] is not None:
            tsa_date = self.__flag['tsa_date'].split('-')
            tsa_dataframe = self.location.meanDay(day=tsa_date[0], month=tsa_date[1], year=tsa_date[2]).copy()
        else:
            tsa_dataframe = self.location.meanDay().copy()
            
        absortance = self.absortance
        tilt = self.tilt
        azimuth = self.azimuth
        
        outside_convection_heat_transfer = config.ho

        if tilt == 0:
            LWR = 3.9
        else:
            LWR = 0.

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

        self.__tsa_solver_version = config.version
        
        return tsa_dataframe

    def __calc_solve(self, AC=False) -> pd.DataFrame:
        """
        Solves the constructive system's inside temperature with the Tsa simulation dataframe.

        Args:
            constructive_system (list): List of tuples from outside to inside with material and width.
            Tsa_dataframe (DataFrame): Predicted sun-air temperature ( Tsa ) per second for the average day DataFrame.

        Returns:
            Ti (DataFrame): Interior temperature for the constructive system.
            ET (float): Energy transfer if energia=True.
            Qcool, Qheat (float): Cooling energy and heating energy values if AC=True.
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
        
        propiedades = config.materials

        cs = set_construction(propiedades, constructive_system)
        k, rhoc, dx = set_k_rhoc(cs, Nx)
        mass_coeff, a_static, b_static, c_static = prepare_static_coefficients(k, rhoc, dx, dt, ho, hi)

        d = np.empty(Nx)
        P = np.empty(Nx)
        Q = np.empty(Nx)
        Tn_aux = np.empty(Nx)
        capacitance_factor = hi * dt / (AIR_DENSITY * AIR_HEAT_CAPACITY * La)

        T = np.full(Nx, SC_dataframe.Tn.mean())
        SC_dataframe['Ti'] = SC_dataframe.Tn.mean()

        SC_dataframe = SC_dataframe.iloc[::dt]
        Tsa_vals = SC_dataframe['Tsa'].to_numpy()
        Ti_vals = SC_dataframe['Ti'].to_numpy(copy=True)
        Ti_new = np.empty_like(Ti_vals)
        n_steps = Tsa_vals.shape[0]

        C = 1
        ET = 0.0

        self.__solve_solver_version = config.version
        
        if AC:  # AC = True
            while C > 5e-4: 
                Told = T.copy()
                Qcool = Qheat = 0.
                for idx in range(n_steps):
                    calculate_coefficients(mass_coeff, T, Tsa_vals[idx], ho, Ti_vals[idx], hi, d)
                    # Llamado de funcion para Acc
                    T, Ti = solve_PQ_AC(a_static, b_static, c_static, d, T, Nx, Ti_vals[idx], hi, La, dt)
                    if (T[Nx-1] > Ti):
                        Qcool += hi*dt*(T[Nx-1]-Ti)
                    if (T[Nx-1] < Ti):
                        Qheat += hi*dt*(Ti-T[Nx-1])
                    Ti_vals[idx] = Ti
                Tnew = T.copy()
                C = abs(Told - Tnew).mean()

            SC_dataframe['Ti'] = Ti_vals
            
            self.__last_solve = 'ac'
            self.__energy_transfer = None
            self.__cooling_energy = Qcool
            self.__heating_energy = Qheat
            
            return SC_dataframe['Ti']

        else:
            while C > 5e-4: 
                Told = T.copy()
                ET_iter = 0.
                for idx in range(n_steps):
                    tint_prev = Ti_vals[idx]
                    calculate_coefficients(mass_coeff, T, Tsa_vals[idx], ho, tint_prev, hi, d)
                    T, tint_new = solve_PQ(a_static, b_static, c_static, d, T, Nx, tint_prev, capacitance_factor, P, Q, Tn_aux)
                    Ti_new[idx] = tint_new
                    if T[Nx-1] > tint_new:
                        ET_iter += hi * (T[Nx - 1] - tint_new) * dt
                Ti_vals[:] = Ti_new
                C = np.abs(Told - T).mean()
                ET = ET_iter

            SC_dataframe['Ti'] = Ti_vals
            
            self.__last_solve = 'temp'
            self.__energy_transfer = ET
            self.__cooling_energy = None
            self.__heating_energy = None
                        
            return SC_dataframe['Ti']

    def __invalidate_cache(self):
        self.__updated = True
        
    @property
    def layers(self):
        return self.__capas
    @layers.setter
    def layers(self, capas:list):
        """
        List of tuples from outside to inside with material and width.
        Example: [('Brick',0.1), ('Insulation',0.05), ('Adobe',0.02)]
        """
        self.__capas = capas
        self.__invalidate_cache()
    
    @property
    def location(self):
        return self.__instance_location
    @location.setter
    def location(self, loc:Location):
        """
        Location object containing climate data.
        """
        self.__instance_location = loc
        self.__invalidate_cache()

    @property
    def tilt(self):
        return self.__tilt
    @tilt.setter
    def tilt(self, angle:float):
        """
        Tilt angle of the surface in degrees.
        """
        if angle != getattr(self, "__tilt", None):
            self.__tilt = angle
            self.__invalidate_cache()
        
    @property
    def azimuth(self):
        return self.__azimuth
    @azimuth.setter
    def azimuth(self, angle:float):
        """
        Azimuth angle of the surface in degrees.
        """
        if angle != getattr(self, "__azimuth", None):
            self.__azimuth = angle
            self.__invalidate_cache()

    @property
    def absortance(self):
        return self.__absortance
    @absortance.setter
    def absortance(self, value:float):
        """
        Surface absortance of the system's external material.
        """
        if value != getattr(self, "__absortance", None):
            self.__absortance = value
            self.__invalidate_cache()

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
