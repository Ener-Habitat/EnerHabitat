import configparser
import os

class Material:
    def __init__(self, k, rho, c):
        self.__k = k
        self.__rho = rho
        self.__c = c
    
    def to_dict(self):
        return {"k": self.k,
                "rho": self.rho,
                "c": self.c
                }
    
    @property
    def k(self):
        return self.__k
    @k.setter
    def k(self, value):
        pass
    
    @property
    def rho(self):
        return self.__rho
    @rho.setter
    def rho(self, value):
        pass
    
    @property
    def c(self):
        return self.__c
    @c.setter
    def c(self, value):
        pass

class Config:
    """
    Global configuration class for EnerHabitat simulations.
    
    Attributes:
        file (str): Path to the materials configuration file.
        La (float): Length of the dummy frame (m).
        Nx (int): Number of discretization elements.
        ho (float): Outdoor convective coefficient (W/m²K).
        hi (float): Indoor convective coefficient (W/m²K).
        dt (float): Time step (seconds). Fixed at 10 s, not configurable.
        AIR_DENSITY (float): Density of air (kg/m³).
        AIR_HEAT_CAPACITY (float): Heat capacity of air (J/kgK).
        
    Methods:
        materials_list(): Returns the list of materials in the configuration file.
        materials_dict(): Returns a dictionary of materials and their properties.
        reset(): Resets configuration parameters to default values.
        info(): Prints the current configuration parameters.
        to_dict(): Returns the current configuration parameters as a dictionary.
    """
    def __init__(self):
        self.reset()
        self.file = "materials.ini"   # Default configuration file path    
        
    def reset(self):
        self.__La = 2.5
        self.__Nx = 200
        self.__ho = 13
        self.__hi = 8.6
        self.__dt = 10
        self.__AIR_DENSITY = 1.1797660470258469
        self.__AIR_HEAT_CAPACITY = 1005.458757

        self.version = 0
        
    def info(self):
        print("<enerhabitat.Config -- Current config Parameters>")
        print(f"Materials file: \t\t\t{self.__materials_file}")
        print(f"La (Length of dummy frame): \t\t{self.La} m")
        print(f"Nx (Number of discretization elements):\t{self.Nx}")
        print(f"ho (Outdoor convective coefficient): \t{self.ho} W/m²K")
        print(f"hi (Indoor convective coefficient): \t{self.hi} W/m²K")
        print(f"dt (Time step): \t\t\t{self.dt} seconds")
        print(f"\nAIR_DENSITY: \t\t\t\t{self.AIR_DENSITY} kg/m³")
        print(f"AIR_HEAT_CAPACITY: \t\t\t{self.AIR_HEAT_CAPACITY} J/kgK")
    
    def to_dict(self):
        return {
            "La": self.La,
            "Nx": self.Nx,
            "ho": self.ho,
            "hi": self.hi,
            "dt": self.dt,
            "AIR_DENSITY": self.AIR_DENSITY,
            "AIR_HEAT_CAPACITY": self.AIR_HEAT_CAPACITY,
        }

    def materials_list(self):
        """
        Returns the list of materials contained in the configuration file

        Returns:
            list: List of materials in the configuration file
        """
        list_materials = list(self.materials.keys())
        return list_materials
    
    def materials_dict(self):
        new_dict = self.materials.copy()
        for material_i in new_dict.keys():
            new_dict[material_i] = self.materials[material_i].to_dict()
        return new_dict
    
    @property
    def file(self):
        try:    
            # Verificar si el archivo existe
            if not os.path.isfile(self.__materials_file):
                raise FileNotFoundError()
            return self.__materials_file
        except FileNotFoundError:
            print(f"Error: {self.__materials_file} not found")
            
    @file.setter
    def file(self, new_file):
        try:    
            # Verificar si el archivo existe
            if not os.path.isfile(new_file):
                raise FileNotFoundError()

            # Actualizar la configuración global si se proporcionó una nueva ruta
            self.__materials_file = new_file
            
            # Leer el .ini y obtener los nuevos materiales
            new_materials_dict = {}
            materials_data = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
            materials_data.read(self.file)

            for material_i in materials_data.sections():
                k = float(materials_data[material_i]['k'])
                rho = float(materials_data[material_i]['rho'])
                c = float(materials_data[material_i]['c'])
                new_materials_dict[material_i] = Material(k, rho, c) 
            
            self.__materials_class = new_materials_dict
            
        except FileNotFoundError:
            print(f"Error: {new_file} not found")    
    
    @property
    def materials(self):
        return self.__materials_class
    @materials.setter
    def materials(self, value):
        pass
    
    @property
    def La(self):
        return self.__La
    @La.setter
    def La(self, value):
        self.__La = value
        self.version += 1
        
    @property
    def Nx(self):
        return self.__Nx
    @Nx.setter
    def Nx(self, value):
        self.__Nx = value
        self.version += 1
        
    @property
    def ho(self):
        return self.__ho
    @ho.setter
    def ho(self, value):
        self.__ho = value
        self.version += 1
        
    @property
    def hi(self):
        return self.__hi
    @hi.setter
    def hi(self, value):
        self.__hi = value
        self.version += 1
        
    @property
    def dt(self):
        return self.__dt
    @dt.setter
    def dt(self, value):
        # dt is fixed at 10 s and is NOT configurable: the indoor-air node is
        # integrated with an explicit step whose stability requires
        # Fo = hi*dt/(rho_air*c_air*La) < 1. At 10 s, Fo ≈ 0.03. The assignment
        # is ignored to avoid non-physical results.
        pass
        
    @property
    def AIR_DENSITY(self):
        return self.__AIR_DENSITY
    @AIR_DENSITY.setter
    def AIR_DENSITY(self, value):
        self.__AIR_DENSITY = value
        self.version += 1
        
    @property
    def AIR_HEAT_CAPACITY(self):
        return self.__AIR_HEAT_CAPACITY
    @AIR_HEAT_CAPACITY.setter
    def AIR_HEAT_CAPACITY(self, value):
        self.__AIR_HEAT_CAPACITY = value
        self.version += 1
    
# Global configuration instance
config = Config()


class Config2D:
    """
    Extra parameters for the 2D solver (joist and filler block). Does not touch
    the 1D side: the rest of the configuration (La, ho, hi, dt, air) is reused
    from ``config``.

    Attributes:
        nx (int): nodes across the cell width (x direction, adiabatic sides).
        ny (int): nodes through the thickness (y direction, outside→inside).
        tol_inner (float): inner-loop tolerance (line-by-line Gauss-Seidel).
        tol_day (float): day-to-day convergence tolerance.
        max_days (int): cap on the convergence iterations.
    """
    def __init__(self):
        self.reset()

    def reset(self):
        self.nx = 80
        self.ny = 160
        self.tol_inner = 1e-10
        self.tol_day = 5e-4
        self.max_days = 60
        # Engine: serial by default. The parallel one (numba prange over rows) is an
        # OPTION (parallel=True): it spreads the rows across the auto-detected cores,
        # but the line-by-line sweep is fine-grained (thousands of prange regions with
        # a barrier) -> it only pays off ~1.3x on fine meshes (>=~150^2) and is SLOWER
        # on small meshes due to thread overhead (see Phase 7). For real volume
        # (sweeping many configurations) prefer parallelizing at the process level,
        # each solve serial. Verified: parallel reproduces serial bit-for-bit.
        self.parallel = False
        self.version = 0

    def info(self):
        print("<enerhabitat.Config2D -- 2D solver parameters>")
        print(f"nx (width):   \t{self.nx}")
        print(f"ny (thickness):\t{self.ny}")
        print(f"tol_inner:   \t{self.tol_inner}")
        print(f"tol_day:     \t{self.tol_day}")
        print(f"max_days:    \t{self.max_days}")
        print(f"parallel:    \t{self.parallel}")

    def to_dict(self):
        return {"nx": self.nx, "ny": self.ny, "tol_inner": self.tol_inner,
                "tol_day": self.tol_day, "max_days": self.max_days,
                "parallel": self.parallel}


# Global 2D configuration instance
config2d = Config2D()
