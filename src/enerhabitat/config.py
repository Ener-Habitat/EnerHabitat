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
        ho (float): Outdoor convective coefficient (W/(m²·K)).
        hi (float): Indoor convective coefficient, vertical surfaces (W/(m²·K)).
        hi_up (float): Indoor coefficient for roofs with **upward** heat flow
            (indoor surface colder than the indoor air; NOM-020: 9.4 W/(m²·K)).
        hi_down (float): Indoor coefficient for roofs with **downward** heat
            flow (indoor surface warmer than the indoor air, stable
            stratification; NOM-020: 6.6 W/(m²·K)).
        hi_flow (bool): If True (default), systems with ``tilt == 0`` (roofs)
            select ``hi_up``/``hi_down`` at every time step from the direction
            of the heat flow at the indoor surface; walls always use ``hi``.
            Set to False to force the fixed ``hi`` on every orientation.
        dt (float): Time step (seconds). Fixed at 10 s, not configurable.
        AIR_DENSITY (float): Density of air (kg/m³).
        AIR_HEAT_CAPACITY (float): Heat capacity of air (J/(kg·K)).
        
    Methods:
        materials_list(): Returns the list of materials in the configuration file.
        materials_dict(): Returns a dictionary of materials and their properties.
        reset(): Restores the numeric defaults (file/materials are kept).
        info(): Prints the current configuration parameters.
        to_dict(): Returns the current configuration parameters as a dictionary.
    """
    def __init__(self):
        self.reset()
        # Default materials file: loaded silently if present in the working
        # directory; otherwise `materials` starts empty ({}) and `config.file`
        # must be set explicitly before solving.
        self.__materials_file = "materials.ini"
        self.__materials_class = {}
        if os.path.isfile(self.__materials_file):
            self.file = self.__materials_file


    def reset(self):
        """Restore the **numeric** defaults (La, Nx, ho, hi, dt, air
        properties). ``file``/``materials`` are user-provided and are kept."""
        self.__La = 2.5
        self.__Nx = 200
        self.__ho = 13
        self.__hi = 8.1
        self.__hi_up = 9.4
        self.__hi_down = 6.6
        self.__hi_flow = True
        self.__dt = 10
        self.__AIR_DENSITY = 1.1797660470258469
        self.__AIR_HEAT_CAPACITY = 1005.458757

        self.version = 0
        
    def info(self):
        print("<enerhabitat.Config -- Current config Parameters>")
        print(f"Materials file: \t\t\t{self.__materials_file}")
        print(f"La (Length of dummy frame): \t\t{self.La} m")
        print(f"Nx (Number of discretization elements):\t{self.Nx}")
        print(f"ho (Outdoor convective coefficient): \t{self.ho} W/(m²·K)")
        print(f"hi (Indoor convective coefficient): \t{self.hi} W/(m²·K)")
        print(f"hi_up / hi_down (roofs, by heat flow): \t{self.hi_up} / {self.hi_down} W/(m²·K)"
              f"  [hi_flow={'on' if self.hi_flow else 'off'}]")
        print(f"dt (Time step): \t\t\t{self.dt} seconds")
        print(f"\nAIR_DENSITY: \t\t\t\t{self.AIR_DENSITY} kg/m³")
        print(f"AIR_HEAT_CAPACITY: \t\t\t{self.AIR_HEAT_CAPACITY} J/(kg·K)")
    
    def to_dict(self):
        return {
            "La": self.La,
            "Nx": self.Nx,
            "ho": self.ho,
            "hi": self.hi,
            "hi_up": self.hi_up,
            "hi_down": self.hi_down,
            "hi_flow": self.hi_flow,
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
        """Path of the configured materials file (no I/O; may not exist if the
        default was never loaded)."""
        return self.__materials_file

    @file.setter
    def file(self, new_file):
        """Load the materials from ``new_file``.

        Raises:
            FileNotFoundError: if ``new_file`` does not exist. The previously
                loaded materials (and path) are kept untouched.
        """
        if not os.path.isfile(new_file):
            raise FileNotFoundError(
                f"materials file not found: {new_file!r} "
                f"(previously loaded materials are kept)")

        # Read the .ini and load the new materials
        new_materials_dict = {}
        materials_data = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
        materials_data.read(new_file)

        for material_i in materials_data.sections():
            k = float(materials_data[material_i]['k'])
            rho = float(materials_data[material_i]['rho'])
            c = float(materials_data[material_i]['c'])
            for prop, val in (("k", k), ("rho", rho), ("c", c)):
                if not (val > 0):
                    raise ValueError(
                        f"material {material_i!r} in {new_file!r}: "
                        f"{prop} must be > 0, got {val}")
            new_materials_dict[material_i] = Material(k, rho, c)

        # Commit only after a successful parse (rollback semantics)
        self.__materials_file = new_file
        self.__materials_class = new_materials_dict
    
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
        if not (value > 0):
            raise ValueError(f"La must be > 0 m, got {value!r}")
        self.__La = value
        self.version += 1

    @property
    def Nx(self):
        return self.__Nx
    @Nx.setter
    def Nx(self, value):
        if int(value) != value or int(value) < 3:
            raise ValueError(f"Nx must be an integer >= 3, got {value!r}")
        self.__Nx = int(value)
        self.version += 1

    @property
    def ho(self):
        return self.__ho
    @ho.setter
    def ho(self, value):
        if not (value > 0):
            raise ValueError(f"ho must be > 0 W/(m2K), got {value!r}")
        self.__ho = value
        self.version += 1
        
    @property
    def hi(self):
        return self.__hi
    @hi.setter
    def hi(self, value):
        if not (value > 0):
            raise ValueError(f"hi must be > 0 W/(m2K), got {value!r}")
        self.__hi = value
        self.version += 1
        
    @property
    def hi_up(self):
        return self.__hi_up
    @hi_up.setter
    def hi_up(self, value):
        if not (value > 0):
            raise ValueError(f"hi_up must be > 0 W/(m2K), got {value!r}")
        self.__hi_up = value
        self.version += 1

    @property
    def hi_down(self):
        return self.__hi_down
    @hi_down.setter
    def hi_down(self, value):
        if not (value > 0):
            raise ValueError(f"hi_down must be > 0 W/(m2K), got {value!r}")
        self.__hi_down = value
        self.version += 1

    @property
    def hi_flow(self):
        return self.__hi_flow
    @hi_flow.setter
    def hi_flow(self, value):
        self.__hi_flow = bool(value)
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
        if not (value > 0):
            raise ValueError(f"AIR_DENSITY must be > 0 kg/m3, got {value!r}")
        self.__AIR_DENSITY = value
        self.version += 1
        
    @property
    def AIR_HEAT_CAPACITY(self):
        return self.__AIR_HEAT_CAPACITY
    @AIR_HEAT_CAPACITY.setter
    def AIR_HEAT_CAPACITY(self, value):
        if not (value > 0):
            raise ValueError(
                f"AIR_HEAT_CAPACITY must be > 0 J/(kgK), got {value!r}")
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
        tol_inner (float): inner-loop tolerance in °C. A time step is accepted
            when BOTH the largest node update of the last sweep and the largest
            scaled residual of the discrete equations (|a_P·T_P − Σa_nb·T_nb −
            b| / a_P) fall below it. 1e-8 °C by default: with 8640 steps/day the
            worst-case accumulated error (8640 × tol_inner ≈ 1e-4 °C) stays
            below ``tol_day``.
        tol_day (float): day-to-day convergence tolerance (°C).
        max_days (int): cap on the day-to-day convergence iterations.
        max_inner (int): cap on inner sweeps per time step; exceeded → the
            solve is flagged ``converged = False``.
    """
    def __init__(self):
        self.reset()

    def reset(self):
        self.nx = 80
        self.ny = 160
        self.tol_inner = 1e-8
        self.tol_day = 5e-4
        self.max_days = 60
        self.max_inner = 10000
        self.version = 0

    def validate(self):
        """Raise ``ValueError`` if any 2D numeric parameter is out of range.
        Called by ``System2D`` at the start of every solve (the attributes are
        plain, so assignment itself is not validated)."""
        for name in ("nx", "ny"):
            v = getattr(self, name)
            if int(v) != v or int(v) < 3:
                raise ValueError(f"config2d.{name} must be an integer >= 3, "
                                 f"got {v!r}")
        for name in ("tol_inner", "tol_day"):
            v = getattr(self, name)
            if not (v > 0):
                raise ValueError(f"config2d.{name} must be > 0, got {v!r}")
        for name in ("max_days", "max_inner"):
            v = getattr(self, name)
            if int(v) != v or int(v) < 1:
                raise ValueError(f"config2d.{name} must be an integer >= 1, "
                                 f"got {v!r}")

    def info(self):
        print("<enerhabitat.Config2D -- 2D solver parameters>")
        print(f"nx (width):   \t{self.nx}")
        print(f"ny (thickness):\t{self.ny}")
        print(f"tol_inner:   \t{self.tol_inner}")
        print(f"tol_day:     \t{self.tol_day}")
        print(f"max_days:    \t{self.max_days}")
        print(f"max_inner:   \t{self.max_inner}")

    def to_dict(self):
        return {"nx": self.nx, "ny": self.ny, "tol_inner": self.tol_inner,
                "tol_day": self.tol_day, "max_days": self.max_days,
                "max_inner": self.max_inner}


# Global 2D configuration instance
config2d = Config2D()
