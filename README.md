# EnerHabitat
Tools for evaluating the thermal performance of structures based on EPW files

## Contents
- [Installation](#installation)
- [Recommended structure](#folder-sturcture)
- [Key concepts](#key-concepts)
- [Quickstart](#quickstart)
- [Main Structures](#main-structures)
  - [Location](#location)
  - [System](#system)
  - [Config](#config-global)
- [Materials](#materials)
- [Dependencies](#dependencies)
- [License](#license)

## Getting started
**enerhabitat** is a Python package for thermal simulation of constructive systems based on meteorological data from EPW files. The values that are calculated include:

- Ta  : Ambient temperature
- Tsa : Sun-air temperature
- Ti  : Interior temperature
- Ig  : Global horizontal irradiance
- Ib  : Direct normal irradiance
- Id  : Diffuse horizontal irradiance
- Is  : Surface irradiance

## Installation

The source code is currently hosted on GitHub at [EnerHabitat](https://github.com/Ener-Habitat/EnerHabitat)

Binary installers for the latest released version are available at the Python Package Index [PyPI](https://pypi.org/project/enerhabitat) 

```bash
pip install enerhabitat
```

If you're working with the [uv](https://docs.astral.sh/uv/) Python package manager you can use the following

```bash
$ uv add enerhabitat
```

### Folder structure

The following shows the basic folder structure recommended
```
├── main.py
├── materials.ini   # Materials properties
└── epw
    ├── ...
    └── example_file.epw
```

## Key concepts
- `Location` reads an EPW file and computes the average day with `meanDay()`.
- `System` uses a `Location` and a list of layers to compute `Tsa`, `solve`, and `solveAC`.
- `config` is a global instance.  Modify `config`'s attributes to change all subsequent computations behavior.

## Quickstart
```python
import enerhabitat as eh

# 1) Global configuration
# eh.config.file = "./materials.ini"

# 2) Location
epw_file = "./epw/example.epw"
loc = eh.Location(epw_file)

# 3) Constructive system
wall = eh.System(location=loc)
wall.layers = [("Adobe", 0.2)]

# 4) Solve
tsa = wall.Tsa()
ti = wall.solve()
```

## Main structures

Load the EnerHabitat package for use
```python
import enerhabitat as eh
```

## Location
Create a location object from an EPW file.

```python
import enerhabitat as eh

epw_file = "./epw/example.epw"
loc = eh.Location(epw_file)
```

### Attributes

The EPW path is stored in`file`, other Location read-only values include:
- `city` : `str` with the city of the EPW header.
- `latitude` : `float` value in degrees.
- `longitude` : `float` in degrees.
- `altitude` : `float` value in meters.
- `timezone` : `pytz.timezone` object.

Read-only attributes are recovered from the `file` header, assign a new EPW file to `file` to modify them

```python
loc.file = "./epw/other.epw"
```

### Methods
- `meanDay()` Calculates the average-day DataFrame with Ta, Ig, Ib, Id, and Tn.
- `copy()` Returns a copy of the instance
- `info()` Prints instance attributes
- `flag()` Dict with the last `meanDay()` computation values 

Examples:
```python
loc.meanDay(month="6", year=2020).info()
loc.info()
print(loc.flag()["date"])
```

## System
Models a constructive system and computes interior temperatures. Crate a `System` from a `Location` object

```python
import enerhabitat as eh

loc = eh.Location("./epw/example.epw")
wall =  eh.System(location = loc)
```

### Attributes
- `location` : Associated `Location` object
- `tilt` : `float`
- `azimuth` : `float`
- `absortance` : `float`
- `layers` : `list` with `(material, width)` tuples
```python
wall.location = loc_2
wall.tilt = 0
wall.azimuth = 45,
wall.absortance = 0.3)

# Layers from outside to inside
wall.layers = [("Adobe", 0.10), ("Acero", 0.05), ("Ladrillo", 0.02)]

# Add and remove layers
wall.add_layer("Mortero", 0.20)
wall.remove_layer(2)
```

Other computation result values are stores as read-only attributes, this include:
- `energy_transfer` : Total energy transfer computed in `solve()` 
- `heating_energy` : Total heating energy need computed in `solveAC()`
- `cooling_energy` : Total cooling energy need computed in `solveAC()`

### Methods
- `Tsa()` computes Tsa and Is based on `Location.meanDay()`
- `solve()` computes Ti for the constructive system's `layers`
- `solveAC()` computes cooling and heating energy with constant Ti
- `copy()` returns a copy of the instance
- `info()` prints attributes
- `flag()` tells whether the cached value was recalculated

Examples:
```python
# Tsa
wall.Tsa().info()

# Computation internal temperature
ti = wall.solve()
energy = wall.energy_transfer

# Computation Air conditioning
wall.solveAC()
c_energy = wall.cooling_energy
h_energy = wall.heating_energy
```


## Config (global)
`config` stores parameters shared by all instances.

Note: this is global. Changing it affects future computations for any `Location` or `System`.

### Attributes:
- `file`: Path to `.ini` file containing material's name and propertys
- `La` : length of the fictional room
- `Nx` : Number of elements to discretize the construction system
- `ho` : Outside convection heat transfer
- `hi` : Inside coonvection heat transfer
- `dt` : Time step inseconds

Example:
```python
import enerhabitat as eh

# Change materials file
eh.config.file = "./materials.ini"

# Update parameters
eh.config.La = 2.0
eh.config.Nx = 300
eh.config.ho = 12
eh.config.hi = 8.3
eh.config.dt = 60
```

Additionally `config` stores materials propertys in the read-only `materials` attribute in the form of a `dict` with the names of materials as keys and their properties `k`, `rho` and `c` class attributes

```python
eh.config.materials

# Prueba con un solo material
adobe = eh.config.materials["Adobe"]

adobe.rho
adobe.c
adobe.k
```

### Methods
- `info()` prints current `config` values
- `to_dict()` returns `config` parameters as a `dict`
- `reset()` restores  default `config` values
- `materials_list()` returns the material names defined in `file`
- `materials_dict()` returns a `dict` of material properties

## Materials
The material's name and their properties are specified in the `materials.ini` configuration file, specifying the material name as the `key` and its values ​​for `k`, `rho` and `c`.

```ini
[concrete]
k   = 1.35  # Overall heat transfer coefficient
rho = 1800  # Density
c   = 1000  # Specific heat

[adobe]
k   = 0.58
rho = 1500
c   = 1480
```

To set a configuration file, use the `config.file` attribute and specify the `path`.

```python
eh.config.file = "./config/new_materials.ini"
```
## Dependencies
- [numba](https://numba.pydata.org/)
- [pvlib](https://pvlib-python.readthedocs.io/en/stable/)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [pytz](https://pypi.org/project/pytz/)

## License
Code released under the [MIT](https://github.com/AltamarMx/enerhabitat/main/LICENSE) license.
