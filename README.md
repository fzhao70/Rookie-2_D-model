# Rookie 2D Atmospheric Model

A simple 2D atmospheric model for undergraduate atmospheric modeling education. This project implements shallow water equations with geostrophic wind initialization, providing a practical introduction to numerical weather prediction concepts.

## Overview

This model simulates atmospheric dynamics on a 2D horizontal plane using the shallow water equations. It includes geostrophic wind calculations and time integration using both Matsuno and Leapfrog schemes. The code is designed for educational purposes, making it accessible for students learning computational atmospheric modeling.

## Features

- **Shallow Water Equations**: Full implementation of 2D shallow water equations including advection, pressure gradient force, and Coriolis effect
- **Multiple Time Integration Schemes**:
  - Matsuno scheme (for initialization)
  - Leapfrog scheme (for time integration)
- **Geostrophic Wind Initialization**: Initialize wind fields from geopotential height using geostrophic balance
- **Finite Difference Methods**: Implements forward, backward, and centered difference schemes with proper boundary handling
- **NetCDF I/O**: Reads and writes NetCDF format files for easy visualization and analysis
- **OpenMP Parallelization**: Optional parallel processing for improved performance
- **Modular Design**: Clean separation of utilities, solvers, and I/O operations

## Physical Model

### Governing Equations

The model solves the shallow water equations on a 2D plane:

```
∂u/∂t = -u·∂u/∂x - v·∂u/∂y - g·∂Z/∂x + f·v
∂v/∂t = -u·∂v/∂x - v·∂v/∂y - g·∂Z/∂y - f·u
∂Z/∂t = -u·∂Z/∂x - v·∂Z/∂y - H·(∂u/∂x + ∂v/∂y)
```

Where:
- `u`, `v`: Horizontal wind components (m/s)
- `Z`: Geopotential height (m)
- `g`: Gravitational acceleration (9.8 m/s²)
- `f`: Coriolis parameter (2Ω·sin(φ))
- `H`: Scale height (8000 m)
- `Ω`: Earth's angular velocity (7.292×10⁻⁵ rad/s)

### Geostrophic Wind Initialization

Wind fields are initialized using geostrophic balance:

```
u_g = -(g/f)·∂Z/∂y
v_g = (g/f)·∂Z/∂x
```

## Requirements

### Software Dependencies

- **Fortran Compiler**: gfortran (recommended) or Intel Fortran (ifort)
- **NetCDF Fortran Library**: netcdf-fortran (version 4.x or later)
- **Make**: GNU Make for building the project
- **OpenMP Support**: Optional, for parallel execution

### Installing Dependencies

#### Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install gfortran libnetcdff-dev make
```

#### macOS (using Homebrew):
```bash
brew install gcc netcdf netcdf-fortran
```

#### CentOS/RHEL:
```bash
sudo yum install gcc-gfortran netcdf-fortran-devel make
```

## Installation

1. **Clone the repository**:
```bash
git clone https://github.com/fzhao70/Rookie-2_D-model.git
cd Rookie-2_D-model
```

2. **Build the project**:
```bash
make
```

This will create two executables:
- `simple_model`: The main 2D atmospheric model
- `geo_wind`: Geostrophic wind calculation utility

### Build Options

- **Build with Intel Fortran**:
  ```bash
  make FC=ifort
  ```

- **Debug build** (with additional checks and debugging symbols):
  ```bash
  make debug
  ```

- **Clean build artifacts**:
  ```bash
  make clean      # Remove executables and object files
  make cleanall   # Also remove output NetCDF files
  ```

- **View help**:
  ```bash
  make help
  ```

## Usage

### 1. Geostrophic Wind Calculation

Calculate geostrophic winds from geopotential height:

```bash
./geo_wind
```

**Input**: `hgt_location_selected.nc` (geopotential height field)
**Output**: `geo_wind.nc` (u and v wind components)

### 2. Simple 2D Model

Run the full atmospheric model simulation:

```bash
./simple_model
```

**Input**: `hgt_location_selected.nc` (geopotential height for initialization)
**Output**: `simple_model_output.nc` (final state: u, v, Z fields)

### Configuration

Model parameters can be adjusted in `simple_model_type.f90`:

- `nx, ny`: Grid dimensions (default: 33×23)
- `nt`: Number of time steps (default: 1000)
- `dt`: Time step in seconds (default: 100.0)

## Input/Output File Format

### Input File: `hgt_location_selected.nc`

NetCDF file containing:
- `hgt(lat, lon)`: Geopotential height field (m)
- `lat(ny)`: Latitude coordinates (degrees or radians)
- `lon(nx)`: Longitude coordinates (degrees or radians)

### Output Files

**`geo_wind.nc`**:
- `u(lat, lon)`: Zonal wind component (m/s)
- `v(lat, lon)`: Meridional wind component (m/s)
- `lat(ny)`: Latitude coordinates
- `lon(nx)`: Longitude coordinates

**`simple_model_output.nc`**:
- `u(lat, lon)`: Final zonal wind field (m/s)
- `v(lat, lon)`: Final meridional wind field (m/s)
- `Z(lat, lon)`: Final geopotential height field (m)
- `lat(ny)`: Latitude coordinates
- `lon(nx)`: Longitude coordinates

## Project Structure

```
Rookie-2_D-model/
├── README.md                  # This file
├── LICENSE                    # MIT License
├── Makefile                   # Build system
├── simple_model_type.f90     # Main model with modules:
│                              #   - type_def: Data structures
│                              #   - utils: Utility functions
│                              #   - solver: Time integration schemes
│                              #   - io: NetCDF input/output
└── geo_wind.f90              # Geostrophic wind calculation program
```

## Code Modules

### `type_def` Module
Defines data structures for grid and meteorological elements.

### `utils` Module
- `message_display()`: Formatted message output
- `distance()`: Haversine distance calculation between coordinates
- `centered_diff()`: Centered finite difference
- `forward_diff()`: Forward finite difference
- `backward_diff()`: Backward finite difference
- `horizontal_diff()`: 2D horizontal derivative with boundary handling

### `solver` Module
- `Matsuno()`: Matsuno time integration scheme
- `Leapfrog()`: Leapfrog time integration scheme

### `io` Module
- `read_input()`: Read NetCDF input files
- `write_output()`: Write NetCDF output files

## Algorithm Details

### Time Integration

1. **First time step**: Uses the Matsuno scheme (a predictor-corrector method) for stability
2. **Subsequent steps**: Uses the Leapfrog scheme for computational efficiency

### Spatial Discretization

- Uses centered differences for interior points
- Applies forward/backward differences at boundaries
- Calculates distances using the Haversine formula for spherical coordinates

## Example Workflow

1. Prepare input data with geopotential height field
2. Run geostrophic wind calculation:
   ```bash
   ./geo_wind
   ```
3. Run the full model simulation:
   ```bash
   ./simple_model
   ```
4. Visualize results using Python/NCL/MATLAB:
   ```python
   import netCDF4 as nc
   import matplotlib.pyplot as plt

   data = nc.Dataset('simple_model_output.nc')
   u = data.variables['u'][:]
   v = data.variables['v'][:]

   # Create wind vector plot
   plt.quiver(u, v)
   plt.show()
   ```

## Performance Notes

- The code includes OpenMP directives for parallel execution (`!$OMP PARALLEL DO`)
- Set `OMP_NUM_THREADS` environment variable to control parallelization:
  ```bash
  export OMP_NUM_THREADS=4
  ./simple_model
  ```

## Limitations

- 2D model (no vertical structure)
- Simplified physics (shallow water approximation)
- No moisture or thermodynamics
- Periodic or fixed boundary conditions
- Educational tool, not for operational forecasting

## Troubleshooting

### NetCDF errors
If you encounter NetCDF library errors, ensure:
1. NetCDF Fortran is properly installed
2. The `nf-config` command is in your PATH
3. Input files exist and are valid NetCDF format

### Compilation errors
- Verify compiler compatibility (tested with gfortran 7.x and later)
- Check NetCDF library paths in the Makefile
- Ensure OpenMP support is available if using parallel features

## Educational Use

This model is designed for teaching:
- Numerical methods for PDEs
- Atmospheric dynamics concepts
- Time integration schemes
- Finite difference methods
- Scientific computing in Fortran
- NetCDF data handling

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes with clear commit messages
4. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Fanghe Zhao**
- Email: zfh1997@mail.ustc.edu.cn
- Institution: USTC-AEMOL (University of Science and Technology of China - Atmospheric and Environmental Modeling Laboratory)

## Acknowledgments

Developed for undergraduate atmospheric modeling education at USTC.

## References

For theoretical background on shallow water equations and numerical methods:
1. Holton, J. R., & Hakim, G. J. (2012). *An Introduction to Dynamic Meteorology* (5th ed.). Academic Press.
2. Durran, D. R. (2010). *Numerical Methods for Fluid Dynamics*. Springer.
3. Mesinger, F., & Arakawa, A. (1976). *Numerical Methods Used in Atmospheric Models*. GARP Publications Series No. 17.

## Citation

If you use this code in your research or teaching, please cite:
```
Zhao, F. (2018). Rookie 2D Atmospheric Model.
GitHub repository: https://github.com/fzhao70/Rookie-2_D-model
```
