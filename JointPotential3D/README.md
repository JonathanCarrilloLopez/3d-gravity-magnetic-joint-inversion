# JointPotential3D

JointPotential3D is a Fortran90 software package for 3D forward and inverse modeling of gravity and magnetic data, including joint inversion through correspondence maps.

The software implements and extends the general joint inversion framework proposed by Carrillo & Gallardo (2018) for integrated interpretation of potential field data in geothermal environments.

## Overview

The software was developed for integrated geophysical characterization of complex subsurface systems, with emphasis on geothermal exploration.

JointPotential3D allows:

- Gravity forward modeling
- Magnetic forward modeling
- Separate inversion of gravity or magnetic data
- Joint inversion of gravity and magnetic data
- Estimation of correspondence relationships between physical properties
- Integrated density and magnetization modeling in 3D

The methodology was applied to geothermal systems such as Los Humeros and Acoculco (Mexico) within the GEMex project.

## Main Features

- 3D prism-based modeling
- Joint inversion using correspondence maps
- Support for prior models and covariance matrices
- Synthetic and real-data applications
- Configurable smoothness constraints
- Gravity and magnetic data integration
- Linux/Unix command-line execution

## Scientific Background

The core of the software is based on a joint inversion methodology that estimates relationships between density and magnetization distributions through correspondence maps, reducing ambiguities commonly associated with independent inversions of potential field data.

Related publication:

Carrillo-Lopez, J., Perez-Flores, M. A., Gallardo, L. A., & Schill, E. (2022).

Joint inversion of gravity and magnetic data using correspondence maps with application to geothermal fields.

Geophysical Journal International, 228(3), 1621–1636.

https://doi.org/10.1093/gji/ggab416

## Compilation

The source code is contained in a single Fortran90 file:

JointPotential3D.f90

Compile using a Fortran compiler such as:

gfortran JointPotential3D.f90 -o InvConj3D

## Running the Program

Execute from the command line:

./InvConj3D

The program automatically reads all required input files from the current working directory.

## Input Files

### Configuration
- Startup.dat

### Geophysical Data
- Gravity.dat
- Magnetic.dat
- Cd1.dat
- Cdd2.dat

### Coordinates
- X_UTM_grav.dat
- Y_UTM_grav.dat
- Z_UTM_grav.dat
- X_UTM_mag.dat
- Y_UTM_mag.dat
- Z_UTM_mag.dat

### Prior Models
- m01pr_density.dat
- m02pr_magnetization.dat

### Covariances
- cov_m01pr.dat
- cov_m02pr.dat
- Cgg.dat

### Mesh
- mesh.dat

## Modeling Modes

The software supports:

1. Forward gravity/magnetic modeling
2. Separate inversion
3. Joint inversion

The process is controlled through the Startup.dat configuration file.

## Output Files

Typical outputs include:

- Predicted gravity data
- Predicted magnetic data
- Density models
- Magnetization models
- Correspondence functions
- RMS fitting statistics

## Applications

JointPotential3D has been applied to:

- Geothermal exploration
- Structural interpretation
- Integrated subsurface characterization
- Density-magnetization relationship analysis
- Potential field inversion research

## Author

Jonathan Carrillo Lopez

## Acknowledgments

Part of this development was carried out within the GEMex project (Cooperation in Geothermal Energy Research Europe-Mexico for development of Enhanced Geothermal Systems and Superhot Geothermal Systems), funded by Horizon 2020 and CONACYT-SENER.

## License

Academic and research use.