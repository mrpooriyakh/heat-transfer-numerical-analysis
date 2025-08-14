# Heat Transfer Numerical Analysis

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![NumPy](https://img.shields.io/badge/NumPy-1.21+-orange.svg)](https://numpy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.4+-green.svg)](https://matplotlib.org/)

A comprehensive numerical analysis toolkit for solving heat transfer problems using finite difference methods. This project implements solutions for 1D, cylindrical, spherical, and 2D heat transfer scenarios with both conduction and convection boundary conditions.

## 🔥 Features

### 🔢 1D Heat Transfer Analysis
- **Bar Conduction**: Temperature distribution in bars with fixed temperature boundaries
- **Bar Convection**: Heat transfer with convective boundary conditions
- Analytical vs. numerical solution comparison
- Error analysis and convergence studies

### 🌀 Cylindrical Coordinate Systems
- **Cylinder Conduction**: Radial heat transfer with temperature boundary conditions
- **Cylinder Convection**: Convective heat transfer in cylindrical geometries
- Heat transfer rate calculations
- High-precision numerical solutions

### 🔮 Spherical Coordinate Systems
- **Sphere Conduction**: Radial heat transfer in spherical geometries
- **Sphere Convection**: Convective boundary conditions on spherical surfaces
- Analytical validation and error quantification
- Heat transfer rate computations

### 📐 2D Heat Transfer
- **2D Plate Conduction**: Temperature distribution in rectangular plates
- Multiple boundary condition support
- 3D visualization capabilities
- Contour plots and cross-sectional analysis

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/mrpooriyakh/heat-transfer-numerical-analysis.git
cd heat-transfer-numerical-analysis

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```python
# 1D Bar Heat Transfer
from src.one_dimensional.bar_conduction import solve_bar_heat_transfer

# Solve 1D conduction problem
x, T_numerical, T_analytical = solve_bar_heat_transfer(
    L=0.2,      # Length (m)
    T0=120,     # Temperature at x=0 (°C)
    TL=50,      # Temperature at x=L (°C)
    k=1.2,      # Thermal conductivity (W/m°C)
    n=500       # Number of nodes
)

# Cylindrical Heat Transfer
from src.cylindrical.cylinder_conduction import solve_cylinder_heat_transfer

radius, temp, T_analytical, error, Q, Q_analytical = solve_cylinder_heat_transfer(
    L=4,        # Length (m)
    r1=0.01,    # Inner radius (m)
    r2=0.05,    # Outer radius (m)
    k=10,       # Thermal conductivity (W/m°C)
    bc1=100,    # Inner temperature (°C)
    bc2=40,     # Outer temperature (°C)
    n=400       # Number of points
)
```

## 📊 Mathematical Formulation

### 1D Heat Conduction Equation
```
d²T/dx² = 0
```

With boundary conditions:
- T(0) = T₀
- T(L) = T_L

### Cylindrical Coordinates
```
d/dr(r dT/dr) = 0
```

### Spherical Coordinates
```
d/dr(r² dT/dr) = 0
```

### 2D Heat Equation
```
∂²T/∂x² + ∂²T/∂y² = 0
```

## 📈 Numerical Methods

### Finite Difference Discretization
- **Central Difference**: Second-order accurate spatial discretization
- **Matrix Solution**: Direct solving using numpy.linalg.solve()
- **Boundary Condition Handling**: Incorporation of Dirichlet and Neumann conditions
- **Convergence Analysis**: Grid refinement studies and error quantification

### Solution Techniques
- **Tridiagonal Systems**: Efficient solution for 1D problems
- **Sparse Matrix Methods**: Optimized for large 2D systems
- **Analytical Validation**: Comparison with exact solutions where available

## 🗂️ Project Structure

```
src/
├── one_dimensional/          # 1D heat transfer problems
│   ├── bar_conduction.py     # Fixed temperature boundaries
│   └── bar_convection.py     # Convective boundaries
├── cylindrical/              # Cylindrical coordinate problems
│   ├── cylinder_conduction.py    # Radial conduction
│   └── cylinder_convection.py    # Convective boundaries
├── spherical/                # Spherical coordinate problems
│   ├── sphere_conduction.py      # Radial conduction
│   └── sphere_convection.py      # Convective boundaries
├── two_dimensional/          # 2D heat transfer
│   └── plate_2d_conduction.py    # 2D rectangular plates
└── utils/                    # Common utilities
    ├── plotting.py           # Visualization functions
    └── numerical_methods.py  # Common numerical tools
```

## 📋 Examples and Results

### 1D Bar Heat Transfer
- **Problem**: Bar with T(0)=120°C, T(L)=50°C
- **Results**: Linear temperature distribution
- **Validation**: Perfect agreement with analytical solution
- **Error**: < 0.001% for n=500 nodes

### Cylindrical Heat Transfer
- **Problem**: Hollow cylinder with inner/outer temperature boundaries
- **Results**: Logarithmic temperature profile
- **Heat Transfer Rate**: Calculated numerically and analytically
- **Error**: < 0.001% for n=400 points

### Spherical Heat Transfer
- **Problem**: Hollow sphere with temperature boundary conditions
- **Results**: Hyperbolic temperature distribution
- **Validation**: Excellent agreement with analytical solution
- **Error**: < 0.001% for high grid resolution

### 2D Plate Heat Transfer
- **Problem**: Rectangular plate with different boundary temperatures
- **Results**: 2D temperature field with smooth transitions
- **Visualization**: 3D surface plots and contour maps
- **Cross-sections**: Temperature profiles along different directions

## 🔧 Key Parameters and Configuration

### Computational Parameters
- **Grid Resolution**: Adjustable number of nodes/points
- **Convergence Tolerance**: Configurable solution accuracy
- **Boundary Conditions**: Support for multiple BC types
- **Material Properties**: Thermal conductivity specification

### Physical Parameters
- **Geometry**: Flexible domain dimensions
- **Thermal Properties**: Conductivity, convection coefficients
- **Boundary Temperatures**: Fixed or convective conditions
- **Heat Transfer Rates**: Automatic calculation and validation

## 📚 Validation and Verification

### Analytical Comparisons
- **1D Problems**: Linear profiles validated against exact solutions
- **Cylindrical**: Logarithmic profiles compared with analytical results
- **Spherical**: Hyperbolic distributions verified analytically
- **2D Cases**: Cross-sectional comparisons with known solutions

### Error Analysis
- **Grid Convergence**: Systematic refinement studies
- **Truncation Error**: Order of accuracy verification
- **Round-off Error**: Numerical precision analysis
- **Relative Errors**: Typically < 0.01% for well-resolved grids

## 🛠️ Dependencies

```python
numpy>=1.21.0          # Numerical computations
matplotlib>=3.4.0      # Plotting and visualization
scipy>=1.7.0           # Scientific computing (optional)
```

## 📖 Usage Examples

### Complete 1D Analysis
```python
from src.one_dimensional.bar_conduction import solve_bar_heat_transfer
import matplotlib.pyplot as plt

# Define problem parameters
L = 0.2      # Length (m)
T0 = 120     # Left boundary temperature (°C)
TL = 50      # Right boundary temperature (°C)
k = 1.2      # Thermal conductivity (W/m°C)
n = 500      # Number of internal nodes

# Solve the problem
x, T_num, T_ana = solve_bar_heat_transfer(L, T0, TL, k, n)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(x, T_num, 'b-', label='Numerical')
plt.plot(x, T_ana, 'r--', label='Analytical')
plt.xlabel('Distance (m)')
plt.ylabel('Temperature (°C)')
plt.legend()
plt.grid(True)
plt.show()
```

### Cylindrical Heat Transfer Analysis
```python
from src.cylindrical.cylinder_convection import solve_cylinder_heat_transfer_convective

# Define parameters
L = 4.0      # Length (m)
r1 = 0.01    # Inner radius (m)
r2 = 0.05    # Outer radius (m)
k = 10       # Thermal conductivity (W/m°C)
ta1 = 250    # Inner gas temperature (°C)
ha1 = 8.17   # Inner convection coefficient (W/m²°C)
ta2 = 35     # Outer gas temperature (°C)
ha2 = 24     # Outer convection coefficient (W/m²°C)
n = 1600     # Number of points

# Solve and visualize
radius, temp, T_analytical, error, Q, Q_analytical = solve_cylinder_heat_transfer_convective(
    L, r1, r2, k, ta1, ha1, ta2, ha2, n)

print(f"Heat Transfer Rate: {Q:.2f} W")
print(f"Maximum Error: {np.max(error):.6f}%")
```

## 🔬 Research Applications

### Academic Use
- **Heat Transfer Courses**: Numerical methods in thermal engineering
- **Computational Methods**: Finite difference implementation
- **Validation Studies**: Benchmark problem solutions
- **Grid Convergence**: Numerical accuracy analysis

### Industrial Applications
- **Thermal Design**: Component temperature prediction
- **Heat Exchanger Analysis**: Performance optimization
- **Thermal Management**: Electronic cooling applications
- **Process Engineering**: Heat transfer equipment design

## 🤝 Contributing

Contributions are welcome! Please follow these guidelines:

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/new-method`
3. **Add your contributions**: New methods, improvements, documentation
4. **Test thoroughly**: Validate against analytical solutions
5. **Update documentation**: Add examples and explanations
6. **Submit a pull request**: Describe your changes clearly

### Development Guidelines
- Follow PEP 8 style conventions
- Add comprehensive docstrings
- Include validation against analytical solutions
- Provide usage examples
- Update README with new features

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Numerical Methods Community** - For established finite difference techniques
- **Heat Transfer Research** - Classical analytical solutions for validation
- **Scientific Python Ecosystem** - NumPy, SciPy, and Matplotlib developers
- **Academic Community** - For continuous advancement in computational heat transfer

## 📞 Contact

**Project Maintainer**: [Your Name]
- GitHub: [@mrpooriyakh](https://github.com/mrpooriyakh)
- Email: [your.email@example.com]

**Project Repository**: [https://github.com/mrpooriyakh/heat-transfer-numerical-analysis](https://github.com/mrpooriyakh/heat-transfer-numerical-analysis)

---

*This project demonstrates comprehensive numerical solutions for heat transfer problems across multiple coordinate systems and boundary condition types.*
