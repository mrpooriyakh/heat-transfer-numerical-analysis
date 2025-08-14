import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 0.2  # Length of bar (m)
T0 = 120  # Temperature at x=0 (°C)
TL = 50   # Temperature at x=L (°C)
k = 1.2   # Thermal conductivity (W/m°C)
A = 15    # Cross-sectional area (m²)
n = 500   # Number of nodes

# Create spatial grid
dx = L/(n+1)
x = np.linspace(0, L, n+2)

# Initialize matrices
a = np.zeros((n, n))
b = np.zeros(n)

# Fill coefficient matrix
for i in range(n):
    if i == 0:
        a[i,i+1] = 1/dx**2
        a[i,i] = -2/dx**2
        b[i] = -T0/dx**2
    elif i == n-1:
        a[i,i-1] = 1/dx**2
        a[i,i] = -2/dx**2
        b[i] = -TL/dx**2
    else:
        a[i,i-1] = 1/dx**2
        a[i,i] = -2/dx**2
        a[i,i+1] = 1/dx**2

# Solve system
T = np.linalg.solve(a, b)

# Combine with boundary conditions
T_numerical = np.concatenate(([T0], T, [TL]))

# Analytical solution
T_analytical = T0 + x*(TL-T0)/L

# Calculate error
error_percentage = abs(T_numerical - T_analytical)/T_analytical * 100

# Calculate heat transfer rate
heat_transfer_rate = -k*A*(TL-T0)/L

# Plotting
plt.figure(figsize=(12, 8))

# Temperature distribution plot
plt.subplot(2, 1, 1)
plt.plot(x, T_numerical, 'b-', label='Numerical')
plt.plot(x, T_analytical, 'r--', label='Analytical')
plt.xlabel('Distance (m)')
plt.ylabel('Temperature (°C)')
plt.title('Temperature Distribution along the Bar')
plt.legend()
plt.grid(True)

# Error plot
plt.subplot(2, 1, 2)
plt.plot(x, error_percentage, 'g-')
plt.xlabel('Distance (m)')
plt.ylabel('Error (%)')
plt.title('Percentage Error')
plt.grid(True)

plt.tight_layout()
plt.show()

print(f"Heat Transfer Rate: {heat_transfer_rate:.2f} W")