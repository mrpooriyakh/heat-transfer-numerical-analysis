import numpy as np
import matplotlib.pyplot as plt

def solve_sphere_heat_transfer(r1, r2, k, bc1, bc2, n):
    """
    Solve heat transfer in sphere with temperature boundary conditions
    """
    # Step size
    m = (r2-r1)/(n+1)
    
    # Generate radius points 
    r = np.array([r1 + i*m for i in range(1, n+1)])
    
    # Initialize coefficient matrix A
    A = np.zeros((n, n))
    
    # Fill matrix elements
    for i in range(n):
        # Diagonal elements
        A[i,i] = -2*r[i]**2/m**2
        
        if i == 0:
            # First row
            A[i,i+1] = r[i]**2/m**2 + r[i]/(m)
        elif i == n-1:
            # Last row
            A[i,i-1] = r[i]**2/m**2 - r[i]/(m)
        else:
            # Middle rows
            A[i,i-1] = r[i]**2/m**2 - r[i]/(m)
            A[i,i+1] = r[i]**2/m**2 + r[i]/(m)
            
    # Initialize and fill boundary condition vector B
    B = np.zeros(n)
    B[0] = -(r[0]**2/m**2 - r[0]/(m))*bc1
    B[-1] = -(r[-1]**2/m**2 + r[-1]/(m))*bc2
    
    # Solve system of equations
    T = np.linalg.solve(A, B)
    
    # Combine with boundary points
    radius = np.concatenate(([r1], r, [r2]))
    temp = np.concatenate(([bc1], T, [bc2]))
    
    # Analytical solution
    T_analytical = r1*r2*(bc1-bc2)/(radius*(r2-r1)) + (r2*bc2-r1*bc1)/(r2-r1)
    Q_analytical = (4*np.pi*k*r2*r1*(bc1-bc2)/(r2-r1))/1e3  # Convert to kW
    
    # Calculate numerical heat transfer rates at each point
    Q = np.zeros(n-1)
    for i in range(1, n-1):
        Q[i-1] = (k*4*np.pi*r[i]**2*(T[i-1]-T[i+1])/(2*m))/1e3
    
    # Calculate error
    error = np.abs(temp - T_analytical)/T_analytical * 100
    
    return radius, temp, T_analytical, error, Q[0], Q_analytical

def plot_results(radius, temp, T_analytical, error):
    """Plot temperature distribution and error"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Temperature plot
    ax1.plot(radius, temp, 'b-', label='Numerical')
    ax1.plot(radius, T_analytical, 'k--', label='Analytical')
    ax1.set_xlabel('Radius (m)')
    ax1.set_ylabel('Temperature (°C)')
    ax1.legend()
    ax1.grid(True)
    
    # Error plot
    ax2.plot(radius, error, 'r-')
    ax2.set_xlabel('Radius (m)')
    ax2.set_ylabel('Absolute Error (%)')
    ax2.grid(True)
    
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    # Problem parameters from the MATLAB example
    r1 = 0.08  # Inner radius (m)
    r2 = 0.16  # Outer radius (m)
    k = 45     # Thermal conductivity (W/m°C)
    bc1 = 200  # Inner surface temperature (°C)
    bc2 = 80   # Outer surface temperature (°C)
    n = 1600   # Number of points
    
    # Solve problem
    radius, temp, T_analytical, error, Q, Q_analytical = solve_sphere_heat_transfer(
        r1, r2, k, bc1, bc2, n)
    
    # Plot results
    fig = plot_results(radius, temp, T_analytical, error)
    plt.show()
    
    print(f"Numerical heat transfer rate: {Q:.3f} kW")
    print(f"Analytical heat transfer rate: {Q_analytical:.3f} kW")
    print(f"Maximum error: {np.max(error):.6f}%")