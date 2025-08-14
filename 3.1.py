import numpy as np
import matplotlib.pyplot as plt

def solve_cylinder_heat_transfer(L, r1, r2, k, bc1, bc2, n):
    # Parameters:
    # L: length of cylinder (m)
    # r1: inner radius (m)
    # r2: outer radius (m)
    # k: thermal conductivity (W/m°C)
    # bc1: temperature at inner radius (°C)
    # bc2: temperature at outer radius (°C)
    # n: number of points between boundaries
    
    m = (r2-r1)/(n+1)  # Step size
    r = np.array([r1 + i*m for i in range(1, n+1)])
    
    # Initialize coefficient matrix A
    A = np.zeros((n, n))
    
    # Fill diagonal elements
    for i in range(n):
        A[i,i] = -2*r[i]/m**2
        
    # Fill off-diagonal elements
    for i in range(n):
        if i == 0:
            A[i,i+1] = r[i]/m**2 + 1/(2*m)
        elif i == n-1:
            A[i,i-1] = r[i]/m**2 - 1/(2*m)
        else:
            A[i,i-1] = r[i]/m**2 - 1/(2*m)
            A[i,i+1] = r[i]/m**2 + 1/(2*m)
            
    # Initialize and fill boundary condition vector B
    B = np.zeros(n)
    B[0] = -(r[0]/m**2 - 1/(2*m))*bc1
    B[-1] = -(r[-1]/m**2 + 1/(2*m))*bc2
    
    # Solve system of equations
    T = np.linalg.solve(A, B)
    
    # Add boundary temperatures
    radius = np.concatenate(([r1], r, [r2]))
    temp = np.concatenate(([bc1], T, [bc2]))
    
    # Calculate heat transfer rate
    Q = k*2*np.pi*r1*L*(bc1-T[0])/(r[0]-r1)
    
    # Calculate analytical solution
    r_analytical = radius
    T_analytical = bc1 + (bc2-bc1)*np.log(r_analytical/r1)/np.log(r2/r1)
    Q_analytical = 2*np.pi*k*L*(bc1-bc2)/np.log(r2/r1)
    
    # Calculate error
    error = np.abs(temp - T_analytical)/T_analytical * 100
    
    return radius, temp, T_analytical, error, Q, Q_analytical

def plot_results(radius, temp, T_analytical, error):
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

# Example usage
if __name__ == "__main__":
    # Problem parameters
    L = 4  # Length of cylinder (m)
    r1 = 0.01  # Inner radius (m)
    r2 = 0.05  # Outer radius (m)
    k = 10  # Thermal conductivity (W/m°C)
    bc1 = 100  # Inner surface temperature (°C)
    bc2 = 40  # Outer surface temperature (°C)
    n = 400  # Number of points
    
    # Solve problem
    radius, temp, T_analytical, error, Q, Q_analytical = solve_cylinder_heat_transfer(
        L, r1, r2, k, bc1, bc2, n)
    
    # Plot results
    fig = plot_results(radius, temp, T_analytical, error)
    plt.show()
    
    print(f"Numerical heat transfer rate: {Q:.2f} W")
    print(f"Analytical heat transfer rate: {Q_analytical:.2f} W")
    print(f"Maximum error: {np.max(error):.6f}%")