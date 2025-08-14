import numpy as np
import matplotlib.pyplot as plt

def solve_cylinder_heat_transfer_convective(L, r1, r2, k, ta1, ha1, ta2, ha2, n):
    """
    Solve heat transfer in cylinder with convective boundary conditions
    """
    m = (r2-r1)/(n+1)  # Step size
    r = np.array([r1 + i*m for i in range(1, n+1)])
    
    # Initialize coefficient matrix A
    A = np.zeros((n, n))
    
    # Fill matrix elements
    for i in range(n):
        A[i,i] = -2*r[i]/m**2
        
        if i == 0:
            A[i,i+1] = (r[i]/m**2 + 1/(2*m)) + (r[i]/m**2 - 1/(2*m))*(1/(1+2*m*ha1/k))
        elif i == n-1:
            A[i,i-1] = (r[i]/m**2 - 1/(2*m)) + (r[i]/m**2 + 1/(2*m))*(1/(1+2*m*ha2/k))
        else:
            A[i,i-1] = r[i]/m**2 - 1/(2*m)
            A[i,i+1] = r[i]/m**2 + 1/(2*m)
    
    # Initialize and fill boundary condition vector B
    B = np.zeros(n)
    B[0] = -(r[0]/m**2 - 1/(2*m))*(1/(1+k/(2*m*ha1)))*ta1
    B[-1] = -(r[-1]/m**2 + 1/(2*m))*(1/(1+k/(2*m*ha2)))*ta2
    
    # Solve system of equations
    T = np.linalg.solve(A, B)
    
    # Calculate boundary temperatures
    T0 = T[1]/(1+2*m*ha1/k) + ta1/(1+k/(2*m*ha1))
    Tn = T[-2]/(1+2*m*ha2/k) + ta2/(1+k/(2*m*ha2))
    
    # Combine all temperatures and radii
    radius = np.concatenate(([r1], r, [r2]))
    temp = np.concatenate(([T0], T, [Tn]))
    
    # Calculate analytical solution
    r_analytical = radius
    T_analytical = ta1 + (ta2-ta1)*(np.log(r_analytical/r1)+k/(ha1*r1))/(k/(ha1*r1) + k/(ha2*r2) + np.log(r2/r1))
    
    # Calculate heat transfer rate
    Q = -k*2*np.pi*r[0]*L*(T[1]-T[0])/(r[1]-r[0])
    Q_analytical = 2*np.pi*k*L*(ta1-ta2)/(k/(ha1*r1) + k/(ha2*r2) + np.log(r2/r1))
    
    # Calculate error
    error = np.abs(temp - T_analytical)/T_analytical * 100
    
    return radius, temp, T_analytical, error, Q, Q_analytical

def plot_results(radius, temp, T_analytical, error):
    """Plot temperature distribution and error"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    ax1.plot(radius, temp, 'b-', label='Numerical')
    ax1.plot(radius, T_analytical, 'k--', label='Analytical')
    ax1.set_xlabel('Radius (m)')
    ax1.set_ylabel('Temperature (°C)')
    ax1.legend()
    ax1.grid(True)
    
    ax2.plot(radius, error, 'r-')
    ax2.set_xlabel('Radius (m)')
    ax2.set_ylabel('Absolute Error (%)')
    ax2.grid(True)
    
    plt.tight_layout()
    return fig

# Example usage
if __name__ == "__main__":
    # Problem parameters
    L = 4  # Length (m)
    r1 = 0.01  # Inner radius (m)
    r2 = 0.05  # Outer radius (m)
    k = 10  # Thermal conductivity (W/m°C)
    ta1 = 250  # Inner gas temperature (°C)
    ha1 = 8.17  # Inner convection coefficient (W/m²°C)
    ta2 = 35  # Outer gas temperature (°C)
    ha2 = 24  # Outer convection coefficient (W/m²°C)
    n = 1600  # Number of points
    
    # Solve problem
    radius, temp, T_analytical, error, Q, Q_analytical = solve_cylinder_heat_transfer_convective(
        L, r1, r2, k, ta1, ha1, ta2, ha2, n)
    
    # Plot results
    fig = plot_results(radius, temp, T_analytical, error)
    plt.show()
    
    print(f"Numerical heat transfer rate: {Q:.2f} W")
    print(f"Analytical heat transfer rate: {Q_analytical:.2f} W")
    print(f"Maximum error: {np.max(error):.6f}%")