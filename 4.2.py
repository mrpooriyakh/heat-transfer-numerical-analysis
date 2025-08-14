import numpy as np
import matplotlib.pyplot as plt

def solve_sphere_heat_transfer_convective(r1, r2, k, ta1, hg1, ta2, hg2, n):
    """
    Solve heat transfer in sphere with convective boundary conditions
    """
    # Step size
    c = (r2-r1)/(n+1)
    
    # Generate radius points
    r = np.array([r1 + i*c for i in range(1, n+1)])
    
    # Initialize matrices
    A = np.zeros((n, n))
    B = np.zeros(n)
    
    for i in range(n):
        # Main diagonal terms
        A[i,i] = -2*r[i]**2/c**2
        
        if i == 0:
            A[i,i+1] = (r[i]**2/c**2 + r[i]/c) + \
                       (r[i]**2/c**2 - r[i]/c)/(1+2*c*hg1/k)
            B[i] = -(r[i]**2/c**2 - r[i]/c)/(1+k/(2*c*hg1))*ta1
        
        elif i == n-1:
            A[i,i-1] = (r[i]**2/c**2 - r[i]/c) + \
                       (r[i]**2/c**2 + r[i]/c)/(1+2*c*hg2/k)
            B[i] = -(r[i]**2/c**2 + r[i]/c)/(1+k/(2*c*hg2))*ta2
        
        else:
            A[i,i-1] = r[i]**2/c**2 - r[i]/c
            A[i,i+1] = r[i]**2/c**2 + r[i]/c
    
    # Solve system
    T = np.linalg.solve(A, B)
    
    # Calculate boundary temperatures
    T0 = T[1]/(1+2*c*hg1/k) + ta1/(1+k/(2*c*hg1))
    Tn = T[-2]/(1+2*c*hg2/k) + ta2/(1+k/(2*c*hg2))
    
    # Combine all points
    radius = np.concatenate(([r1], r, [r2]))
    temp = np.concatenate(([T0], T, [Tn]))
    
    # Analytical solution
    T_analytical = ((ta1-ta2)/((1/r2 - 1/r1) - k*(1/(hg1*r1**2) + 1/(hg2*r2**2))))*\
                  (-1/radius + (1/2)*(k*(1/(hg1*r1**2)-1/(hg2*r2**2))+(1/r1+1/r2))) + \
                  (ta1+ta2)/2
    
    # Calculate heat transfer rates
    Q = -k*4*np.pi*r[0]**2*(T[1]-T[0])/(c)
    Q_analytical = -(4*np.pi*k*(ta1-ta2)/(1/r2 - 1/r1 - k/(hg1*r1**2) - k/(hg2*r2**2)))
    
    # Calculate error
    error = np.abs(temp - T_analytical)/T_analytical * 100
    
    return radius, temp, T_analytical, error, Q/1000, Q_analytical/1000

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
    # Problem parameters
    r1 = 80e-2  # Inner radius (m)
    r2 = 90e-2  # Outer radius (m)
    k = 45      # Thermal conductivity (W/m°C)
    ta1 = 250   # Inner gas temperature (°C)
    ta2 = 25    # Outer gas temperature (°C)
    hg1 = 9     # Inner convection coefficient (W/m²°C)
    hg2 = 24    # Outer convection coefficient (W/m²°C)
    n = 1600    # Number of points
    
    # Solve problem
    radius, temp, T_analytical, error, Q, Q_analytical = solve_sphere_heat_transfer_convective(
        r1, r2, k, ta1, hg1, ta2, hg2, n)
    
    # Plot results
    fig = plot_results(radius, temp, T_analytical, error)
    plt.show()
    
    print(f"Numerical heat transfer rate: {Q:.3f} kW")
    print(f"Analytical heat transfer rate: {Q_analytical:.3f} kW")
    print(f"Maximum error: {np.max(error):.6f}%")