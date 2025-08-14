import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def solve_2d_heat_transfer(n, x0, xf, y0, yf, bc1, bc2, bc3, bc4):
    """
    Solve 2D heat transfer problem using finite difference method
    """
    # Calculate step sizes
    h = (xf - x0)/(n+1)  # x direction
    k = (yf - y0)/(n+1)  # y direction
    
    # Initialize matrices
    a = np.zeros((n*n, n*n))
    b = np.zeros(n*n)

    # Fill coefficient matrix
    for i in range(n*n):
        # Main diagonal elements
        a[i,i] = -2/h**2 - 2/k**2
        
        # Connect to point above
        if i >= n:
            a[i,i-n] = 1/h**2
            
        # Connect to point below
        if i < n*(n-1):
            a[i,i+n] = 1/h**2
            
        # Connect to point on left (if not at left edge)
        if i % n != 0:
            a[i,i-1] = 1/k**2
            
        # Connect to point on right (if not at right edge)
        if (i+1) % n != 0:
            a[i,i+1] = 1/k**2

    # Fill boundary condition vector
    for i in range(n*n):
        row = i // n
        col = i % n
        
        # Left boundary
        if col == 0:
            b[i] -= bc2/k**2
            
        # Right boundary
        if col == n-1:
            b[i] -= bc3/k**2
            
        # Bottom boundary
        if row == n-1:
            b[i] -= bc1/h**2
            
        # Top boundary
        if row == 0:
            b[i] -= bc4/h**2

    # Solve system
    t = np.linalg.solve(a, b)

    # Reshape solution into grid form
    T = np.zeros((n+2, n+2))
    
    # Fill interior points
    for i in range(n):
        for j in range(n):
            T[i+1,j+1] = t[i*n + j]
    
    # Fill boundary points
    T[0,:] = bc4  # Top boundary
    T[-1,:] = bc1  # Bottom boundary
    T[:,0] = bc2  # Left boundary
    T[:,-1] = bc3  # Right boundary
    
    return T

def plot_results_enhanced(T):
    """
    Plot the temperature distribution with multiple views
    """
    n = T.shape[0]
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, y)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(15, 10))
    
    # 3D surface plot
    ax1 = fig.add_subplot(221, projection='3d')
    surf = ax1.plot_surface(X, Y, T, cmap='hot')
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_zlabel('Temperature (°C)')
    ax1.set_title('3D Temperature Distribution')
    plt.colorbar(surf, ax=ax1, label='Temperature (°C)')
    
    # X-Y view (top view)
    ax2 = fig.add_subplot(222)
    c2 = ax2.contourf(X, Y, T, levels=20, cmap='hot')
    plt.colorbar(c2, ax=ax2, label='Temperature (°C)')
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('y (m)')
    ax2.set_title('Top View (X-Y plane)')
    
    # X-Z view (side view)
    ax3 = fig.add_subplot(223)
    mid_y_idx = T.shape[0] // 2
    ax3.plot(x, T[mid_y_idx,:], 'r-', label=f'y = {y[mid_y_idx]:.2f}m')
    ax3.set_xlabel('x (m)')
    ax3.set_ylabel('Temperature (°C)')
    ax3.set_title('Side View (X-Z plane)')
    ax3.grid(True)
    ax3.legend()
    
    # Y-Z view (front view)
    ax4 = fig.add_subplot(224)
    mid_x_idx = T.shape[1] // 2
    ax4.plot(y, T[:,mid_x_idx], 'b-', label=f'x = {x[mid_x_idx]:.2f}m')
    ax4.set_xlabel('y (m)')
    ax4.set_ylabel('Temperature (°C)')
    ax4.set_title('Front View (Y-Z plane)')
    ax4.grid(True)
    ax4.legend()
    
    plt.tight_layout()
    plt.show()

# Example usage
if __name__ == "__main__":
    # Problem parameters
    n = 3  # number of internal points in each direction
    x0, xf = 0, 1  # x domain
    y0, yf = 0, 1  # y domain
    
    # Boundary conditions
    bc1 = 35  # temperature at T(x,0)
    bc2 = 60  # temperature at T(0,y)
    bc3 = 100  # temperature at T(L,y)
    bc4 = 70  # temperature at T(x,B)
    
    # Solve problem
    T = solve_2d_heat_transfer(n, x0, xf, y0, yf, bc1, bc2, bc3, bc4)
    
    # Display numerical results
    print("\nTemperature Distribution:")
    print(T)
    
    # Plot results with enhanced visualization
    plot_results_enhanced(T)