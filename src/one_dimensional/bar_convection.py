import numpy as np
import matplotlib.pyplot as plt

def solve_heat_transfer_convective(x0, xL, k, area, ha1, ta1, ha2, ta2, n):
    m = (xL - x0)/(n + 1)
    x = np.linspace(x0, xL, n + 2)
    
    # Initialize matrices
    A = np.zeros((n, n))
    B = np.zeros(n)
    
    # Main diagonal elements
    for i in range(n):
        A[i,i] = -2/m**2
        
    # First row
    A[0,1] = 1/m**2 + 1/m**2 * (1/(1+2*m*ha1/k))
    B[0] = -(1/m**2 * (1/(1+k/(2*m*ha1))))*ta1
    
    # Middle rows
    for i in range(1, n-1):
        A[i,i-1] = 1/m**2
        A[i,i+1] = 1/m**2
        
    # Last row
    A[n-1,n-2] = 1/m**2 + 1/m**2 * (1/(1+2*m*ha2/k))
    B[n-1] = -(1/m**2 * (1/(1+k/(2*m*ha2))))*ta2
    
    # Solve system
    T_inner = np.linalg.solve(A, B)
    
    # Calculate temperatures including boundaries
    T = np.zeros(n + 2)
    T[0] = (T_inner[1]/(1+2*m*ha1/k))+(ta1/(1+(k/(2*m*ha1))))
    T[1:-1] = T_inner
    T[-1] = (T_inner[-2]/(1+2*m*ha2/k))+(ta2/(1+(k/(2*m*ha2))))
    
    # Analytical solution
    T_analytical = np.array([ta1 - (ta1-ta2)*(k/ha1+xi)/(k/ha1+k/ha2+xL) for xi in x])
    
    return x, T, T_analytical

# Test case
x0, xL = 0, 0.4
k = 2.3
area = 20
ha1, ta1 = 30, 80
ha2, ta2 = 24, 15
n = 20

x, T, T_analytical = solve_heat_transfer_convective(x0, xL, k, area, ha1, ta1, ha2, ta2, n)

plt.figure(figsize=(10, 8))
plt.subplot(2, 1, 1)
plt.plot(x, T, 'b-', label='Numerical')
plt.plot(x, T_analytical, 'r--', label='Analytical')
plt.xlabel('Distance (m)')
plt.ylabel('Temperature (°C)')
plt.legend()
plt.grid(True)

plt.subplot(2, 1, 2)
error = np.abs(T - T_analytical)*100/T_analytical
plt.plot(x, error, 'k-')
plt.xlabel('Distance (m)')
plt.ylabel('Absolute Error (%)')
plt.grid(True)
plt.tight_layout()
plt.show()