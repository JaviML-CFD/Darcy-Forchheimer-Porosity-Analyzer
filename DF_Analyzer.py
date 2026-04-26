import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

class DarcyForchheimerZone:
    """
    Faithful Python replica of OpenFOAM's explicitPorositySource 
    using DarcyForchheimerCoeffs and axesRotation.
    """
    def __init__(self, d, f, e1, e2, nu=1.5e-5):
        # Kinematic viscosity (default: air at standard conditions)
        self.nu = nu
        
        # Darcy (D) and Forchheimer (F) local diagonal tensors
        self.d = np.array(d, dtype=float)
        self.f = np.array(f, dtype=float)
        
        # --- OpenFOAM axesRotation logic ---
        e1 = np.array(e1, dtype=float)
        e2 = np.array(e2, dtype=float)
        
        # 1. Normalize e1
        self.e1 = e1 / np.linalg.norm(e1)
        
        # 2. Orthogonalize e2 with respect to e1 (Gram-Schmidt), then normalize
        e2_ortho = e2 - np.dot(e2, self.e1) * self.e1
        self.e2 = e2_ortho / np.linalg.norm(e2_ortho)
        
        # 3. e3 is the cross product of e1 and e2
        self.e3 = np.cross(self.e1, self.e2)
        
        # Transformation Tensor (E): Converts Global -> Local
        self.E = np.vstack((self.e1, self.e2, self.e3))
        
        # Transpose Tensor (E^T): Converts Local -> Global
        self.E_T = self.E.T

    def evaluate_sink(self, U_global):
        """
        Calculates the momentum sink for a given global velocity vector.
        Returns: S_global, S_local, U_local
        """
        U_global = np.array(U_global, dtype=float)
        magU = np.linalg.norm(U_global)
        
        if magU == 0:
            return np.zeros(3), np.zeros(3), np.zeros(3)

        # Transform global velocity to local porous coordinates
        U_local = np.dot(self.E, U_global)
        
        # Apply Darcy-Forchheimer equation locally (component-wise)
        # S = -(nu * D + 0.5 * |U| * F) * U_local
        viscous_term = self.nu * self.d
        inertial_term = 0.5 * magU * self.f
        
        S_local = -(viscous_term + inertial_term) * U_local
        
        # Transform the resulting sink back to global coordinates
        S_global = np.dot(self.E_T, S_local)
        
        return S_global, S_local, U_local
    
##################################################################
# 1. Porous Zone Setup
##################################################################

d_coeffs = [-33262, 100000, 100000] 
f_coeffs = [14.4844, 100000, 100000]

# The fixed rotated axis (e.g., pitched slightly up in Z)
e1_vector = [0.9848, 0, 0.1736] 
e2_vector = [0, 1, 0]

# Initialize the zone
fan_zone = DarcyForchheimerZone(d=d_coeffs, f=f_coeffs, e1=e1_vector, e2=e2_vector)

##################################################################
# 2. Velociity Sweep
##################################################################

velocities = np.linspace(0.1, 0.5, 50)
pressure_drops_global_x = []
pressure_drops_global_z = []
pressure_drops_total = []

for v in velocities:

    U_test = [0, v, 0] 
    
    S_glob, S_loc, U_loc = fan_zone.evaluate_sink(U_test)
    
    pressure_drops_global_x.append(S_glob[0])
    pressure_drops_global_z.append(S_glob[2])
    pressure_drops_total.append(np.linalg.norm(S_glob))

##################################################################
# 3. Plot the results
##################################################################

fig, ax = plt.subplots(figsize=(10, 6)) # Explicitly create figure and axis

ax.plot(velocities, pressure_drops_global_x, label='S_x component')
ax.plot(velocities, pressure_drops_global_z, label='S_z component', linestyle='--')
ax.plot(velocities, pressure_drops_total, label='Total Magnitude |S|', color='black', linewidth=2)

ax.yaxis.set_major_locator(ticker.MaxNLocator(20)) 
ax.xaxis.set_major_locator(ticker.MaxNLocator(15))

ax.set_xlabel('Bulk Velocity (m/s)')
ax.set_ylabel('Momentum Sink/Source (Pa/m)')
ax.set_title('OpenFOAM Darcy-Forchheimer Replication')
ax.legend()
ax.grid(True)

ax.axhline(0, color='red', linewidth=1.5, linestyle=':')

plt.show()