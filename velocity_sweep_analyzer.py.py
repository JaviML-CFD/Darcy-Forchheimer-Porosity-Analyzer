import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from porous_engine import DarcyForchheimerZone
    
##################################################################
# 1. Porous Zone Setup
##################################################################
d_coeffs = [-33262, 33262, 33262] 
f_coeffs = [14.4844, 1450, 1450]

# The fixed rotated axis (e.g., pitched slightly up in Z)
e1_vector = [np.cos(np.deg2rad(4)), 0, np.sin(np.deg2rad(4))] 
e2_vector = [0, 1, 0]

# Initialize the zone
fan_zone = DarcyForchheimerZone(d=d_coeffs, f=f_coeffs, e1=e1_vector, e2=e2_vector)

##################################################################
# 2. Velociity Sweep
##################################################################
velocities = np.linspace(0, 6, 50)
pressure_drops_global_x = []
pressure_drops_global_y = []
pressure_drops_global_z = []
pressure_drops_total = []

for v in velocities:

    U_test = [0, v, 0] 
    
    S_glob, S_loc, U_loc = fan_zone.evaluate_sink(U_test)
    
    pressure_drops_global_x.append(S_glob[0])
    pressure_drops_global_y.append(S_glob[1])
    pressure_drops_global_z.append(S_glob[2])
    pressure_drops_total.append(np.linalg.norm(S_glob))

##################################################################
# 3. Plot the results
##################################################################
fig, ax = plt.subplots(figsize=(10, 6)) # Explicitly create figure and axis

color_x = '#1f77b4'    # Muted Blue
color_y = '#ff7f0e'    # Safety Orange
color_z = '#d62728'    # Brick Red
color_total = '#2ca02c' # Forest Green 

# Main Plotting
ax.plot(velocities, pressure_drops_global_x, label='$S_x$ Component', color=color_x, linewidth=2)
ax.plot(velocities, pressure_drops_global_y, label='$S_y$ Component', color=color_y, linestyle='-.', linewidth=2)
ax.plot(velocities, pressure_drops_global_z, label='$S_z$ Component', color=color_z, linestyle='--', linewidth=2)
ax.plot(velocities, pressure_drops_total, label='Total Magnitude $|S|$', color=color_total, linewidth=2.5)

# 0 Axis
ax.axhline(0, color='black', linewidth=1.2, linestyle='-')

ax.yaxis.set_major_locator(ticker.MaxNLocator(20)) 
ax.xaxis.set_major_locator(ticker.MaxNLocator(15))

ax.set_xlabel('Bulk Velocity (m/s)', fontsize=11, fontweight='bold')
ax.set_ylabel('Momentum Sink/Source (Pa/m)', fontsize=11, fontweight='bold')
ax.set_title('OpenFOAM Darcy-Forchheimer Replication', fontsize=14, pad=15)

ax.legend(loc='upper left', framealpha=0.9, edgecolor='black')
ax.grid(True, linestyle=':', alpha=0.7)



plt.tight_layout()
plt.show()