import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tickerE
from porous_engine import DarcyForchheimerZone

##################################################################
# 1. Porous Zone Setup 
##################################################################
d_coeffs = [-33262, 33262, 33262] 
f_coeffs = [14.4844, 1450, 1450]

pitch_rad = np.deg2rad(4)
e1_vector = [np.cos(pitch_rad), 0, np.sin(pitch_rad)] 
e2_vector = [0, 1, 0]

fan_zone = DarcyForchheimerZone(d=d_coeffs, f=f_coeffs, e1=e1_vector, e2=e2_vector)

##################################################################
# 2. Yaw Angle Sweep (Misalignment Testing at 5 m/s)
##################################################################
v_mag = 5.0 
z_noise = 0.5 

angles_deg = np.linspace(0, 90, 100)
angles_rad = np.deg2rad(angles_deg)

pressure_drops_global_x = []
pressure_drops_global_y = []
pressure_drops_global_z = []
pressure_drops_total = []

for angle in angles_rad:
    Ux = v_mag * np.cos(angle)
    Uy = v_mag * np.sin(angle)
    Uz = z_noise 
    
    U_test = [Ux, Uy, Uz] 
    
    S_glob, S_loc, U_loc = fan_zone.evaluate_sink(U_test)
    
    pressure_drops_global_x.append(S_glob[0])
    pressure_drops_global_y.append(S_glob[1])
    pressure_drops_global_z.append(S_glob[2])
    pressure_drops_total.append(np.linalg.norm(S_glob))

##################################################################
# 3. Plot the Results
##################################################################
fig, ax = plt.subplots(figsize=(10, 6))

color_x = '#1f77b4'    # Muted Blue
color_y = '#ff7f0e'    # Safety Orange
color_z = '#d62728'    # Brick Red
color_total = '#2ca02c' # Forest Green 

# Main Plotting
ax.plot(angles_deg, pressure_drops_global_x, label='$S_x$ Component', color=color_x, linewidth=2)
ax.plot(angles_deg, pressure_drops_global_y, label='$S_y$ Component', color=color_y, linestyle='-.', linewidth=2)
ax.plot(angles_deg, pressure_drops_global_z, label='$S_z$ Component', color=color_z, linestyle='--', linewidth=2)
ax.plot(angles_deg, pressure_drops_total, label='Total Magnitude $|S|$', color=color_total, linewidth=2.5)

# 0 Axis
ax.axhline(0, color='black', linewidth=1.2, linestyle='-')

ax.set_xlabel('Incoming Flow Angle (0° = Perfect, 90° = Fully Sideways)', fontweight='bold')
ax.set_ylabel('Momentum Sink/Source (Pa/m)', fontweight='bold')
ax.set_title('Fan Pressure Drop Sensitivity vs. Sideways Flow (v= 5 m/s)', pad=15)

ax.legend(loc='upper left')
ax.grid(True, linestyle='--', alpha=0.7)
ax.axhline(0, color='black', linewidth=1.5, linestyle='-')

plt.tight_layout()
plt.show()