# Darcy-Forchheimer Porosity Analyzer 

A Python-based **algebraic model** (single-point calculator) that aims to replicate OpenFOAM's `explicitPorositySource` tensor transformations. 

This tool is designed to help CFD engineers debug mass flow choking and unphysical pressure spikes caused by poorly tuned transverse resistance coefficients in complex aerodynamic models (e.g., automotive cooling packs, radiators, and pitched fans).

## 🛑 The Problem: The "Brick Wall" Bug
Many academic tutorials suggest using astronomically high coefficients (e.g., $10^7$) for transverse directions to enforce flow directionality. While this works in simple straight ducts, applying these "brick wall" coefficients to real-world models with sideways cross-flows (like an engine bay) causes the Navier-Stokes solver to artificially kill the kinetic energy, destroying mass conservation.

Industry standards (such as those in commercial solver documentation) recommend using a **100x to 1000x multiplier** on the Forchheimer (inertial) term to smoothly guide the flow without crashing the solver. This repository provides the tools to test and visualize those coefficients *before* running a computationally expensive simulation.

---

## 🛠️ The Tools

### 1. `porous_engine.py` (The Core Engine)
A faithful Python replication of the OpenFOAM C++ source code. It handles the local-to-global tensor transformations (`axesRotation`) and computes the exact momentum sink vector ($S_x, S_y, S_z$) for any given fluid velocity.

### 2. `fan_curve_analyzer.py` (Velocity Sweep)
Tests your coefficients across a range of bulk velocities (e.g., 0.1 to 1.0 m/s). 
* Isolates the intended streamwise resistance ($S_x$) from the artificial transverse resistance ($S_y, S_z$).
* Plots the total pressure drop magnitude $|S|$.
* Handles **negative Darcy coefficients** (often used to simulate active fans/pumps), clearly showing the zero-crossing where the zone stops pumping and starts acting as a brake.

### 3. `yaw_sweep_analyzer.py` (Angle Sweep)
Tests your fan's sensitivity to sideways cross-flow.
* Locks the total velocity magnitude (e.g., 5 m/s) and sweeps the incoming flow angle from 0° (perfectly aligned) to 90° (fully sideways).
* Reveals if your transverse Darcy/Forchheimer coefficients are set too high, allowing you to tune them to the industry-standard 100x multiplier before running your 3D mesh.

--- 
⚠️ Important Note: Sink Magnitude (Pa/m) vs. Actual Pressure Drop (Pa)It is critical to note that the output of these tools is the Momentum Sink Magnitude ($|S|$), which has units of Pascals per meter (Pa/m). It represents the resistance density.To calculate the actual total pressure drop ($\Delta P$) your fluid will experience in the 3D OpenFOAM simulation, you must multiply the tool's output by the physical thickness ($L$) of your meshed cellZone:$$ \Delta P = |S| \times L $$Example: If the analyzer predicts a total magnitude of 5000 Pa/m at a 20° cross-flow, and your modeled porous zone is 0.1 meters (10 cm) thick, your actual pressure drop across the zone will be 500 Pa.
--- 

**Running the Analyzers:**
Simply open either analyzer script, edit the `d_coeffs`, `f_coeffs`, and `pitch` variables in the Setup section to match your OpenFOAM `fvOptions` dictionary, and run:

```bash
python fan_curve_analyzer.py
# or
python yaw_sweep_analyzer.py
```

