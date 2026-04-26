# Darcy-Forchheimer Porosity Analyzer for OpenFOAM

A Python-based standalone analytical solver that replicates OpenFOAM's `explicitPorositySource` tensor mathematics. 

This tool is designed to try to help CFD engineers debug, visualize, and tune Darcy-Forchheimer porous media coefficients (e.g., for modeling cooling fans, radiators, or flow straighteners) without needing to run computationally expensive 3D simulations.

## The Problem: The "Brick Wall" Pressure Drop
When setting up directional porous zones in OpenFOAM, it is standard practice to apply artificially high resistance coefficients in the transverse directions ($y$ and $z$) to prevent cross-flow. 

However, if the bulk flow is slightly misaligned with the porous zone's local coordinate system (e.g., due to a physical pitch angle or numerical noise), OpenFOAM takes that slight transverse velocity and multiplies it by the massive transverse resistance coefficient. 

Because the viscous Darcy term is linear, it does not drop off at near-zero velocities:
$$S = - \left( \nu D + \frac{1}{2} |U| F \right) U$$

This results in a massive, unphysical pressure spike across the porous zone, even as the bulk velocity approaches zero.

## The Solution
This script replicates the C++ source code architecture of OpenFOAM's `axesRotation` and local-to-global tensor transformations. By providing your exact $D$ and $F$ vectors, along with your geometric flow vectors ($e_1$ and $e_2$), you can plot the exact momentum sink curve and verify that your coordinate rotation safely passes the bulk flow while blocking transverse noise.

## Features
* **Faithful OpenFOAM Replication:** Translates OpenFOAM's local coordinate transformation matrix ($E$) into `numpy` array operations.
* **Component-wise Sink Evaluation:** Calculates both Viscous (linear) and Inertial (quadratic) pressure drops.
* **Misalignment Testing:** Sweep through velocity magnitudes and inject artificial transverse velocity to see exactly when and where your porous zone will cause a pressure spike.
