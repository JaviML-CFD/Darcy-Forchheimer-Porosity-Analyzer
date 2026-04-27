import numpy as np

class DarcyForchheimerZone:
    """
    Faithful Python replica of OpenFOAM's explicitPorositySource 
    using DarcyForchheimerCoeffs and axesRotation.
    """
    def __init__(self, d, f, e1, e2, nu=1.5e-5):
        self.nu = nu
        self.d = np.array(d, dtype=float)
        self.f = np.array(f, dtype=float)
        
        # OpenFOAM axesRotation logic
        e1 = np.array(e1, dtype=float)
        e2 = np.array(e2, dtype=float)
        
        self.e1 = e1 / np.linalg.norm(e1)
        e2_ortho = e2 - np.dot(e2, self.e1) * self.e1
        self.e2 = e2_ortho / np.linalg.norm(e2_ortho)
        self.e3 = np.cross(self.e1, self.e2)
        
        self.E = np.vstack((self.e1, self.e2, self.e3))
        self.E_T = self.E.T

    def evaluate_sink(self, U_global):
        U_global = np.array(U_global, dtype=float)
        magU = np.linalg.norm(U_global)
        
        if magU == 0:
            return np.zeros(3), np.zeros(3), np.zeros(3)

        U_local = np.dot(self.E, U_global)
        viscous_term = self.nu * self.d
        inertial_term = 0.5 * magU * self.f
        
        S_local = -(viscous_term + inertial_term) * U_local
        S_global = np.dot(self.E_T, S_local)
        
        return S_global, S_local, U_local