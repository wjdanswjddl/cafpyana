import numpy as np

# CONSTANTS

# RECOMBINATION
MODA = 0.906
B90 = 0.203
R = 1.25
Wion = 1e3 / 4.237e7
def ellipsoid_beta(phi, B90=B90, R=R):
    return B90 / np.sqrt(np.sin(phi)**2 + np.cos(phi)**2/R**2)

# TODO: check these

# SBND detector parameters
LAr_density_gmL_sbnd = 1.38434
Efield_sbnd = 0.5                           

# ICARUS detector parameters
LAr_density_gmL_IC = 1.390
Efield_IC = 0.4938

def recombination_cor(dQdx, phi, E=0.5, rho=1.39, A=MODA, B90=B90, R=R):
    alpha = A
    beta = ellipsoid_beta(phi, B90, R) / (rho * E)

    dEdx = (np.exp(dQdx*Wion*beta)- alpha) / beta

    return dEdx

def recombination(dEdx, phi, E=0.5, rho=1.39, A=MODA, B90=B90, R=R):
    alpha = A
    beta = ellipsoid_beta(phi, B90, R) / (rho * E)

    dQdx = np.log(alpha + dEdx*beta) / (Wion * beta)
    return dQdx
