# Documentation on how to use derivedParameters.py can be found in the
# NFDI4Ing knowledge base:
# https://nfdi4ing.pages.rwth-aachen.de/knowledge-base/how-tos/how_to_set_dependent_parameters_in_a_parameter_study_using_pyfoam_with_derivedparameters_py/

def computeDeltaT(values):
    import math
    import sys

    # Give explicitly prescribed time step sizes precedence over
    # computed values
    if "DELTA_T" in values:
        return values["DELTA_T"]

    
    domainLength = values["DOMAIN_LENGTH"]
    resolution = values["RESOLUTION"]
    h = domainLength/resolution

    # rho_ambient might be expressed in terms of rho_droplet
    rho_ambient = 0.0
    rho_droplet = 0.0
    sigma = 0.0
    if values["FLUID_PAIRING"] == "water-air":
        rho_ambient = 1.0
        rho_droplet = 1000.0
        sigma = 72.75e-3
    elif values["FLUID_PAIRING"] == "oil_novec7500-water": 
        rho_ambient = 1000.0
        rho_droplet = 1614.0
        sigma = 49.5e-3
    elif values["FLUID_PAIRING"] == "gearoil-air":
        rho_ambient = 1.0
        rho_droplet = 888.0
        sigma = 30.0e-3
    else:
        sys.exit("Error: unknown FLUID PAIRING:" + values["FLUID_PAIRING"])

    # Coefficient to scale the capillary time step size
    sf = values["SCALE_CAPILLARY_DELTA_T"]

    # The computation follows eq. (18) in Popinet 2009 and
    # eq. in (43) Denner & van Wachem 2015, respectively.
    capillary_delta_t = 1e15
    if sigma > 0.0:
        capillary_delta_t = sf*math.sqrt((rho_droplet + rho_ambient)*math.pow(h,3.0)/(2.0*math.pi*sigma))

    return capillary_delta_t

DELTA_T = computeDeltaT(locals())

# This dummy is required to avoid a PyFoam error if DELTA_T is given explicitly in the
# parameter file
dummy = 0