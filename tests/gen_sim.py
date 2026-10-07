from petkinetics.simulation import simulation as sim_func
import numpy as np

"""
Python script for testing tac_simulation functions

"""

# ── File paths ────────────────────────────────────────────────────────────
framebased_bif      = r'path to frame based input function'
interpolated_bif    = r'path to interpolated input function'

base = r'path to output folder'

if __name__ == '__main__':

    organ     = 'liver'
    
    methods   = dict(
        Patlak = ['Patlak_tZero', 'Patlak_tStar_linear'],
        Logan  = ['dLogan_tStar', 'dLogan_tZero'])

    interval   = 20
    datapoints = 8
    noise_values = [0.01]
    tstart_values = [3600]

    print(f'\nRunning {organ} simulations...')
    Tissue_values = sim_func.define_tissue_values(organ)

    K1_values = np.linspace(0.1, 2.0,  interval)  # mL·cm⁻³·min⁻¹   | Delivery rate constant
    k2_values = np.linspace(0.1, 2.0, interval)   # min⁻¹           | Efflux rate constant
    k3_values = np.linspace(0.005, 0.1, interval) # min⁻¹           | Phosphorylation rate
    k4_values = np.linspace(0.0001, 0.1, interval) # min⁻¹          | Dephosphorylation rate (irreversible: k4=0)

    Vb_values = [0.]

    for vb in Vb_values:
        sim_func.run(
            K1_values, k2_values, k3_values, k4_values, 
            vb, 
            methods,
            datapoints, organ,
            noise_values, tstart_values, 
            framebased_bif, interpolated_bif,
            base, genFig=False)
        