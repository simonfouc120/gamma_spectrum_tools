# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 09:34:15 2023

@author: sf270338
"""

import math 
import matplotlib.pyplot as plt
import numpy as np

detectors_Z_avg = {
    "gaseous_detectors": {
        "Argon": 18,
        "Argon_CO2_mix": 17.4
    },
    "scintillators": {
        "NaI_Tl": 32,
        "CsI_Tl": 54,
        "BGO": 27.16,
        "LaBr3": 40.5
    },
    "semiconductors": {
        "HPGe": 32,
        "CdTe": 50,
        "CdZnTe": 45.5
    }
}

print(detectors_Z_avg)



def compton_diff(E0) : 
    """
    This function calculates the energy of the photon and the electron after a Compton diffusion.
    It also plots the energy of the photon and the electron as a function of the angle of deviation.

    Parameters:
    E0 (float): Energy of the incident photon in keV

    Returns:
    E_retro (float): Energy of the photon after retro-diffusion
    E_fc (float): Energy of the electron after retro-diffusion
    """
    
    E_photon = np.array([])
    E_electron = np.array([])
    for teta in range(360):
        En = float(E0 / (1 + (E0 / 511) * (1 - np.cos((teta * 2 * np.pi) / 360))))
        E_photon = np.append(E_photon, En)
        E_electron = np.append(E_electron, float(E0 - En))
        if teta == 180:
            E_retro = En
            E_fc = np.float16(E0 - En)

    plt.figure("Compton Diffusion")
    plt.plot(E_photon, label="Energy of the Compton photon")
    plt.plot(E_electron, label="Energy of the Compton electron")
    plt.title("Compton Diffusion for incident Energy of " + str(E0) + " keV")
    plt.ylabel("Energy [keV]")
    plt.xlabel("Angle of deviation [°]")
    plt.legend()
    plt.show()
    return E_retro, E_fc
    
# On définit une fonction qui calcule la probabilité qu'un photon atteigne la surface du détecteur
import math

def probabilite_photon_atteint(S, d):
    """
    Calcule la probabilité qu'un photon atteigne une surface S à une distance d
    d'une source isotrope émettant des photons dans toutes les directions.

    Paramètres :
    S (float) : surface du détecteur (en mètres carrés)
    d (float) : distance entre la source et le détecteur (en mètres)

    Retourne :
    float : probabilité que le photon atteigne le détecteur
    """
    P = S / (4 * np.pi * d**2)
    return P




def estimated_proportion(energy, Z):
    """
    Approximate estimation of the proportions of photoelectric effect, Compton scattering,
    and pair production based on photon energy and atomic number Z of the material.

    Parameters:
    energy (float): photon energy in MeV
    Z (int): atomic number of the material

    Returns:
    dict: approximate proportions of the three effects
    """
    # Estimating proportions based on energy range
    if energy < 0.1:  # < 100 keV
        photo_prop = 0.8 + 0.2 * (Z / 100)  # Increases with Z
        compton_prop = 0.2
        pair_prop = 0.0
    elif energy < 1.022:  # Between 100 keV and 1 MeV
        photo_prop = 0.2 * (Z / 100)
        compton_prop = 0.8 - 0.2 * (Z / 100)
        pair_prop = 0.0
    else:  # E > 1.022 MeV (pair production possible)
        photo_prop = 0.1 * (Z / 100)
        compton_prop = 0.5 - 0.1 * (Z / 100)
        pair_prop = 0.4 + 0.6 * (Z / 100)  # Increases strongly with Z

    # Normalize to ensure the sum is 1
    total = photo_prop + compton_prop + pair_prop
    proportions = {
        'photoelectric': photo_prop / total,
        'compton': compton_prop / total,
        'pair_production': pair_prop / total
    }
    return proportions