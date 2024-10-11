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
        "CdZnTe": 45.5,
        "Si": 14        
    }
}


EMCE = 511  # keV


def calculate_energy(E0, teta):
    """
    Calculate the energy of the photon after Compton scattering.

    Parameters:
    E0 (float): Energy of the incident photon in keV
    teta (float): Angle of deviation in degrees

    Returns:
    float: Energy of the photon after scattering
    """
    return float(E0 / (1 + (E0 / 511) * (1 - np.cos(np.deg2rad(teta)))))

def compton_diff(E0, E_photon = np.array([]), E_electron = np.array([])) :
    """
    This function calculates the energy of the photon and the electron after a Compton diffusion.

    Parameters:
    E0 (float): Energy of the incident photon in keV

    Returns:
    E_retro (float): Energy of the photon after retro-diffusion
    E_fc (float): Energy of the electron after retro-diffusion
    """
    
    for teta in range(360):
        En = calculate_energy(E0, teta)
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

def escape_peaks(E0) : 
    # add conventionnal comments
    E_single_escape = E0 - EMCE
    E_double_escape = E0 - EMCE * 2
    return E_single_escape, E_double_escape
    
def probabilite_photon_atteint(S, d):
    """
    Calculates the probability that a photon reaches a surface S at a distance d
    from an isotropic source emitting photons in all directions.

    Parameters:
    S (float): surface area of the detector (in square meters)
    d (float): distance between the source and the detector (in meters)

    Returns:
    float: probability that the photon reaches the detector
    """
    P = S / (4 * np.pi * d**2)
    return P


def estimated_proportion(energy, Z, plot=True):
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
    if energy < 100 :  # < 100 keV
        photo_prop = 0.8 + 0.2 * (Z / 100)  # Increases with Z
        compton_prop = 0.2
        pair_prop = 0.0
    elif energy < 1022:  # Between 100 keV and 1 MeV
        photo_prop = 0.2 * (Z / 100)
        compton_prop = 0.8 - 0.2 * (Z / 100)
        pair_prop = 0.0
    else:  # E > 1.022 MeV (pair production possible)
        photo_prop = 0.1 * (Z / 100)
        compton_prop = 0.5 - 0.1 * (Z / 100)
        pair_prop = 0.4 + 0.6 * (Z / 100)  # Increases strongly with Z

    total = photo_prop + compton_prop + pair_prop
    proportions = {
        'photoelectric': photo_prop / total,
        'compton': compton_prop / total,
        'pair_production': pair_prop / total
    }

    # Plot the proportions of photon-electron interactions
    if plot:
        plt.figure("Proportions of Photon-Electron Interactions")
        plt.pie(proportions.values(), labels=proportions.keys(), autopct='%1.1f%%')
        plt.title("Proportions of Photon-Electron Interaction Effects")
        plt.show()
    return proportions

def theorical_spectrum(proportion, E0):
    energies = np.linspace(0, E0 * 1.2, 1000)
    theorical_spectrum = np.zeros_like(energies)
    theorical_spectrum[np.argmin(np.abs(energies - E0))] = 1
    if proportion['compton'] != 0:
        Eretro = calculate_energy(E0, 180)
        Efc = E0 - Eretro
        theorical_spectrum[np.argmin(np.abs(energies - Efc))] = 1
        theorical_spectrum[np.argmin(np.abs(energies - Eretro))] = 1        
        
    if proportion['pair_production'] != 0:
        E_single_escape, E_double_escape = escape_peaks(E0)
        theorical_spectrum[np.argmin(np.abs(energies - E_single_escape))] = 1
        theorical_spectrum[np.argmin(np.abs(energies - E_double_escape))] = 1
    # show the spectrum
    plt.figure("Theorical Spectrum")
    plt.plot(energies, theorical_spectrum, label="Theorical Spectrum", color = "black")
    plt.axvline(x=E0, color='r', linestyle='--', label=f'Photoelectric Effect: {E0:.2f} keV')
    if proportion['compton'] != 0:
        plt.axvline(x=Efc, color='g', linestyle='--', label=f'Compton Front: {Efc:.2f} keV')
        plt.axvline(x=Eretro, color='b', linestyle='--', label=f'Retrodiffusion: {Eretro:.2f} keV')
    if proportion['pair_production'] != 0:
        plt.axvline(x=E_single_escape, color='m', linestyle='--', label=f'Single Escape Peak: {E_single_escape:.2f} keV')
        plt.axvline(x=E_double_escape, color='c', linestyle='--', label=f'Double Escape Peak: {E_double_escape:.2f} keV')
    plt.title("Theorical Spectrum")
    plt.ylabel("Intensity")
    plt.xlabel("Energy [keV]")
    plt.legend()
    plt.show()
    return theorical_spectrum

def plot_IRM():
    # plot energy of photoelectric effect and Compton scattering in function of incident energy of the photon
    energies = np.linspace(0, 2000, 2000)
    # create a 2D array with the axis x = incident energy and y = energy of the several interactions
    IRM = np.zeros((len(energies), 3))
    for i, energy in enumerate(energies):
        proportion = estimated_proportion(energy, 32, plot=False)
        IRM[i, 0] = energy
        # retrodiffusion
        IRM[i, 1] = calculate_energy(energy, 180) if proportion['compton'] != 0 else energy
        # Compton front
        IRM[i, 2] = energy - IRM[i, 1] if proportion['compton'] != 0 else IRM[i, 1]
        
    plt.figure("IRM")
    plt.plot(IRM[:, 0], IRM[:, 1], label="Compton Back scattering")
    plt.plot(IRM[:, 0], IRM[:, 2], label="Compton front")
    plt.plot(IRM[:, 0], IRM[:, 0], label="Photoelectric Effect")
    plt.title("Incident Energy vs Interaction Energy")
    plt.ylabel("Interaction Energy [keV]")
    plt.xlabel("Incident Energy [keV]")
    plt.legend()
    plt.show()
    # add intensity of the spectrum
    # add pair creation into the spectrum
    # Z dependency of the proportion calculation
    

    return IRM


    
    