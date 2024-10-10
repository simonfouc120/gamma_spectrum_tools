# -*- coding: utf-8 -*-

""" 
@author : Simon Foucambert
@Date   : 09/10/2024
@Nom    : spectrum_tools
"""

from library_spectrum_tools import *

S = 1.0  # detector surface in m²
d = 2.0  # Distance between source and detector in m
energy = 100  # Unit : keV
Z = 32  # Atomic number of the material


if __name__ == "__main__":
    geometric_eff = probabilite_photon_atteint(S, d)
    print("Geometric efficacity is {:.2f} %".format(geometric_eff * 100))
    E_retro, E_fc = compton_diff(energy)
    print("Energy of retrodiffusion : {:.2f} keV".format(E_retro))
    print("Energy of Compton front : {:.2f} keV".format(E_fc))
    proportion = estimated_proportion(energy, Z)
    print("Proportion of photoelectric effect : {:.2f}".format(proportion['photoelectric']))
    print("Proportion of Compton effect : {:.2f}".format(proportion['compton']))
    print("Proportion of pair production : {:.2f}".format(proportion['pair_production']))
    theorical_spectrum(proportion, energy)

# to do : modify intensity of the spectrum
# to do : add pair creation into the spectrum
# to do : search library to improve proportion calculation
