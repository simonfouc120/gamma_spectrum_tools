# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 09:34:15 2023

@author: sf270338
"""

import math 
import matplotlib.pyplot as plt


def compton_diff(E0,fig) : 
    E_photon=list()
    E_electron=list()
    for teta in range (360):
        En=(float(E0/(1+(E0/511)*(1-math.cos((teta*2*math.pi)/360)))))
        E_photon.append(En)
        E_electron.append(float(E0-En))
        if teta == 180 : 
            print("Photon incident d'énergie : ", E0 , "keV")
            print("Energie correspondante à la rétrodiffusion :", En , "keV")
            print("Energie correspondante au front Compton : ", float(E0-En),"keV")
    plt.figure(fig)
    plt.plot(E_photon, label = "Energie du photon Compton")
    plt.plot(E_electron,label="Energie de l'électron Compton")
    plt.title("Comton Diff for Energy of " + str(E0) +" keV")
    plt.ylabel("Energy [keV]")
    plt.xlabel("Angle of déviation [°]")
    plt.legend()
    plt.show
    

compton_diff(59,1)
compton_diff(300,2)
compton_diff(1333,3)
