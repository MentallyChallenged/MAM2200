import numpy as np
import matplotlib.pyplot as plt
# FikspunktIterasjon.py: Skript som beregner nullpunkt ved fikspunktiterasjon.
# Definer funksjonen g(x) i ligninga g(x)=x som vi skal bestemme nullpunktet til:
def g(x):
    return (1/np.exp(x))
# Sett startpunktet x0 :
x0=1
# Definerer hjelpevariabel "avvik":
avvik=1
antall = 0
while avvik>0.0001:
# Regn ut neste punkt i iterasjonen:
    xneste=g(x0)
# Regn ut avviket mellom to påfølgende punkter:
    avvik=np.abs(x0-xneste)
    x0=xneste
# Skriv til skjerm
    print(xneste)
    antall +=1
print(antall)