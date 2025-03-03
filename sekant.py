import numpy as np
import math

# Definer funksjonen f(x):
def f(x):
    return x*np.exp(x)-1
# Startpunktene for sekantmetoden:
x0 = 1
x1 = 0.6
# Toleranse for konvergens:
toleranse = 0.0001
# Teller antall iterasjoner:
antall = 0
# Sekantiterasjon:
while True:
    # Beregn neste punkt med sekantmetoden:
    xneste = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    # Beregn avvik:
    avvik = np.abs(xneste - x1)
    # Oppdater punktene:
    x0 = x1
    x1 = xneste
    # Tell iterasjonen:
    antall += 1
    # Avslutt hvis avviket er mindre enn toleransen:
    if avvik < toleranse:
        break
    # Skriv ut mellomresultater (valgfritt):
    print(f"Iterasjon {antall}: x = {xneste:.6f}")

# Skriv sluttresultatet:
print(f"Løsning funnet etter {antall} iterasjoner: x = {xneste:.6f}")
