import numpy as np
import math

# Halvering.py: Skript som beregner nullpunkt ved halveringsmetoden.
# Definer funksjonen som vi skal bestemme nullpunktet til:
def f(x):
    return np.exp(x-np.sqrt(x)) - 1
# Sett startintervallet [a,b] hvor vi vet at det er nøyaktig ett nullpunkt:.
# math.sqrt(x+1) + x**(1/3)-2
a=0
b=3
antall = 0
# Sjekk om det er ulike fortegn i endepunktene:
if f(a)*f(b) > 0:
# Skriv beskjed til skjerm:
    print('Ikke ulike fortegn i endepunktene. Stopper.')
else:
# Fortsetter så lenge intervallet [a,b] ikke er smalt nok
    while b-a >= 0.00001:
# Regn ut midtpunkt i intervallet og bestem det nye intervallet:
        c=(a+b)/2
        antall += 1
        if f(a)*f(c) > 0:
            a=c # Velger (c,b)
        else:
            b=c # Velger (a,c)
        print(c)
print(f"Antall iterasjoner: {antall}")
