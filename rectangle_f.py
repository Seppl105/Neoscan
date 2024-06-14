import numpy as np
import matplotlib.pyplot as plt

r_innen = 430 # [mm] innerer Radius
r_aussen = 646 # [m] äußerer Radius

t = 0.36 # [mm] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
mt = np.array([0.12,0.23,0.01])
mE = np.array([500, 500, 450, 450, 400, 400])


# Calculate Number of Windings
windings = lambda x, y, z: (x - y) / z
numberOfWindings = int(windings(r_aussen, r_innen, t))  # should be 600 Windings with given measurements
print("Anzahl der Windungen: ", numberOfWindings)


def calcJumpspot(r_innen, materialWidths, numberOfWindings):
    '''Calculates all radii (results) at which the jumps occur.'''
    thisSpot = r_innen
    results = []
    for i in range(numberOfWindings):
        for width in materialWidths:
            thisSpot += width
            results.append(thisSpot)
            results.append(thisSpot)
    return np.array(results)


jumps = calcJumpspot(r_innen, mt, numberOfWindings)
jumps = np.insert(jumps, 0, 430)
jumps = jumps[:-1]
print(jumps)

numberOfMaterials = len(mE)
print(numberOfMaterials)

# Three YOUNG's Modulus for each winding
E = np.ones(numberOfMaterials * numberOfWindings)
for i in range(numberOfMaterials):
    E[i::numberOfMaterials] = mE[i]
print(E)
print(len(E))
print(len(jumps))

plt.figure(figsize=(8, 6))
plt.plot(jumps, E, label='E(r) in N/mm^2', color='b')
# Lable und Titel hinzufügen
plt.xlabel('r')
plt.ylabel('E(r)')
plt.title('E(r)')
plt.grid(True)
plt.legend()
plt.show()