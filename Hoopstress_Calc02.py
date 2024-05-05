import matplotlib.pyplot as plt
import numpy as np


# Definierte Werte
#   Konstanten des Versuchsaufbau
nu = 0.3                # [-] possion's ratio
#p = 1 - nu # [-] const.

#   Spannungsrandbedingungen
#s_z0 = 0    *  (10**6)  # [Pa] !!!Wert frei gewählt
#s_ra = 10   *  (10**6)  # [Pa] !!!Wert frei gewählt
#s_ri = 10   *  (10**6)  # [Pa] !!!Wert frei gewählt

r_i = 0.430 # [m] innerer radius
r_a = 0.646 # [m] äußerer radius
r = np.linspace(r_i,r_a,216) # array mit diskreten Radien

j = 114.2   * (1000**2) # [A/m^2] Stromdichte
b_za = -2.5             # [T] magnetische Flussdichte außen
b_zi = 15               # [T] magnetische Flussdichte innen
b_0 = b_za - (b_za-b_zi)/(r_a-r_i) * r_a # [T] absolutes Glied der Geradengleichung für B(r)

#def calcBFeld(r, r_a, r_i, b_za, b_zi, b_0):
#    return ((b_za - b_zi)/(r_a - r_i))  *  r + b_0


def calcIntegral(n, r, r_End, r_Start, j, b_za, b_zi, b_0):
    return j * (  (b_za - b_zi)/(r_End - r_Start) * (1 / (n+2) ) * (r**(n+2) - r_i**(n+2)) 
                 +  b_0                     * (1/  (n+1) ) * (r**(n+1) - r_i**(n+1)))
    
def calcStresses(s_z0, s_ra, s_ri):
    ### Berechne s_z
    s_z = s_z0  # from eq (4)
    
    ### Berechne s_phii
    s_phii = (  2 / (1 - (r_i**2/r_a**2))  ) * ( s_ra  +  ((1 + nu) / 2) * calcIntegral(0, r_a, r_a, r_i, j, b_za, b_zi, b_0)
                                                +  (1/r_a**2) * (1 - (1 + nu)/2 ) * calcIntegral(2,r_a, r_a, r_i, j, b_za, b_zi, b_0)
                                                -  (s_ri  *  (1/2)  *  (1 + (r_i**2/r_a**2))))
    
    ###Berechne s_r(r)
    s_rArray = (  (1/2) * s_phii * (1 - r_i**2/r**2) 
                - ((1 + nu) / 2) * calcIntegral(0, r, r_a, r_i, j, b_za, b_zi, b_0)
                -  (1/r**2) * (1 - (1 + nu)/2 ) * calcIntegral(2, r, r_a, r_i, j, b_za, b_zi, b_0)
                +  (s_ri  *  (1/2)  *  (1 + (r_i**2/r**2))))
    
    ### Berechne s_phi
    s_phiArray = ( (s_ri - nu * s_z0)  
                 +  nu * s_z  
                 -  (nu + 1) * calcIntegral(0, r, r_a, r_i, j, b_za, b_zi, b_0)
                 + s_ri
                 - s_rArray )
    
    return[s_rArray, s_phiArray]

sArray = []
sArray.insert(len(sArray), calcStresses(0, 0, 0))
sArray.insert(len(sArray), calcStresses(0, 0, 5 * 10**6))
sArray.insert(len(sArray), calcStresses(0, 0, -5 * 10**6))

fig, (axs0, axs1) = plt.subplots(2,3, figsize=(12,15),sharex="col", sharey="row")
axs0[0].set_ylabel('radial stress [Pa]')
axs1[0].set_ylabel('hoop stress [Pa]')

i=0
for ax in axs0:
    ax.plot(r, sArray[i][0], label='radial stresses over radius', linewidth=2)
    ax.grid(True)
    #ax.set_ylim([-100, 100])
    i += 1



i=0
for ax in axs1:
    ax.plot(r, sArray[i][1], label='hoop stresses over radius', linewidth=2)
    ax.grid(True)
    ax.set_xlabel('radius [m]')
    i += 1
    
plt.show()

    