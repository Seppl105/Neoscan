import matplotlib.pyplot as plt
import numpy as np

# Definierte Werte
#   Konstanten des Versuchsaufbau
nu = 0.3 # [-] possion's ratio
p = 1 + nu # [-] const.

#   Spannungsrandbedingungen
#s_z0 = 0 #  [MPa] !!!Wert frei gewählt
#s_ra = 10 # [MPa] !!!Wert frei gewählt
#s_ri = 10 # [MPa] !!!Wert frei gewählt

r_i = 430 # [mm] innerer radius
r_a = 646 # [mm] äußerer radius
r = np.linspace(r_i,r_a,216) # array mit diskreten Radien

j = 114.2 # [A/mm^2] Stromdichte
b_za = -2.5 # [T] magnetische Flussdichte außen
b_zi = 15 # [T] magnetische Flussdichte innen
b_0 = b_za - (b_za-b_zi)/(r_a-r_i) * r_a # [T] absolutes Glied der Geradengleichung für B(r)

# j = 100
# b_za = -1
# b_zi = 10
# b_0 = 31.9


# Funktionen
def term(n,r_2, r_1):
    #wiederkehrender Term
    return (1/n) * (r_2**n - r_1**n)



def calcStresses(s_z0, s_ra, s_ri):
    # Berechne sigma_z

    s_z = s_z0 # from eq (4)

    ### Berechne sigma_phi_i

    a = j * (b_za - b_zi)/(r_a - r_i) * term(2,r_a, r_i)  +  j * b_0 * term(1,r_a, r_i) # Integral (I)

    b = j * (b_za - b_zi)/(r_a - r_i) * term(4,r_a, r_i)  +  j * b_0 * term(3,r_a, r_i) # Integral (II)

    s_phi_i = (  (2/(1-(r_i**2 / r_a**2)))    *    ((p/2) * a + 1/r_a**2 * (1 - p/2) * b)  + 1000 * s_ra - 1/2 * (1 + r_i**2/r_a**2) * 1000 * s_ri) / 1000 # [MPa = N/mm^2]

    print(s_phi_i)

    ### Berechne sigma_r(r)

    s_r = 1/2 * s_phi_i * (1 - r_i**2/r**2) - 1/2 * p * ((j * (b_za - b_zi)/(r_a - r_i) * term(2,r, r_i) + j * b_0 * term(1,r, r_i)) / 1000) - (1/r**2 * (1 - p/2) * (j * (b_za - b_zi)/(r_a - r_i) * term(4,r, r_i) + j * b_0 * term(3,r, r_i))) / 1000 + 1/2 * (1 + r_i**2/r**2) * s_ri

    # print(s_r)
    
    ### Berechne sigma_phi

    s_phi = s_phi_i - p * ((j * (b_za - b_zi)/(r_a - r_i) * term(2,r, r_i   ) + j * b_0 * term(1,r, r_i)) / 1000) - s_r

    return [s_r, s_phi]

#f = -75 + 13 * (r - 420)**(1/3)
sArray = []
sArray.insert(len(sArray), calcStresses(0, 0, 0))
sArray.insert(len(sArray), calcStresses(0, 0, 5))
sArray.insert(len(sArray), calcStresses(0, 0, -5))

fig, (axs0, axs1) = plt.subplots(2,3, figsize=(12,15),sharex="col", sharey="row")
axs0[0].set_ylabel(r'radial stress [MPa]')
axs1[0].set_ylabel(r'hoop stress [MPa]')

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
    ax.set_xlabel('radius [mm]')
    #ax.set_ylim([100, 450])
    i += 1

# fig, (ax1, ax2) = plt.subplots(1,2)
# ax1.plot(r, s_phi, label='hoop stresses over radius', linewidth=2)
# ax1.plot(r, (s_phi), label='hoop stresses over radius', linewidth=2)
# ax1.grid()
# ax1.set_xlabel('r, [mm]')
# ax1.set_ylabel(r'Sigma_phi, [MPa]')
# ax1.set_title('Hoop Stresses over radius')
# #ax1.legend
# #plt.show()

# ax2.plot(r, s_r, label='radial stresses over radius', linewidth=2)
# ax2.plot(r, (b_za-b_zi)/(r_a-r_i) * r + b_0)
# ax2.grid()
# ax2.set_xlabel('r, [mm]')
# ax2.set_ylabel(r'Sigma_r, [MPa]')
# ax2.set_title('Radial Stresses over radius')
# #ax2.legend

plt.show()




#print(np.linspace(0,len(sArray) - 1,len(sArray)))
# for i in np.linspace(0,len(sArray) - 1,len(sArray)):
#     axs[i, 0].plot(r, sArray[i][0], label='hoop stresses over radius', linewidth=2)
#     axs[i, 1].plot(r, sArray[i][1], label='radial stresses over radius', linewidth=2)