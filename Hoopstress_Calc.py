import matplotlib.pyplot as plt
import numpy as np

# Berechne sigma_z

s_z = 0 #

# Berechne sigma_phi_i

def term(n,r_2, r_1):
    """nützlicher Term"""
    return (1/n) * (r_2**n - r_1**n)

r_i = 430 # [mm] innerer radius
r_a = 646 # [mm] äußerer radius

r = np.linspace(r_i,r_a,216)

nu = 0.3 # [-] possion's ratio
p = 1 - nu # [-] const.

j = 114.2 # [A/mm^2] Stromdichte
b_za = -2.5 # [T] magnetische Flussdichte außen
b_zi = 15 # [T] magnetische Flussdichte innen
b_0 = 49.84 # [T] absolutes Glied der Geradengleichung für B(r)

# j = 100
# b_za = -1
# b_zi = 10
# b_0 = 31.9

a = j * (b_za - b_zi)/(r_a - r_i) * term(2,r_a, r_i) + j * b_0 * term(1,r_a, r_i)

b = j * (b_za - b_zi)/(r_a - r_i) * term(4,r_a, r_i) + j * b_0 * term(3,r_a, r_i)

s_phi_i = ( (2/(1-(r_i**2 / r_a**2))) * ((p/2) * a + 1/r_a**2 * (1 - p/2) * b) ) / 1000 # [MPa = N/mm^2]

print(s_phi_i)

# Berechne sigma_r(r)

s_r = 1/2 * s_phi_i * (1 - r_i**2/r**2) - 1/2 * p * ((j * (b_za - b_zi)/(r_a - r_i) * term(2,r, r_i) + j * b_0 * term(1,r, r_i)) / 1000) - (1/r**2 * (1 - p/2) * (j * (b_za - b_zi)/(r_a - r_i) * term(4,r, r_i) + j * b_0 * term(3,r, r_i))) / 1000

print(s_r)

s_phi = s_phi_i - p * ((j * (b_za - b_zi)/(r_a - r_i) * term(2,r, r_i   ) + j * b_0 * term(1,r, r_i)) / 1000) - s_r

f = -75 + 13 * (r - 420)**(1/3)

plt.plot(r, s_phi, label='hoop stresses over radius', linewidth=2)
plt.plot(r, (s_phi+f), label='hoop stresses over radius', linewidth=2)
plt.grid()
plt.xlabel('r, [mm]')
plt.ylabel(r'Sigma_phi, [MPa]')
plt.title('Hoop Stresses over radius')
plt.legend
plt.show()


plt.plot(r, s_r, label='radial stresses over radius', linewidth=2)
plt.grid()
plt.xlabel('r, [mm]')
plt.ylabel(r'Sigma_r, [MPa]')
plt.title('Radial Stresses over radius')
plt.legend
plt.show()
