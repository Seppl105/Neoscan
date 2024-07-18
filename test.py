import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime # Nur für die Bennenung der Grafiken
import math as m
#from delete.materials import calcMaterialTanhCoefficients
#from delete.materials import calcMaterialValue

# # import numpy as np

# # array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# # print(array[-5:])


# # Er = np.ones(length) * 100 * 10**9 # E Modul
# # print(len(Er[:5]))
# # print(len(Er[-5:]))
# # print(len(Er[-50:]))



# #############################################################################

# # Eigenschaften der Spule
# r_i = 430 # [mm] innerer Radius
# r_a = 646 # [mm] äußerer Radius
# t = 0.36 # [mm] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
# t_con = 0.12 # [mm] Dicke des Bandleiters
# t_cow = 0.23 # [mm] Dicke des Cowindings
# t_ins = 0.01 # [mm] Dicke der Isolation
# materialsWidth = [t_con, t_cow, t_ins]

# E_con = 500 #* 10**9 # E-Modul Conductor
# E_cow = 450 #* 10**9 # E-Modul Cowinding
# E_ins = 400 #* 10**9 # E-Modul Insulation
# materialsE = [500, 450, 400]

# # Anteil der drei Materialien pro Wicklung bestimmen
# anteilConductor = t_con / t
# anteilCowinding = t_cow / t
# anteilInsulation = t_ins / t

# #for i in range(5):
# #    print(i)


# r,E,cof  = calcMaterialTanhCoefficients(r_i,r_a,t,[t_con,t_cow + t_ins],[E_con,E_cow],steigung=1000)
# print("Erste Berechnung abgeschllossen")
# # Plot E(r) für alle Windungen
# plt.figure(figsize=(8, 6))
# plt.plot(r, E, label='E(r) in N/mm^2', color='b')
# plt.plot(r, calcMaterialValue(r, cof, materialsE[0]))
# plt.scatter(434, calcMaterialValue(434, cof, materialsE[0]))
# plt.scatter(434.3245, calcMaterialValue(434.3245, cof, materialsE[0]))
# plt.scatter(432.23, calcMaterialValue(432.23, cof, materialsE[0]))
# plt.scatter(434.32, calcMaterialValue(436.32, cof, materialsE[0]))
# # Lable und Titel hinzufügen
# plt.xlabel('r')
# plt.ylabel('E(r)')
# plt.title('E(r) durch Tangens Hyperbolicus genähert')
# plt.grid(True)
# plt.legend()
# plt.show()


# #############################################################################

# # # Sample data points (you can replace this with your actual data)
# # length = 100
# # data = np.ones(length) * 100 * 10**9 # E Modul
# # data[:int(0.6*length)] = 150* 10**9



# # #data = np.sin(0.03* np.linspace(0,length,length)) + np.cos(0.09* np.linspace(0,length,length))
# # #print(data)
# # n = len(data)
# # t = np.linspace(0, length, n, endpoint=False)  # Time vector (assuming sampling over 1 second)
# # tManyValues = np.linspace(0, length, length*5)

# # # Compute the FFT
# # coefficients = np.fft.fft(data)

# # # Frequency components
# # frequencies = np.fft.fftfreq(n)

# # # Create the function from the Fourier coefficients
# # def fourier_series(t, coefficients, frequencies):
# #     # Initialize the result with the mean value
# #     result = np.real(coefficients[0]) / n
    
# #     for k in range(1, n//2):  # We only need the first half of the frequencies due to symmetry
# #         amplitude = 2 * np.abs(coefficients[k]) / n  # Get the amplitude
# #         phase = np.angle(coefficients[k])  # Get the phase
# #         result += amplitude * np.cos(2 * np.pi * frequencies[k] * t + phase)
        
# #     return result


# # Reconstruct the function
# #reconstructed_function = fourier_series(tManyValues, coefficients, frequencies)

# # a = range(10)
# # print(a)
# # print(a[1:])
# # print(a[:-1])

# # # Plot the original data and the reconstructed function
# # plt.figure(figsize=(12, 6))
# # plt.plot(t, data, label='Original Data')
# # plt.plot(tManyValues, reconstructed_function, label='Reconstructed Function', linestyle='--')
# # plt.legend()
# # plt.show()
# b = np.ones(5)
# x = np.ones(7)
# c = np.array(range(5))
# print(b.size)
# print(x.size)
# print(np.outer(x, b))
# print(np.outer(x, b) + c)
# print(np.sum(np.outer(x, b) + c, axis=1))

# # currentTime = datetime.now()

# # # currentTime = str(currentTime.strftime("%H-%M%S"))[0:5]
# # # pictureName = f"Bilder\Graphen{currentTime}.png"

# # pictureName = f"Bilder\Graphen{currentTime.year}-{currentTime.month:02d}-{currentTime.day:02d}-{currentTime.hour:02d}-{currentTime.minute:02d}"

# # plt.plot(np.linspace(0, 10, 100), np.linspace(0, 10, 100))
# # plt.savefig(pictureName)
# # plt.show()
a = np.array([4,4,3])
if a.any():
    print("True")
else:
    print("false")

from math import log10, floor # floor returns the greatest integer

number = 3.2345 * 10**26
print(number/ 10**(abs(floor(log10(abs(number))))))

#print(number / 10**(int(number * 10**26 / 10)))
t = [3, 3, 3, 3, 3]
print([0] * (len(t) - 1))

#print(m.log10(int(0.000000000000000000000000000000000000000000000001)))

print(t[:-1])