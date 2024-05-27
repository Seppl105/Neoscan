import numpy as np
import matplotlib.pyplot as plt


length = 100

Er = np.ones(length) * 100 * 10**9 # E Modul
Er[:int(0.6*length)] = 150* 10**9

# # imaginärer Teil verschwindet
# Er[:5] = 150* 10**9 # Test für length = 100
# Er[-4:] = 150* 10**9  # Test für length = 100

# # realer Teil verschwindet
# Er = np.zeros(length)
# Er[:5] = 150 * 10**9 # Test für length = 100
# #Er[-50:] = -100 * 10**9 # Test für length = 100
# Er[-4:] = -150 * 10**9  # Test für length = 100

# ############################### Warum geht das nicht?
# Er[:5] = 150 * 10**9 # Test für length = 100
# Er[-50:] = -100 * 10**9 # Test für length = 100
# Er[-5:] = -150 * 10**9  # Test für length = 100

# Fuktion (/Fourierreihe) aus den Fourierkoeffizienten bestimmen
def inverseFourier(r, coefficients, frequencies):
    n = len(coefficients)
    result = np.real(coefficients[0]) / n # konstanter Anteil
    
    for k in range(1, n//2):  # aus symmetriegründen nur die erste Hälfte der Freuqenzen zu betrachten
        amplitude = 2 * np.abs(coefficients[k]) / n  # Amplitude berechnen
        phase = np.angle(coefficients[k])  # Phase berechnen
        result += amplitude * np.cos(2 * np.pi * frequencies[k] * r + phase)
        
    return result

fftEr = np.fft.fft(Er)  # aus dem fft package die fft funktion

anzahlRowPlots = 2
anzahlColumnPlots = 3

plt.subplot(anzahlRowPlots, anzahlColumnPlots, 1)
plt.plot(Er, label="Original Fuction")
plt.xlabel("Radius in m")
plt.ylabel("E-Modul in N/m^2")

#plt.subplot(anzahlRowPlots, anzahlColumnPlots, 1)
plt.plot(np.linspace(0, length, length*5),inverseFourier(np.linspace(0, length, length*5), fftEr, np.fft.fftfreq(length)), label='Reconstructed Function', linestyle='--')
plt.legend()


#plt.show()


# fftEr = np.fft.fft(Er)  # aus dem fft package die fft funktion
# # rrster Plot Realteil
# plt.subplot(1, 3, 2)
# plt.plot(np.fft.fftshift(np.real(fftEr)))   # fftshift() ist sinnvoll da die erste Hälfte der Werte für die aufsteigende Frequenz berechnet wird und dann zur negativen maximalen Frequenz gesprungen wird und dann steigt es wieder zu 0
# # jetzt startet die Frequenz bei der größt negativen und steigt dann zur positiven maximalen
# plt.title("Plot of real part of FT of E(r)")
# plt.ylabel("Magnitude")
# plt.xlabel("Frequency Index")
# # zweiter Plot Imaginärteil
# plt.subplot(1, 3, 3)
# plt.plot(np.fft.fftshift(np.imag(fftEr)))
# plt.title("Plot of imaginary part of FT of E(r)")
# plt.ylabel("Magnitude")
# plt.xlabel("Frequency Index")


# Plot in polar form
# rrster Plot Realteil
plt.subplot(anzahlRowPlots, anzahlColumnPlots, 2)
plt.plot(np.fft.fftshift(np.abs(fftEr)))   # fftshift() ist sinnvoll da die erste Hälfte der Werte für die aufsteigende Frequenz berechnet wird und dann zur negativen maximalen Frequenz gesprungen wird und dann steigt es wieder zu 0
# jetzt startet die Frequenz bei der größt negativen und steigt dann zur positiven maximalen
plt.title("Magnitude of the FT of E(r)")
plt.ylabel("Magnitude")
plt.xlabel("Frequency Index")
# zweiter Plot Imaginärteil
plt.subplot(anzahlRowPlots, anzahlColumnPlots, 3)
plt.plot(np.fft.fftshift(np.angle(fftEr)))
plt.title("Plot Angle of the FT of E(r)")
plt.ylabel("Magnitude")
plt.xlabel("Frequency Index")


plt.subplot(anzahlRowPlots, anzahlColumnPlots, 4)
ErTest = (np.fft.fftshift(np.real(fftEr)) * np.cos(2 * np.pi * 100 *  np.fft.fftfreq(len(Er)))
         + np.fft.fftshift(np.imag(fftEr)) * np.sin(2 * np.pi * 100 * np.fft.fftfreq(len(Er))))
plt.plot(ErTest)



plt.subplot(anzahlRowPlots, anzahlColumnPlots, 5)
ErIfft = np.fft.ifft(fftEr)
plt.plot(ErIfft)

plt.show()