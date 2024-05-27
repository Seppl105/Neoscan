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

anzahlPlots = 5

plt.subplot(1, anzahlPlots, 1)
plt.plot(Er)
plt.xlabel("Radius in m")
plt.ylabel("E-Modul in N/m^2")
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
fftEr = np.fft.fft(Er)  # aus dem fft package die fft funktion
# rrster Plot Realteil
plt.subplot(1, anzahlPlots, 2)
plt.plot(np.fft.fftshift(np.abs(fftEr)))   # fftshift() ist sinnvoll da die erste Hälfte der Werte für die aufsteigende Frequenz berechnet wird und dann zur negativen maximalen Frequenz gesprungen wird und dann steigt es wieder zu 0
# jetzt startet die Frequenz bei der größt negativen und steigt dann zur positiven maximalen
plt.title("Magnitude of the FT of E(r)")
plt.ylabel("Magnitude")
plt.xlabel("Frequency Index")
# zweiter Plot Imaginärteil
plt.subplot(1, anzahlPlots, 3)
plt.plot(np.fft.fftshift(np.angle(fftEr)))
plt.title("Plot Angle of the FT of E(r)")
plt.ylabel("Magnitude")
plt.xlabel("Frequency Index")

plt.subplot(1, anzahlPlots, 4)
ErTest = (np.fft.fftshift(np.real(fftEr)) * np.cos(2 * np.pi * 100 *  np.fft.fftfreq(len(Er)))
         + np.fft.fftshift(np.imag(fftEr)) * np.sin(2 * np.pi * 100 * np.fft.fftfreq(len(Er))))
plt.plot(ErTest)

plt.subplot(1, anzahlPlots, 5)
ErResult = np.fft.ifft(fftEr)
plt.plot(ErResult)

plt.show()