# import numpy as np

# array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# print(array[-5:])


# Er = np.ones(length) * 100 * 10**9 # E Modul
# print(len(Er[:5]))
# print(len(Er[-5:]))
# print(len(Er[-50:]))

import numpy as np
import matplotlib.pyplot as plt

# Sample data points (you can replace this with your actual data)
length = 100
data = np.ones(length) * 100 * 10**9 # E Modul
data[:int(0.6*length)] = 150* 10**9



#data = np.sin(0.03* np.linspace(0,length,length)) + np.cos(0.09* np.linspace(0,length,length))
#print(data)
n = len(data)
t = np.linspace(0, length, n, endpoint=False)  # Time vector (assuming sampling over 1 second)
tManyValues = np.linspace(0, length, length*5)

# Compute the FFT
coefficients = np.fft.fft(data)

# Frequency components
frequencies = np.fft.fftfreq(n)

# Create the function from the Fourier coefficients
def fourier_series(t, coefficients, frequencies):
    # Initialize the result with the mean value
    result = np.real(coefficients[0]) / n
    
    for k in range(1, n//2):  # We only need the first half of the frequencies due to symmetry
        amplitude = 2 * np.abs(coefficients[k]) / n  # Get the amplitude
        phase = np.angle(coefficients[k])  # Get the phase
        result += amplitude * np.cos(2 * np.pi * frequencies[k] * t + phase)
        
    return result

# Reconstruct the function
reconstructed_function = fourier_series(tManyValues, coefficients, frequencies)

# a = range(10)
# print(a)
# print(a[1:])
# print(a[:-1])

# # Plot the original data and the reconstructed function
# plt.figure(figsize=(12, 6))
# plt.plot(t, data, label='Original Data')
# plt.plot(tManyValues, reconstructed_function, label='Reconstructed Function', linestyle='--')
# plt.legend()
# plt.show()
b = np.ones(5)
x = np.ones(7)
c = np.array(range(5))
print(b.size)
print(x.size)
print(np.outer(x, b))
print(np.outer(x, b) + c)
print(np.sum(np.outer(x, b) + c, axis=1))