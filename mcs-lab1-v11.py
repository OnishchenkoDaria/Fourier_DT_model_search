import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

def FDT_analysis():
  yt, N, T, t_values, delta_t = set_metrics()
  #calling Discrete Fourier Transform function
  cy = Fourier_DT_function(yt, N)

  peaks = finding_peaks(cy, N)
  
  #Finding the frequencies values of the peaks
  delta_f = 1 / T
  frequency = delta_f * peaks
  print(f"Frequencies with the largest contribution: {frequency}")

def finding_peaks(cy, N):
  #partening the original spectre into half
  half_N = math.ceil(N / 2)
  cy_half = cy[:half_N]

  #visualizing the half 
  plt.plot(np.abs(cy_half))
  plt.title("Half of the FDT graph with true values")
  plt.show()

  #finding absolute spectre peaks
  # _ - secondory unnecessary for analysis return values
  peaks, _ = find_peaks(np.abs(cy_half))
  print(f"Max indexes: {peaks}")
  return peaks

def set_metrics():
  #reading observations from the file, initializing y array
  yt = read_numbers_from_file("f11.txt")
  N = len(yt)
  
  #initializing constants
  T = 5

  #initializing t array
  t_values = np.linspace(0, T, N) # N - for the t array to be divided eaually to the y array
  delta_t = t_values[1] - t_values[0]

  #visualizing the frequency plot
  plt.plot(yt, t_values)
  plt.title('y(t) graph')
  plt.show()
  return yt, N, T, t_values, delta_t 

def Fourier_DT_function(yt, N):
  
  #initializing empty complex array for the results
  cy = np.zeros(N, dtype=complex) # from 0 to N-1
  
  for k in range(N):  
    sum_value = 0
    for m in range(N):  #the Sum volume loop
       sum_value += yt[m] * np.exp(-1j * 2 * np.pi * k * m / N)  
       
       # Alternative way to find the exp value
       # e^(-i*2*pi*k*m / N)
       # e^i*f = cos(f) + i*sin(f)
       # f = -2*pi*k*m / N 
    cy[k] = sum_value / N

  #visualization
  plt.plot(np.abs(cy))
  plt.title('Graph of absolute FDT values')
  plt.show()
  return cy

def read_numbers_from_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
        numbers = list(map(float, content.split()))
    return numbers

if __name__ == "__main__":
  FDT_analysis()  