import math
import matplotlib.pyplot as plt
import numpy as np

def FDT_analysis():
  yt, N, T, t_values, delta_t = set_metrics()
  #calling Discrete Fourier Transform function
  cy = Fourier_DT_function(yt, N)

  half_N = math.ceil(N / 2)
  peaks = finding_peaks(cy, N)
  
  #Finding the frequencies values of the peaks
  delta_f = 1 / T
  frequency = [delta_f * peak for peak in peaks]
  print(f"Frequencies with the largest contribution: {frequency}")
  a_coefficient = solve_coefficients(half_N, t_values, frequency, yt)
  [round(a) for a in a_coefficient]
  print(f"a_j ccoefficient: {a_coefficient}")
  #OLS_error(half_N, t_values, frequency, yt, a_coefficient
  show_model(a_coefficient, frequency)

def model_function(t, a, f):
    poly_part = a[0] * pow(t, 3) + a[1] * pow(t, 2) + a[2] * t
    
    sin_sum = 0
    #number of total a coefficients
    if isinstance(f, list):
      k = len(f) + 4
    elif isinstance(f, float):
      k = 4
    
    for i in range(3, k-2):
        sin_sum += a[i] * np.sin(2 * math.pi * f[i-3] * t)   
    return poly_part + sin_sum + a[k-1]

#display the aproximate function
def show_model(a, f):
    if isinstance(f, list):
      k = len(f) + 4
    elif isinstance(f, float):
      k = 4
      
    sin_sum = ""
    for i in range(3, k-2):
      temp = "sin(2* pi * " + str(f[i-3]) + " * t"
      sin_sum  += temp
    print(f"{a[0]} * t^3 + {a[1]} * t^2 + {a[2]} * t + {sin_sum} + {a[k-1]}")

""" def OLS_error(N, t_values, f, yt, a):
    F_acc = 0 
    if isinstance(f, list):
      k = len(f) + 4
    elif isinstance(f, float):
      k = 4  
    for j in range(N):
        model_value = model_function(t_values[j], a, f)
        result = model_value - yt[j]
        F_acc += 0.5 * pow(result, 2) """

#function to compute partial derivatives
def solve_coefficients(N, t_values, f, yt):
    if isinstance(f, list):
      k = len(f) + 4
    elif isinstance(f, float):
      k = 4     
    A = np.zeros((k, k))  #coefficients of the system
    B = np.zeros(k)  #vector
    
    for j in range(N):
        t_j = t_values[j]
        y_j = yt[j]
        
        poly_terms = np.array([pow(t_j, 3), pow(t_j, 2), t_j, 1])
        
        sin_terms = np.array([math.sin(2 * math.pi * f[i-3] * t_j) for i in range(3, k-1)])
        
        #full model terms (polynomial + sinusoidal)
        terms = np.concatenate([poly_terms, sin_terms, [1]])
        
        #update matrix A and vector B for each equation
        for i in range(k):
            for l in range(k):
                A[i, l] += terms[i] * terms[l]
            B[i] += terms[i] * y_j

    a_coefficients = np.linalg.solve(A, B)
    rounded_coefficients = np.round(a_coefficients).astype(int)

    equation = model_function(t_values, rounded_coefficients, f)
    plt.plot(t_values, yt)
    plt.plot(t_values, equation)
    plt.show()
    return rounded_coefficients

def finding_peaks(cy, N):
  #partening the original spectre into half
  half_N = math.ceil(N / 2)
  cy_half = np.abs(cy[:half_N])

  #finding absolute spectre peaks
  peaks = []
  for i in range(1, len(cy_half)-1):    
    if round(cy_half[i], 4) > round(cy_half[i+1], 4) and round(cy_half[i], 4) > round(cy_half[i-1], 4):
       peaks.append(i)
  print(f"Max indexes: {peaks}")

  #visualizing the half 
  plt.plot(cy_half)
  plt.title("Half of the FDT graph with true values")
  plt.scatter(peaks, cy_half[peaks], color='red')
  plt.show()

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
  plt.plot(t_values, yt)
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