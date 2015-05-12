import numpy as np
import matplotlib.pyplot as plt
class String:
  
  def __init__(self,note='C2',*args):
    
    self.delta_t = 1/(16e3)
    
    # Choose constrains
    if str.lower(note) == 'c2':
      # Bass note
      b_1 = 0.5
      b_2 = 6.25e-9
      self.N = 100
      self.L = 1.9
      self.Ms = 35
      self.epsilon = 7.5e-9
      self.T = 750
      self.rho = self.Ms/self.L
      self.c = np.sqrt(self.T/self.rho)
    if str.lower(note) == 'c4':
      # Midrange note
      b_1 = 0.5
      b_2 = 6.25e-9
      self.N = 50
      self.L = 0.62
      self.epsilon = 3.82e-5
    if str.lower(note) == 'c7':
      # Treble note
      b_1 = 0.5
      b_2 = 2.6e-10
      self.N = 16
      self.L = 0.09
      self.epsilon = 8.67e-4
      


    self.delta_x = self.N/self.L
    kappa_squared = self.epsilon*(self.c**2)*(self.L**2)
    mu = kappa_squared/((self.delta_x**2)*(self.c**2))

    D = 1+b_1*self.delta_t
    r = self.c*self.delta_t/self.delta_x
    if r> 1:
        raise ValueError("r is larger than 1, simulation is unstable. r should be smaller than 1.")
    if self.delta_t > (self.delta_x)/self.c:
        raise ValueError("CCFL condition broken")
    
    nu = 2*b_2*(self.delta_t)/(self.delta_x**2)
    

    self.a_1 = (-(r**2)*mu)/D
    self.a_2 = (r**2+4*(r**2)*mu+nu)/D
    self.a_3 = (2-2*(r**2)-6*(r**2)*mu-2*nu)/D
    self.a_4 = (-1+b_1*self.delta_t+2*nu)/D
    self.a_5 = -nu/D

    # y_[plus/minus]_[n/2n] are string displacement vectors
    self.y_plus_n = np.zeros((self.N),dtype = float)
    self.y = np.zeros((self.N),dtype = float)
    self.y_minus_n = np.zeros((self.N),dtype = float)

  def time_evolution(self):

    F = 0.0000
    init_velocity = 1.0
    hammer_length = 0.01
    hammer_position = 0.15*self.L
    hammer_displacement = 0.0
    hammer_displacement_minus_n = -init_velocity*self.delta_t
    hammer_mass = 1.0
    K = 1.
    p = 2.
    g = np.zeros((self.N),dtype = float)

    n = 0
    beginwindow = int(np.floor((hammer_position-hammer_length/2)*self.N/self.L))
    endwindow = int(np.ceil((hammer_position+hammer_length/2)*self.N/self.L))

    for i in range(beginwindow,endwindow):
        g[i] = 1.0/(endwindow-beginwindow)


    for t in range(0,2001):
        # F = K/hammer_mass * np.abs(hammer_displacement - self.y)**p        
        for i in range(2,self.N-2):
          self.y_plus_n[i] = self.a_1*(self.y[i+2]+self.y[i-2])+self.a_2*(self.y[i+1]+self.y[i-1]) \
          + self.a_3*self.y[i]+self.a_4*self.y_minus_n[i]+self.a_5*(self.y_minus_n[i+1]+self.y_minus_n[i-1])
          + (self.delta_t**2*self.N*F*g[i])/self.Ms
        # hammer_displacement = 2*hammer_displacement - hammer_displacement_minus_n + F*delta_t**2
        if t == 1:
            print(t)
            F = 0.5

        if t == 5:
            F = 0

        # Update new string heights to old ones
        self.y_minus_n[:] = self.y[:]
        self.y[:] = self.y_plus_n[:]

        # plotting. Only used for testing
        x = np.linspace(0,self.L,self.N)
        plt.plot(x,self.y)
        plt.axis([0,self.L,-0.01,0.01])
        plt.show()
