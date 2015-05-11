import numpy as np

class String:
  
  def __init__(self,note='C2',*args):
    self.c = 343
    self.delta_t = 0.1
    
    # Choose constrains
    if str.lower(note) == 'c2':
      # Bass note
      b_1 = 0.5
      b_3 = 6.25e-9
      self.N = 100
      self.L = 1.9
      self.epsilon = 7.5e-9
    if str.lower(note) == 'c4':
      # Midrange note
      b_1 = 0.5
      b_3 = 6.25e-9
      self.N = 50
      self.L = 0.62
      self.epsilon = 3.82e-5
    if str.lower(note) == 'c7':
      # Treble note
      b_1 = 0.5
      b_3 = 2.6e-10
      self.N = 16
      self.L = 0.09
      self.epsilon = 8.67e-4
      

    self.D = 1+b_1*self.delta_t+(2*b_3/self.delta_t)
    self.delta_x = self.N/self.L
    self.r = self.c*self.delta_t/self.delta_x
    
    self.a_1 = (2-2*self.r**2+(b_3/self.delta_t)-6*self.epsilon*self.N**2*self.r**2)/self.delta_t
    self.a_2 = (-1+b_1*self.delta_t+2*b_3/self.delta_t)/self.D
    self.a_3 = ((self.r**2)*(1+4*self.epsilon*self.N**2))/self.D
    self.a_4 = (b_3/self.delta_t-self.epsilon*self.N**2*self.r**2)/self.D
    self.a_5 = (-b_3/(self.delta_t*self.D))

    # y_[plus/minus]_[n/2n] are string displacement vectors
    self.y_plus_n = np.zeros((self.N),dtype = float)
    self.y = np.zeros((self.N),dtype = float)
    self.y_minus_n = np.zeros((self.N),dtype = float)
    self.y_minus_2n = np.zeros((self.N),dtype = float)

  def time_evolution(self):
    F = 0
    hammer_length = 0.01
    hammer_position = 0.15*self.L
    Ms = 1
    g = np.zeros((self.N),dtype = float)
    n = 0
    beginwindow = int(np.floor((hammer_position-hammer_length/2)*self.N/self.L))
    endwindow = int(np.ceil((hammer_position+hammer_length/2)*self.N/self.L))
    
    for i in range(beginwindow,endwindow):
        g[i] = 1.0/(endwindow-beginwindow)
    
    for i in range(2,self.N-2):
      print(i)
      self.y_plus_n[i] = self.a_1*self.y[i] + self.a_2*self.y_minus_n[i] + self.a_3*(self.y[i+1]+self.y[i-1]) \
      + self.a_4*(self.y[i+2]+self.y[i-2]) + self.a_5*(self.y_minus_n[i+1]+self.y_minus_n[i-1]+self.y_minus_2n[i]) \
      + (self.delta_t**2*self.N*F*g[i])/Ms
    