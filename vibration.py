import numpy as np

class String:
  
  def __init__(self,note='C2',*args):

    self.c = 343
    
    N,L,b_1,b_3,F,Ms,g
    
    self.N = N
    self.delta_t = 0.1
    self.epsilon = 1
    self.D = 1+b_1*self.delta_t+(2*b_3/self.delta_t)
    self.delta_x = N/L
    self.r = self.c*self.delta_t/self.delta_x
    
    # Choose constrains
    if str.lower(note) == 'c2':
      # Bass note
      self.N = N
      self.delta_t = 0.1
      self.epsilon = 1
      self.D = 1+b_1*self.delta_t+(2*b_3/self.delta_t)
      self.delta_x = N/L
      self.c = 343
      self.r = self.c*self.delta_t/self.delta_x
    if str.lower(note) == 'c4':
      
    if str.lower(note) == 'c7':
      

    self.a_1 = (2-2*self.r**2+(b_3/self.delta_t)-6*self.epsilon*N**2*self.r**2)/self.delta_t
    self.a_2 = (-1+b_1*self.delta_t+2*b_3/self.delta_t)/self.D
    self.a_3 = ((self.r**2)*(1+4*self.epsilon*N**2))/self.D
    self.a_4 = (b_3/self.delta_t-self.epsilon*N**2*self.r**2)/self.D
    self.a_5 = (-b_3/(self.delta_t*self.D))

    # y_[plus/minus]_[n/2n] are string displacement vectors
    self.y_plus_n = np.zeros((N),dtype = float)
    self.y = np.zeros((N),dtype = float)
    self.y_minus_n = np.zeros((N),dtype = float)
    self.y_minus_2n = np.zeros((N),dtype = float)

  def time_evolution(self):
    F = 0
    hammer_length = 0.01
    hammer_position = 0.15*L
    Ms = 1
    g = np.zeros((N),dtype = float)
    n = 0
    beginwindow = np.floor((hammer_position-hammer_length/2)*N/L)
    endwindow = np.ceil((hammer_position+hammer_length/2)*N/L)
    
    for i in range(beginwindow,endwindow):
        g[i] = 1.0/(endwindow-beginwindow)
    
    for i in range(2,self.N-2):
        self.y_plus_n[i] = self.a_1*self.y[i] + self.a_2*self.y_minus_n[i] + self.a_3*(self.y[i+1]+self.y[i-1]) \
        + self.a_4*(self.y[i+2]+self.y[i-2]) + self.a_5*(self.y_minus_n[i+1]+self.y_minus_n[i-1]+self.y_minus_2n[i]) \
        + (self.delta_t**2*self.N*F*g[i])/Ms
    