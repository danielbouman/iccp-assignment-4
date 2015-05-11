import numpy as np

class String:
  
  def __init__(self,N,L,b_1,b_3):
    self.N = N
    self.D = 1+b_1*delta_t+(2*b_3/delta_t)
    self.delta_x = N/L
    self.r = c*delta_t/delta_x
    self.delta_t = 0.1
    self.epsilon = 1
    # create constants 

    self.a_1 = (2-2*r**2+(b_3/delta_t)-6*epsilon*N**2*r**2)/delta_t
    self.a_2 = (-1+b_1*delta_t+2*b_3/delta_t)/D
    self.a_3 = ((r**2)*(1+4*epsilon*N**2))/D
    self.a_4 = (b_3/delta_t-epsilon*N**2*r**2)/D
    self.a_5 = (-b_3/(delta_t*D))

    # y_[plus/minus]_[n/2n] are string displacement vectors
    self.y_plus_n = np.zeros((N),dtype = float)
    self.y = np.zeros((N),dtype = float)
    self.y_minus_n = np.zeros((N),dtype = float)
    self.y_minus_2n = np.zeros((N),dtype = float)


  def displacement(self):
    

  def time_evolution(self):
    F = 0
    init_velocity = 1.0
    hammer_length = 0.01
    hammer_position = 0.15*L
    hammer_displacement = 0.0
    hammer_displacement_minus_n = -init_velocity*delta_t
    hammer_mass = 1.0
    K = 1.0
    p = 2.0
    Ms = 1.0
    g = np.zeros((N),dtype = float)
    n = 0
    beginwindow = np.floor((hammer_position-hammer_length/2)*N/L)
    endwindow = np.ceil((hammer_position+hammer_length/2)*N/L)    
    
    for i in range(beginwindow,endwindow):
        g[i] = 1.0/(endwindow-beginwindow)
    
    F = K/hammer_mass * np.abs(hammer_displacement - self.y)**p    
    for i in range(2,N-2):
        self.y_plus_n[i] = self.a_1*self.y[i] + self.a_2*self.y_minus_n[i] + self.a_3*(self.y[i+1]+self.y[i-1])]
        + self.a_4*(self.y[i+2]+self.y[i-2]) + self.a_5*(self.y_minus_n[i+1]+self.y_minus_n[i-1]+self.y_minus_2n[i])
        + (self.delta_t**2*self.N*F*g)/Ms
    hammer_displacement = 2*hammer_displacement - hammer_displacement_minus_n
    + F*delta_t**2
    

