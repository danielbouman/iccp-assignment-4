import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.io.wavfile import write
import timeevolution as te

class String:
  
  def __init__(self,note,duration):
    self.note = note
    self.duration = int(duration*4*44e3)
    self.delta_t = 1/(4*44e3)
    zeta_b = 1000
    zeta_l = 1e20

    # Choose paramters
    if str.lower(note) == 'c2':
      # Bass note
      self.L = 1.92
      self.Ms = 35e-3
      self.T = 750
      b_1 = 0.25
      b_2 = 7.5e-5
      self.epsilon = 7.5e-6
      self.N = 521
    if str.lower(note) == 'c4':
      # Midrange note
      self.L = 0.62
      self.Ms = 39.3e-3
      self.T = 670
      b_1 = 1.1
      b_2 = 2.7e-4
      self.epsilon = 3.82e-5
      self.N = 140
    if str.lower(note) == 'c7':
      # Treble note
      self.L = 0.09
      self.Ms = 0.467e-3
      self.T = 750
      b_1 = 9.17
      b_2 = 2.1e-3
      self.epsilon = 867e-4
      self.N = 23
    
    self.rho = self.Ms/self.L
    self.c = np.sqrt(self.T/self.rho)
    self.delta_x = self.L/self.N
    self.kappa_squared = self.epsilon*(self.c**2)*(self.L**2)
    mu = self.kappa_squared/((self.delta_x**2)*(self.c**2))

    D = 1+b_1*self.delta_t
    self.r = self.c*self.delta_t/self.delta_x
    if self.r> 1:
        raise ValueError("r is larger than 1, simulation is unstable. r should be smaller than 1.")
    if self.delta_t > (self.delta_x)/self.c:
        raise ValueError("CFL condition broken")
    
    nu = 2*b_2*(self.delta_t)/(self.delta_x**2)

    self.a_1 = (-(self.r**2)*mu)/D
    self.a_2 = (self.r**2+4*(self.r**2)*mu+nu)/D
    self.a_3 = (2-2*(self.r**2)-6*(self.r**2)*mu-2*nu)/D
    self.a_4 = (-1+b_1*self.delta_t+2*nu)/D
    self.a_5 = -nu/D

    b_R_denom = 1 + b_1*self.delta_t+zeta_b*self.r
    self.b_R1 = (2-2*self.r**2*mu - 2*self.r**2)/b_R_denom
    self.b_R2 = (4*self.r**2*mu - 2*self.r**2)/b_R_denom
    self.b_R3 = (-2*self.r**2*mu)/b_R_denom
    self.b_R4 = (-1+b_1*self.delta_t + zeta_b*self.r)/b_R_denom
    self.b_RF = (self.delta_t**2/self.rho)/b_R_denom
    
    b_L_denom = 1 + b_1*self.delta_t+zeta_l*self.r
    self.b_L1 = (2-2*self.r**2*mu - 2*self.r**2)/b_L_denom
    self.b_L2 = (4*self.r**2*mu - 2*self.r**2)/b_L_denom
    self.b_L3 = (-2*self.r**2*mu)/b_L_denom
    self.b_L4 = (-1+b_1*self.delta_t + zeta_l*self.r)/b_L_denom
    self.b_LF = (self.delta_t**2/self.rho)/b_L_denom

    # y_[plus/minus]_[n/2n] are string displacement vectors
    # self.y_plus_n = np.zeros((self.N),dtype = float)
    # self.y_minus_n = np.zeros((self.N),dtype = float)
    self.y = np.zeros((self.N),dtype = float)
    self.time_evolved_string = np.zeros((3,self.duration),dtype = float)

  def time_evolution_f90(self):
         F = 0.
         init_velocity = 1.0
         hammer_length = round(0.1/self.delta_x)
         hammer_center_position = round(1/8*self.N)
         hammer_displacement = 0.
         hammer_displacement_minus_n = -init_velocity*self.delta_t
         hammer_mass = 1.
         K = 1.
         p = 2.
         g = np.zeros((self.N),dtype = float)
         n = 0
         beginwindow = int(np.floor((hammer_center_position-hammer_length/2)))
         endwindow = int(np.floor(hammer_center_position+hammer_length/2))
         g[beginwindow:endwindow-1] = 1.0/(endwindow-beginwindow)
         bridgeposition = 0.8   # value between 0 and 1
         
         self.time_evolved_string = te.timeevolution.time_evolution_bridge(self.y,g,self.N,self.duration,self.Ms,bridgeposition,self.b_L1,self.b_L2,self.b_L3,self.b_L4,self.b_LF,self.a_1,self.a_2,self.a_3,self.a_4,self.a_5,self.b_R1,self.b_R2,self.b_R3,self.b_R4,self.b_RF,self.delta_t)
          
  def saveSound(self,extra_name=''):
    scaled_data = np.int16(self.time_evolved_string[0,:]/np.max(np.abs(self.time_evolved_string[0,:])) * 32767)
    write(self.note+'_string'+''+'.wav', int(4*44e3), scaled_data)
    if extra_name != '':
      write('string_'+extra_name+'.wav', int(4*44e3), scaled_data)
      
  def animate(self,saveAnimation=False):
    fig, ax = plt.subplots()
    self.xAxis = np.linspace(0,self.L,self.N)
    
    line, = ax.plot(self.xAxis, np.sin(self.xAxis))

    def animate(i):
        line.set_ydata(self.time_evolved_string[:,i])  # update the data
        return line,

    #Init only required for blitting to give a clean slate.
    def init():
        line.set_ydata(np.ma.array(self.xAxis, mask=True))
        return line,

    ani = animation.FuncAnimation(fig, animate, np.arange(1, self.duration), init_func=init,
        interval=5, blit=True,repeat=False)
    plt.axis([0,self.L,-0.015,0.015])
    if saveAnimation == True:
      # Set up formatting for the movie files
      Writer = animation.writers['ffmpeg']
      writer = Writer(fps=15, metadata=dict(artist='Bouman and Goodenough'))
      print("Saving animation...")
      ani.save('vibration.mp4', writer=writer)
      print("Done. Animation saved as vibration.mp4")
    else:
      plt.show()
    return