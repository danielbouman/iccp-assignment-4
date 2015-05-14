""" Simulation for a piano
Piano strings are simulated using finite differences. The wave equation
is extended with stiffness terms and damping terms.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.io.wavfile import write
import timeevolution as te

class String:

  """ All the available notes are defined here """
  def __init__(self,note,duration):
    self.note = note
    self.delta_t = 1/(4*44e3)
    self.duration = int(duration/self.delta_t)
    
    if str.lower(note) == 'c2':
      # Bass note
      self.L = 1.92
      self.Ms = 3.5e-3
      self.T = 750
      self.b_1 = 0.25
      self.b_2 = 7.5e-5
      self.epsilon = 7.5e-6
      self.N = 523
    if str.lower(note) == 'c4':
      # Midrange note
      self.L = 0.62
      self.Ms = 3.93e-3
      self.T = 650.8
      self.b_1 = 1.1
      self.b_2 = 2.7e-4
      self.epsilon = 3.82e-5
      self.N = 140
    if str.lower(note) == 'd4':
      # Midrange note
      self.L = 0.543
      self.Ms = 3.86e-3
      self.T = 703.5
      self.b_1 = 1.25
      self.b_2 = 2.7e-4
      self.epsilon = 3.82e-5
      self.N = 140
    if str.lower(note) == 'e4':
      # Midrange note
      self.L = 0.466
      self.Ms = 3.81e-3
      self.T = 750.8
      self.b_1 = 1.5
      self.b_2 = 4.7e-4
      self.epsilon = 4.12e-5
      self.N = 100
    if str.lower(note) == 'f4':
      # Midrange note
      self.L = 0.458
      self.Ms = 3.51e-3
      self.T = 756.8
      self.b_1 = 1.7
      self.b_2 = 4.9e-4
      self.epsilon = 4.29e-5
      self.N = 70
    if str.lower(note) == 'g4':
      # Midrange note
      self.L = 0.406
      self.Ms = 3.31e-3
      self.T = 801.8
      self.b_1 = 1.8
      self.b_2 = 5.7e-4
      self.epsilon = 4.92e-5
      self.N = 90
    if str.lower(note) == 'a4':
      # Midrange note
      self.L = 0.3576
      self.Ms = 3.062e-3
      self.T = 809.9
      self.b_1 = 1.3
      self.b_2 = 2.5e-4
      self.epsilon = 5.12e-5
      self.N = 49
    if str.lower(note) == 'b4':
      # Midrange note
      self.L = 0.2876
      self.Ms = 3.068e-3
      self.T = 819.0
      self.b_1 = 1.3
      self.b_2 = 2.5e-4
      self.epsilon = 5.12e-5
      self.N = 45
    if str.lower(note) == 'g4_2':
      # Midrange note
      self.L = 0.522
      self.Ms = 2.41e-3
      self.T = 735.8
      self.b_1 = 1.8
      self.b_2 = 5.7e-4
      self.epsilon = 4.92e-5
      self.N = 40
    if str.lower(note) == 'g4_3':
      # Midrange note
      self.L = 0.522
      self.Ms = 2.41e-3
      self.T = 735.8
      self.b_1 = 3.8
      self.b_2 = 11.7e-4
      self.epsilon = 14.92e-5
      self.N = 40
    if str.lower(note) == 'g4_4':
      # Midrange note
      self.L = 0.522
      self.Ms = 2.41e-3
      self.T = 735.8
      self.b_1 = 2.8
      self.b_2 = 8.2e-4
      self.epsilon = 9.92e-5
      self.N = 40
    if str.lower(note) == 'g4_5':
      # Midrange note
      self.L = 0.522
      self.Ms = 2.41e-3
      self.T = 735.8
      self.b_1 = 0.8
      self.b_2 = 2.2e-4
      self.epsilon = 1.92e-5
      self.N = 40
    if str.lower(note) == 'g4_6':
      # Midrange note
      self.L = 0.522
      self.Ms = 2.41e-3
      self.T = 735.8
      self.b_1 = 0.4
      self.b_2 = 1.1e-4
      self.epsilon = 0.42e-5
      self.N = 40
    if str.lower(note) == 'g4_7':
      # Midrange note
      self.L = 0.522
      self.Ms = 2.41e-3
      self.T = 735.8
      self.b_1 = 0.2
      self.b_2 = 0.1e-4
      self.epsilon = 2.42e-4
      self.N = 40
    if str.lower(note) == 'g4_8':
      # Midrange note
      self.L = 0.406
      self.Ms = 3.31e-3
      self.T = 801.8
      self.b_1 = 1.8
      self.b_2 = 5.7e-4
      self.epsilon = 8.92e-4
      self.N = 49
    if str.lower(note) == 'g4_9':
      # Midrange note
      self.L = 0.406
      self.Ms = 3.31e-3
      self.T = 801.8
      self.b_1 = 0.8
      self.b_2 = 5.7e-4
      self.epsilon = 4.92e-5
      self.N = 100
    if str.lower(note) == 'g4_10':
      # Midrange note
      self.L = 0.406
      self.Ms = 3.31e-3
      self.T = 801.8
      self.b_1 = 1.8
      self.b_2 = 35.7e-4
      self.epsilon = 4.92e-5
      self.N = 100
    if str.lower(note) == 'g5':
      # Treble note
      self.L = 0.165
      self.Ms = 2.53e-3
      self.T = 973.8
      self.b_1 = 1.8
      self.b_2 = 35.7e-4
      self.epsilon = 4.92e-5
      self.N = 55
    if str.lower(note) == 'c7':
      # Treble note
      self.L = 0.09
      self.Ms = 0.467e-3
      self.T = 750
      self.b_1 = 9.17
      self.b_2 = 2.1e-3
      self.epsilon = 867e-4
      self.N = 5
    if str.lower(note) == 'g3':
      # Midrange note
      self.L = 0.925
      self.Ms = 4.9e-3
      self.T = 670.8
      self.b_1 = 1.8
      self.b_2 = 5.7e-4
      self.epsilon = 4.92e-5
      self.N = 150
    
    # Define general parameters
    self.rho = self.Ms/self.L
    self.c = np.sqrt(self.T/self.rho)
    self.delta_x = self.L/self.N
    self.kappa_squared = self.epsilon*(self.c**2)*(self.L**2)
    self.mu = self.kappa_squared/((self.delta_x**2)*(self.c**2))
    self.r = self.c*self.delta_t/self.delta_x
    
    # Stability checks
    if self.r > 1:
        raise ValueError("Courant–Friedrichs–Lewy (CFL) condition is not met. \n Choose smaller time step or larger time steps.")
    if self.delta_t > self.delta_x/self.c:
        raise ValueError("Upper bound for the time step reached.\n Choose larger time steps or smaller spacial steps.")
    
    # Preallocation of string deflection arrays
    self.y = np.zeros((self.N),dtype = float)
    self.time_evolved_string = np.zeros((3,self.duration),dtype = float)
  
  """ Time evolution of the string. The discretized differential equation is calculated in a Fortran module. """
  def time_evolution_f90(self,damping=1,dampT=0):
    
    zeta_b = 1000   # Left-end normalized impedance
    zeta_l = 1e20   # Bridge normalized impedance
    
    # Hammer
    F = 15.
    durationF = 90
    
    hammer_length = round(0.1/self.delta_x)
    hammer_center_position = round(1/8*self.N)
    
    beginwindow = int(np.floor((hammer_center_position-hammer_length/2)))
    endwindow = int(np.floor(hammer_center_position+hammer_length/2))
    g = np.zeros((self.N),dtype = float)
    g[beginwindow:endwindow-1] = 1.0/(endwindow-beginwindow)
    
    # Sound bridge
    bridgeposition = 0.8   # value between 0 and 1
    
    # String damper
    b_1Damped = self.b_1*damping
    b_2Damped = self.b_2*damping
    dampIndex = dampT/self.delta_t
    
    # Actual time evolution of the string
    self.time_evolved_string = te.timeevolution.time_evolution_bridge(self.y,g,self.N,self.duration,self.Ms,bridgeposition,self.delta_t,F,durationF,self.b_1,self.b_2,self.delta_x,self.r,self.mu,zeta_b,zeta_l,self.rho,b_1Damped,b_2Damped,dampIndex)

  """ Export note to a wave data for creating a chord. """
  def getWave(self):
    return self.time_evolved_string[0,:]
  
  """ Export note to a .wav file. """
  def exportSound(self,extra_name=''):
    scaled_data = np.int16(self.time_evolved_string[0,:]/np.max(np.abs(self.time_evolved_string[0,:])) * 32767)
    if extra_name != '':
      print('Saved as: '+self.note+'_string_'+extra_name+'.wav')
      write(self.note+'_string_'+extra_name+'.wav', int(4*44e3), scaled_data)
    else:
      print(self.note+'_string'+''+'.wav')
      write(self.note+'_string'+''+'.wav', int(4*44e3), scaled_data)

  """ Plot an animation of the string deflection. """
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
  
class Chord:
  def createChord(self,*args):
    if len(args) == 1:
      self.notesCombined = args[0]
    if len(args) == 2:
      self.notesCombined = np.add(args[0],args[1])
    if len(args) == 3:
      self.notesCombined = np.add(args[0],np.add(args[1],args[2]))
      
  def exportSound(self,name=''):
    scaled_chord = np.int16(self.notesCombined/np.max(np.abs(self.notesCombined)) * 32767)
    write('chord_'+name+'.wav', int(4*44e3), scaled_chord)