""" Simulation for a piano
Piano strings are simulated using finite differences. The wave equation
is extended with stiffness terms and damping terms.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.io.wavfile import write
import timeevolution as te
import csv

class String:

  """ All the available notes are defined here """
  def __init__(self,note,duration):
    self.note = note
    self.delta_t = 1/(4*44e3)
    self.duration = int(duration/self.delta_t)
    
    # Import note parameters from csv file
    with open('notes.csv', mode='r') as infile:
      reader = csv.reader(infile)
      notesDict = {rows[0]:rows[1:] for rows in reader}
      
    try:
      self.L = float(notesDict[note][0])
      self.Ms = float(notesDict[note][1])
      self.T = float(notesDict[note][2])
      self.b_1 = float(notesDict[note][3])
      self.b_2 = float(notesDict[note][4])
      self.epsilon = float(notesDict[note][5])
      self.N = float(notesDict[note][6])
    except KeyError:
      raise KeyError('There are no parameters found for the note '+note+'. Add them to notes.csv')

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
    
    zeta_b = 1000   # Bridge normalized impedance
    zeta_l = 1e20   # Left-end normalized impedance
    
    # Hammer
    F = 15.
    durationF = 352
    
    initialhammerheight = 0.005
    initialhammervelocity = 0.8
    self.hammerK = 4.5e9
    self.hammerP = 2.5
    self.hammerM = 2.97e-3
    hammer_length = round(0.3/self.delta_x)
    hammer_center_position = round(1/7*self.N)
    
    beginwindow = int(np.floor((hammer_center_position-hammer_length/2)))
    endwindow = int(np.floor(hammer_center_position+hammer_length/2))
    g = np.zeros((self.N),dtype = float)
    for i in range(beginwindow,endwindow+1):
        g[i] = (0.005/(0.03*0.03))*i*i*(float(endwindow)-float(beginwindow))*self.delta_x
    
    # Sound bridge
    bridgeposition = 0.8   # value between 0 and 1
    
    # String damper
    b_1Damped = self.b_1*damping
    b_2Damped = self.b_2*damping
    dampIndex = dampT/self.delta_t
    
    # Actual time evolution of the string
    self.time_evolved_string = te.timeevolution.time_evolution_bridge(self.y,g,self.N,self.duration,self.Ms,bridgeposition,self.delta_t,F,durationF,self.b_1,self.b_2,self.delta_x,self.r,self.mu,zeta_b,zeta_l,self.rho,b_1Damped,b_2Damped,dampIndex,initialhammerheight,initialhammervelocity,self.hammerK,self.hammerP,self.hammerM)

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
