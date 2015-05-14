# piano wire simulation
from vibration import String
from vibration import Chord
import numpy as np

string1 = String('c4',3)
string1.time_evolution_f90()
note1 = string1.getWave()
# string1.exportSound()

string2 = String('e4',3)
string2.time_evolution_f90()
note2 = string2.getWave()
# string2.exportSound()

string3 = String('g4',3)
string3.time_evolution_f90()
note3 = string3.getWave()
# string2.exportSound()

chord = Chord()
chord.createChord(note1,note2,note3)
chord.exportSound('c-majeur')