""" Electric piano simulation
"""
from vibration import String
from vibration import Chord
import numpy as np

string1 = String('c4',2)            # Define string with: note, duration in seconds
string1.time_evolution_f90()        # Start time evolution with: damper (optional), damper delay
# note1 = string1.getWave()         # Get wave when creating a chord
string1.exportSound('new_hammer2')     # Save as .wav file with: extra text in filename (optional)

# chord = Chord()
# chord.createChord(note1,note2,note3)
# chord.exportSound('c-majeur')