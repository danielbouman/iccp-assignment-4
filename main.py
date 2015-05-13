# piano wire simulation
from vibration import String

string1 = String('c7',1)
string1.time_evolution_f90()
# string1.animate(saveAnimation=False)
string1.saveSound()