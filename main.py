# piano wire simulation
from vibration import String

# Variables
N = 1000
L = 1.9
b_1 = 6.25e-9
b_3 = 0.5

string1 = String('c2',5)
string1.time_evolution_f90()
# string1.animate(saveAnimation=False)
# string1.saveSound('c2test',amp=10000)