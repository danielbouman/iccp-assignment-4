# piano wire simulation
from vibration import String

# Variables
N = 100
L = 1.9
b_1 = 6.25e-9
b_3 = 0.5

string1 = String(N,L,b_1,b_3)

string1.time_evolution()