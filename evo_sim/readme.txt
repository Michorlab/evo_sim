main_rf_bp is program to run gradual branching process simulation
arguments are for all scenarios and currently the only way to
run specific scenarios is to zero out particular variables

To Do:
Create a wrapper or method to require user to only 
include parameters of interest: either make input a dictionary
type (python wrapper), overload functions, or define some file type
(JSON) which reads parameters.

Option for sampling as well as how many of each type to sample

Include punctuated model into this as well - maybe make an argument for
the following scenarios:

gradual or punctuated
baseline or mutation distribution
none or epistasis 1 (x mutations then burst) or epistasis 2
constant mutation rate or time-burst
combinations (2x2x3x2 = 24 possible models, but some unecessary)

