import os,sys
import numpy as np
# to run:
# python prog.py d50(m) particle density
# python3 calc_settlingV.py 2e-6 2700 

# definitions
ps = int(sys.argv[2]) # grain density km.m-3

p = 1027 # water density kg.m-3

v = 1.36e-6 # kinematic viscosity m2.s-1
print('v: ', v)

s = ps/p # ratio of grain density to water density
print('s: ',s)

d = float(sys.argv[1]) # d50 grain diameter (normally given in um), so ensure entered into the model in m
print(d)

# calculate D*
var1 = 9.81*(s-1)
var2 = v**2
var3 = var1/var2
var4 = var3**(1/3)

D = var4*d # D*

# three equations for Van Rijn:

if D**3 <= 16.187:
    ws = (v*(D**3))/(18*d)

elif D**3 > 16.187 and D**3 <= 16187:
    ws = ((10*v)/d)*((((1+0.01*(D**3)))**(1/2))-1.0)

elif D**3 > 16187:
    ws = (1.1*v*(D**(1.5)))/d

print('Settling Velocity: ', ws)
