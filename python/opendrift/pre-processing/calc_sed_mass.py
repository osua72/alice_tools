import os,sys
import numpy as np

# to run:
# python prog.py [dredge/disp: 1 = dredge, 2 = disp] [dredge type: 1 TSHD, 2, TSHD OF, 3 BHD] [nb_parts] volume2remove[site] [#of sites] [sed: 1 - sand, 2: silt, 3: clay] [sed_frac %]
# python3 calc_sed_mass.py 1 189400 3 0.44 1

# definitions

dredgeordisp = int(sys.argv[1]) #1 = dredge, 2 = disp

dredge = int(sys.argv[2]) # 1 = TSHD; 2 = BHD

nb_parts = int(sys.argv[3]) # of particles

vol2rem = int(sys.argv[4]) #m3

sites = int(sys.argv[5]) # number of sites to dredge per area

sed = int(sys.argv[6])

sed_frac = float(sys.argv[7]) #%


if dredge == 1:
    dredge_rate = 117600 #m3/week
    dredge_capacity = 635
    cycle_time = 3.5 # h
    loading_time = (cycle_time - 1.0)*3600 # loading time in seconds
elif dredge == 2:
    dredge_rate = 50000 # m3/week
    dredge_capacity = 1860
    cycle_time = 24 # h
    loading_time = (cycle_time - 1.0)*3600 # loading time in seconds

if dredgeordisp == 1:
    vol2rem_site = vol2rem/sites
    print('Volume to remove per site (m3)', vol2rem_site)

    loading_time = (cycle_time - 1.0)*3600 # loading time in seconds
    print('Loading Time (s): ', loading_time)

    placement_duration = 10*60
    print('Placement Duration (s): ', placement_duration)
    
    if sed == 1:
        dry_weight = 1600 # kg/m3
    else: dry_weight = 500 # kg/m3
    print('Dry Weight (kg/m3)', dry_weight)

    weeks = vol2rem_site/dredge_rate
    print('Days per site: ', weeks * 7)

    cycles = (weeks*(7*24*3600))/loading_time
    print('# Cycles', cycles)

    total_loading_time = cycles*loading_time
    print('Total loading time: ', total_loading_time)

    Drate = dredge_rate*(7*24*3600) # m3/s
    Prate = vol2rem_site/total_loading_time # m3/s

    vol_per_cycle = np.round(loading_time*Prate) # m3
    print('Volume to dispose per cycle (m3)', vol_per_cycle) 

elif dredgeordisp == 2:
    if sed == 1:
        dry_weight = 1600 # kg/m3
    else: dry_weight = 500 # kg/m3
    print('Dry Weight (kg/m3)', dry_weight)
    #    
    vol_per_cycle = dredge_capacity
    print('Volume to dispose per cycle (m3)', vol_per_cycle)    

mass_per_cycle = vol_per_cycle*dry_weight*sed_frac # kg 
print('mass of sand removed per cycle: ',mass_per_cycle) # kg

if dredgeordisp == 1:
    if dredge == 2:
        pdh0 = 3/100 #
        ss_mass_1 = mass_per_cycle*pdh0
        mass_hopper_1 = mass_per_cycle-ss_mass_1
        R0 = (2.5*60-2*60)/(2.5*60)
        print(R0)
        fsett = 0.25
        ftrap = 0.05
        mass_OF = R0*(1-fsett)*(1-ftrap)*mass_hopper_1
        pp0 = 20/100
        ss_mass_2 = mass_OF*pp0 # OF plume kg
        ss_mass_3 = (1-pp0)*mass_OF # density current kg
        mass_sum = ss_mass_1+ss_mass_2+ss_mass_3
    if dredge == 1:
        pdh0 = 3/100 #
        ss_mass_1 = mass_per_cycle*pdh0
        mass_sum = ss_mass_1

    if dredge == 3:
        pdh0 = 4/100 #
        ss_mass_1 = mass_per_cycle*pdh0 # bucket drip fraction
        mass_sum = ss_mass_1
        

if dredgeordisp == 2:
    if dredge == 1:
        pdh0 = 10/100 #%
        ss_mass_1 = mass_per_cycle*pdh0

    if dredge == 2:
        pdh0 = 5/100 #%
        ss_mass_1 = mass_per_cycle*pdh0 
    
    mass_sum = ss_mass_1

print('Total mass released per cycle: ', mass_sum) 

mass_per_part = mass_sum/nb_parts
print('Mass per particle (kg): ', mass_per_part) 







