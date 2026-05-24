

#GM - reduced for display purposes
mu = 398600
# radius of central body
Re = 6371.01#kilometers


## Begin with initial orbit parameters
# input on dashboard
#periapse
P_o = 195.4#Kilometers
#apiapse
A_o = 264.9
#diamter of ellipse
D_o = P_o +A_o + 2*Re
# eccentricity
e0 = (Re + A_o - Re - P_o)/(2*Re +A_o +P_o)



## final orbit parameters
# input on dashboard
#periapse
P_f = 195.4
#apiapse
A_f = 300
#diamter of final ellipse
D_f = P_f +A_f + 2*Re
# eccentricity
ef = (Re + A_f - Re - P_f)/(2*Re +A_f +P_f)


# Momentum and Velocity/DV calculations Periapse
H_Po = sqrt((Re+P_o)*(1+e0)*mu) 
H_Pf = sqrt((Re+P_f)*(1+ef)*mu) ## momentum at P-o for mid-transfer elliptic orbit
V_Po = H_Po/(Re+P_o)
V_Pf = H_Pf/(Re+P_f)### velocity at periapse_o! for mid-hohmman transfer.
DV_Pof = V_Pf - V_Po

# Momentum and Velocity/DV calculations Apiapse
H_Ao = sqrt((Re+A_o)*(1+e0)*mu) 
H_Af = sqrt((Re+A_f)*(1+ef)*mu) 
V_Ao = H_Ao/(Re+ A_o)
V_Af = H_Af/(Re+ A_f)
DV_Aof = V_Af - V_Ao


Thrust = 6197#N/s
M_dot = 2.1#kg/s
#initial Mass of spacecraft
M_o = 7259.72 #Kg

# from periapse to transition path
BT10 = M_o / ( Thrust/(DV_Pof*1000) +2.1  )
BT10


BT1 = M_o/( ( (Thrust-  (M_o-2.1*BT10)*(DV_Pof/BT10)) / (DV_Pof*1000) )  +2.1  ) 
BT1


## Begin with circ. orbit parameters
#periapse
P_o = 195.4#Kilometers
#apiapse
A_o = 300
#diamter of ellipse
D_o = P_o +A_o + 2*Re
# eccentricity
e0 = (Re + A_o - Re - P_o)/(2*Re +A_o +P_o)



## final orbit parameters
# input on dashboard
#periapse
P_f = 300
#apiapse
A_f = 300
#diamter of final ellipse
D_f = P_f +A_f + 2*Re
# eccentricity
ef = (Re + A_f - Re - P_f)/(2*Re +A_f +P_f)


# Momentum and Velocity/DV calculations Periapse
H_Po = sqrt((Re+P_o)*(1+e0)*mu) 
H_Pf = sqrt((Re+P_f)*(1+ef)*mu) ## momentum at P-o for mid-transfer elliptic orbit
V_Po = H_Po/(Re+P_o)
V_Pf = H_Pf/(Re+P_f)### velocity at periapse_o! for mid-hohmman transfer.
DV_Pof = V_Po - V_Pf## reversed for circuliarization

# Momentum and Velocity/DV calculations Apiapse
H_Ao = sqrt((Re+A_o)*(1+e0)*mu) 
H_Af = sqrt((Re+A_f)*(1+ef)*mu) 
V_Ao = H_Ao/(Re+ A_o)
V_Af = sqrt(mu/(Re+A_f))
DV_Aof = V_Ao - V_Af



M_1 = 7232.85#M_o - M_dot*BT1
# from apiapse to circularization
BT20 = M_1 / ( Thrust/((DV_Aof)*1000) +2.1  )
BT20

BT2 = M_1/ ( ( (Thrust-  ((M_1)-2.1*BT20)*(DV_Aof/BT20)) / (DV_Aof*1000) )  +2.1  ) 
BT2## correction


## mass after circularization
M_2 = M_1 - M_dot*BT2


# time = 0
# 
# V_t = (M_o*Thrust*time  )/( M_o - M_dot*time  )
# 
# V_select = 1350# input
# 
# e_n = (Re+P_o)*V_select^2/mu -1
# 
# A_n = ( e_n*P_o + 2*e_n*Re+P_o )/( 1 - e_n )
# 
# #build new ellipse for display
# a_n = (P_o + A_n +2*Re)/2
# 
# c_n = a_n - P_o - Re
# 
# b_n = sqrt(a_n^2 - c_n^2  )
# 
# T_n = 2*pi*sqrt(a_n^3/mu)
# 
# w_n = 2*pi/T_n
# 
# t_n = seq(0 , T_n, ts)
# 
# #plot( a_n*cos(w_n*t_n)+c_n, b_n*sin(w_n*t_n)       )
# #points(Ex1,Ey1, col = 'blue')
# 
# 
# 







