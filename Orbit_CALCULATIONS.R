


#GM - reduced for display purposes
mu = 398600
# radius of central body
Re = 6371.01#kilometers


## Begin with initial orbit parameters
# input on dashboard
#periapse
P_o = 195.2#Kilometers
#apiapse
A_o = 250
#diamter of ellipse
D_o = P_o +A_o + 2*Re
# eccentricity
e0 = (Re + A_o - Re - P_o)/(2*Re +A_o +P_o)



## final orbit parameters
# input on dashboard
#periapse
P_f = 192.2
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



#########################################
#########################################
## Build the ellipse path - initial orbit

#semi-major axis
a = (P_o + A_o +2*Re )/2

# distance from center to focal point
c = a - P_o - Re

#semi-minor axis
b = sqrt(a^2 - c^2)

#Period
T_o = 2*pi*sqrt(a^3/mu)

# angular speed
w = 2*pi/T_o
# time-steps
#Number of orbits
N = 10
ts = 1## GLOBAL
t = seq(0 , N*T_o, ts)


# elliptical path of orbit(ideal ellipse)
## apply w*t later...
Ex = a*cos(w*t) +c
Ey = b*sin(w*t)


###############################################
#### Build Model Orbit with orbit equation
r0 = (P_o + Re)/(1 + e0*cos(w*t) )
d = A_o -P_o
n1 = (a+c-d)/ (  (P_o + Re)/(1 + e0*cos(w*0) )  )
n2 = n1
# model coords.
Ex1 = n1*r0*cos(w*t) +d
Ey1 = n2*r0*sin(w*t)


# radial velocity
r_dot = 0
for(i in 0:length(Ex1)){
  r_dot[i] = (r0[i+1]-r0[i])/ts
}#ENDOF 'for' R-dot
#plot(t,r_dot) - analytical solution r'
f_dot = (P_o+Re)*e0*w*sin(w*t)/(1+e0*cos(w*t))^2
#plot(t,f)
#plot(t,r_dot, lwd = 5)
#points(t,f_dot, col='red', lwd = 1)


# radial acceleration - assuming numerical is exact since r_dot is exact
r_ddot = 0
for(i in 0:length(Ex1)){
  r_ddot[i] = (r_dot[i+1]-r_dot[i])/ts
}#ENDOF 'for' R-dot
#plot(t,r_ddot)


### Generate Model Data Frames
Orbit_0 = data.frame(t,P_o,A_o,e0,a,b,c,H_Po,V_Po,H_Ao,V_Ao,w,Ex,Ey )
MODEL_0 = data.frame(t,P_o,A_o,e0,a,b,c,H_Po,V_Po,H_Ao,V_Ao,w,Ex1,Ey1,r_dot,r_ddot )



###############################
###############################
##### Build Final Orbit Ellipse

#semi-major axis
a1 = (P_f + A_f +2*Re )/2

# distance from center to focal point
c1 = a1 - P_f - Re

#semi-minor axis
b1 = sqrt(a1^2 - c1^2)

#Period
T_f = 2*pi*sqrt(a1^3/mu)

# angular speed
w1 = 2*pi/T_f
# time-span for final orbit
tf = seq(0 , N*T_f, ts)


Exf = a*cos(w1*tf) +c
Eyf = b*sin(w1*tf)


###############################################
#### Build FINAL Model Orbit with orbit equation
rf = (P_f + Re)/(1 + ef*cos(w1*tf) )
d1 = A_f -P_f
n3 = (a1+c1-d1)/ (  (P_f + Re)/(1 + ef*cos(w1*0) )  )
n4 = n3
# model coords.
Ex2 = n3*rf*cos(w1*tf) +d1
Ey2 = n4*rf*sin(w1*tf)

#plot(Ex2,Ey2)


#FINAL
# radial velocity
r_dot2 = 0
for(i in 0:length(Ex2)){
  r_dot2[i] = (rf[i+1]-rf[i])/ts
}#ENDOF 'for' R-dot
#plot(t,r_dot) - analytical solution r'
f_dot2 = (P_f+Re)*ef*w1*sin(w1*tf)/(1+ef*cos(w1*tf))^2
#plot(tf,f2)
#plot(tf,r_dot2, lwd = 5)
#points(tf,f_dot2, col='red', lwd = 1)


#FINAL
# radial acceleration - assuming numerical is exact since r_dot is exact
r_ddot2 = 0
for(i in 0:length(Ex2)){
  r_ddot2[i] = (r_dot2[i+1]-r_dot2[i])/ts
}#ENDOF 'for' R-dot
#plot(tf,r_ddot2)

# Final Orbit Data Frames
Orbit_f = data.frame(tf,P_f,A_f,ef,a1,b1,c1,H_Pf,V_Pf,H_Af,V_Af,w1,Exf,Eyf )
MODEL_f = data.frame(tf,P_f,A_f,ef,a1,b1,c1,H_Pf,V_Pf,H_Af,V_Af,w1,Ex2,Ey2,r_dot2,r_ddot2 )




