GlowScript 3.1 VPython


G = 6.67430e-11
ME = 5.972e24
#GM
mu = 398600#km^3/s^2
Re = 6371.01#km


## Begin with initial orbit parameters
# input on dashboard
#periapse
P_o = 60*Re#Kilometers
#apiapse
A_o = 60*Re
#diamter of ellipse
D_o = P_o +A_o + 2*Re
# eccentricity
e0 = (Re + A_o - Re - P_o)/(2*Re +A_o +P_o)


## final orbit parameters
# input on dashboard
#periapse
P_f = P_o
#apiapse
A_f = A_o
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


#semi-major axis
a = (P_o + A_o +2*Re )/2

# distance from center to focal point
c = a - P_o - Re

#semi-minor axis
b = sqrt(a**2 - c**2)

#Period
T_o = 2*pi*sqrt(a**3/mu)

# angular speed
w = 2*pi/T_o
# time-steps
#Number of orbits
N = 0.05
ts = 1## GLOBAL
#t = range(0,T_o,ts)


# elliptical path of orbit(ideal ellipse)
## apply w*t later...
#Ex = a*cos(w*t) +c
#Ey = b*sin(w*t)

Ex1 = []
Ey1 = []
Ez1 = []

t = range(0,N*T_o,ts)
print(t[0])
print(len(t))

for i in t:
  r0 = (P_o + Re)/(1 + e0*cos(w*(t[i]-1)) )
  d = A_o -P_o
  n1 = (a+c-d)/ (  (P_o + Re)/(1 + e0*cos(w*0) )  )
  n2 = n1
  # model coords.
  Ex1.append(n1*r0*cos(w*t[i]) +d)
  Ey1.append( n2*r0*sin(w*t[i]))
  Ez1.append(0)

print("length of Ex1 = ",len(Ex1))
print("T_o = " , T_o)

#lineLAT = curve( pos = vec(Re*cos(range(0,2*pi,0.1)),0,Re*sin(range(0,2*pi,0.1))) , color = color.purple)
earth = sphere(pos=vector(0,0,0),radius=Re, texture=textures.earth,make_trail=True)
craft = sphere(radius=1, color=color.yellow,make_trail=True)
shng =  sphere(radius=10, color=color.yellow,texture = textures.flower)
shng.pos = vector(Re*sin(115*pi/180)*cos(37*pi/180),Re*sin(115*pi/180)*sin(37*pi/180),Re*cos(115*pi/180) )
#moon =  sphere(radius=Re/6, color=color.white,make_trail=True)
g1 = graph(xtitle="LAT", ytitle="LON", width=500, height=250)
Cgraph = gcurve(color=color.blue, label = "TR", dot=True)

g2 = graph(xtitle="t", ytitle="theta_dot", width=500, height=250)
Vgraph = gcurve(color=color.purple, label = "TR", dot=True)

g3 = graph(xtitle="t", ytitle="psi_dot", width=500, height=250)
Agraph = gcurve(color=color.red, label = "TR", dot=True)

g4 = graph(xtitle="t", ytitle="Alt.", width=500, height=250)
Rgraph = gcurve(color=color.red, label = "TR", dot=True)



#lineLAT.pos = Re*vector(cos(range(0,2*pi,0.1)),0,sin(range(0,2*pi,0.1)))       

# initial Conditions

M_dot = 2.1#kg/s
#initial Mass of spacecraft
M_o = 7280 #Kg


#initial position/motion of craft(3rd body)
bR0 = Re
LON = -150*pi/180
LAT = -44*pi/180# HI ADIZ @ LON -150, LAT -44

bR0x = bR0*sin(LON)*cos(LAT)
bR0y = bR0*sin(LON)*sin(LAT)
bR0z = bR0*cos(LON)
craft.pos = vector(bR0x,bR0y,bR0z)

Thrust = 110# kN
bRv0 = 0
bRa0 = Thrust/M_o

AV_lat = LAT 
AV_lon = LON 

bRv0x = bRv0*sin(AV_lon)*cos(AV_lat)
bRv0y = bRv0*sin(AV_lon)*sin(AV_lat)
bRv0z = bRv0*cos(AV_lon)

bRa0x = bRa0*sin(AV_lon)*cos(AV_lat)
bRa0y = bRa0*sin(AV_lon)*sin(AV_lat)
bRa0z = bRa0*cos(AV_lon)

## acceleration wrt central body
AEx = [-(mu/(bR0**3))*bR0x]
AEy = [-(mu/(bR0**3))*bR0y]
AEz = [-(mu/(bR0**3))*bR0z]
#magnitude of acceleration
AE = [sqrt( (AEx[0]**2) + (AEy[0]**2) )]

MAE0 = AE*sqrt(AEx[0]**2 + AEy[0]**2 + AEz[0]**2)



print("MAE0 = ",MAE0)

## accel wrt moon
m2b = 0.012*ME# not used
# radius from moon body
mmR = [sqrt(  (Ex1[0] - bR0x)**2 + (Ey1[0] - bR0y)**2 +(Ez1[0] - bR0z)**2 )]
AMR = [(0.0123031469*mu)/(mmR[0]**3)]

MAMR0 = AMR[0]*mmR[0]

print("MAMR0 =", MAMR0 )


# moon component acceleration
AMRx = [AMR[0]*(bR0x-Ex1[0]) ]
AMRy = [AMR[0]*(bR0y-Ey1[0]) ]
AMRz = [AMR[0]*(bR0z-Ez1[0]) ]

mlist = [M_o]

Xlist = [bR0x]
Ylist = [bR0y]
Zlist = [bR0z]

print(Xlist)
print(Ylist)
print(Zlist)

Vxlist = [bRv0x]
Vylist = [bRv0y]
Vzlist = [bRv0z]

Axlist = [bRa0x]
Aylist = [bRa0y]
Azlist = [bRa0z]

theta = [LAT]
psi = [LON]

imp_count = 0
rng = 0

i=1
while i < 1510:
  rate(100)
###############################################
#### Build Model Orbit with orbit equation
  
  
  #moon.pos = vector(Ex1[i-1], Ey1[i-1],0)
  
  # time of fuel expiration
  if i > 591.2 : 
    Thrust = 0 #half the mass fo craft is fuel m - mdot*t = m/2
    
    
  if i > 120 : 
   #92.772, -336.168, thrust = 77, burntime = 77
    AV_lon =  115*pi/180
    AV_lat =  37*pi/180
  
  #x = Xlist[i-1] + 0.5*( AEx[i-1] + AMRx[i-1] + Axlist[i-1] + Thrust/(mlist[i-1] - M_dot*ts)*sin(AV_lon)*cos(AV_lat) )*(ts**2) + Vxlist[i-1]*ts + Thrust*ts*ts*sin(AV_lon)*cos(AV_lat)/(mlist[i-1]-M_dot*ts) 
  #y = Ylist[i-1] + 0.5*( AEy[i-1] + AMRy[i-1] + Aylist[i-1] + Thrust/(mlist[i-1] - M_dot*ts)*sin(AV_lon)*sin(AV_lat) )*(ts**2) + Vylist[i-1]*ts + Thrust*ts*ts*sin(AV_lon)*sin(AV_lat)/(mlist[i-1]-M_dot*ts)
  #z = Zlist[i-1] + 0.5*( AEz[i-1] + AMRz[i-1] + Azlist[i-1] + Thrust/(mlist[i-1] - M_dot*ts)*cos(AV_lon)             )*(ts**2) + Vzlist[i-1]*ts + Thrust*ts*ts*cos(AV_lon)            /(mlist[i-1]-M_dot*ts)
  
  x = Xlist[i-1] + 0.5*( AEx[i-1] + AMRx[i-1] + Axlist[i-1] + Thrust/(mlist[i-1] - M_dot*ts)*sin(AV_lon)*cos(AV_lat) )*(ts**2) + Vxlist[i-1]*ts 
  y = Ylist[i-1] + 0.5*( AEy[i-1] + AMRy[i-1] + Aylist[i-1] + Thrust/(mlist[i-1] - M_dot*ts)*sin(AV_lon)*sin(AV_lat) )*(ts**2) + Vylist[i-1]*ts 
  z = Zlist[i-1] + 0.5*( AEz[i-1] + AMRz[i-1] + Azlist[i-1] + Thrust/(mlist[i-1] - M_dot*ts)*cos(AV_lon)             )*(ts**2) + Vzlist[i-1]*ts 
  
  
  
  Xlist.append( x )
  Ylist.append( y )
  Zlist.append( z )
  craft.pos =  vector(Xlist[i-1],Ylist[i-1],Zlist[i-1])
  
  R = sqrt(Xlist[i]**2 +Ylist[i]**2+Zlist[i]**2)
  
  mlist.append(mlist[i-1]-M_dot*ts)
  
  AEx.append(-(mu/R**3)*Xlist[i])
  AEy.append(-(mu/R**3)*Ylist[i])
  AEz.append(-(mu/R**3)*Zlist[i])
  
  mmR.append(sqrt(  (Ex1[i] - Xlist[i])**2 + (Ey1[i] - Ylist[i])**2  +(Ez1[i] - Zlist[i])**2  ))
  AMR.append((0.0123031469*mu)/(mmR[i]**3))
  AMRx.append(AMR[i]*(Ex1[i]-Xlist[i]))
  AMRy.append(AMR[i]*(Ey1[i]-Ylist[i]))
  AMRz.append(AMR[i]*(Ez1[i]-Zlist[i]))
  
  
  Vxlist.append( Vxlist[i-1] + Axlist[i-1]*ts)
  Vylist.append( Vylist[i-1] + Aylist[i-1]*ts)
  Vzlist.append( Vzlist[i-1] + Azlist[i-1]*ts)
  
  Axlist.append(AEx[i-1] + AMRx[i-1]  + Thrust/(mlist[i-1] - M_dot*ts)*sin(AV_lon)*cos(AV_lat))
  Aylist.append(AEy[i-1] + AMRy[i-1]  + Thrust/(mlist[i-1] - M_dot*ts)*sin(AV_lon)*sin(AV_lat))
  Azlist.append(AEz[i-1] + AMRz[i-1]  + Thrust/(mlist[i-1] - M_dot*ts)*cos(AV_lon))
  
  
  theta.append(acos(Zlist[i-1]/(R))*180/pi)#LAT
  psi.append(acos(Xlist[i-1]/(R*sin(acos(Zlist[i-1]/(R)))))*180/pi)#LON
  
  rng = rng + sqrt((Xlist[i]-Xlist[i-1])**2+(Ylist[i]-Ylist[i-1])**2+(Zlist[i]-Zlist[i-1])**2)
  
  if i > floor(M_o/M_dot*ts/10) : 
    Rgraph.label = "NO FUEL! : impact count_" + imp_count + ": Rng: " + round(rng,2)
    
  else:
    Rgraph.label = "TR : impact count_" + imp_count + ": Rng: " + round(rng,2)
    
  if R < Re:
    Rgraph.label = "IMPACT! : impact count_" + imp_count + ": Rng: " + round(rng,2)
    imp_count = imp_count +1
    shng.radius += 3.5
  elif R <Re & (i > floor(M_o/M_dot*ts/10)):
    Rgraph.label = "TR: impact count_" + imp_count+ ": Rng: " + round(rng,2)
    shng.radius += 3.5
  Rgraph.plot(i,R-Re)  
  Cgraph.plot(theta[i-1],psi[i-1])
  Vgraph.plot(i,theta[i]-theta[i-1])
  Agraph.plot(i,psi[i]-psi[i-1])
  
  i = i+1
 
print("END PROGRAM")



