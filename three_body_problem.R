



G = 6.67430*10^-11
ME = 5.972*10^24

# initial Conditions
Thrust = 000#N
M_dot = 2.1#kg/s
#initial Mass of spacecraft
M_o = 7272 #Kg

#initial position of craft
bR0theta = 0
bR0 = 150000+Re
bR0x = bR0*cos(bR0theta)
bR0y = bR0*sin(bR0theta)


bRv0 = 0
bRa0 = Thrust/M_o
AV_theta = 90*pi/180

bRv0x = bRv0*cos(AV_theta)
bRv0y = bRv0*sin(AV_theta)
bRa0x = bRa0*cos(AV_theta)
bRa0y = bRa0*sin(AV_theta)


## acceleration wrt central body
AEx = -(mu/bR0^3)*bR0x
AEy = -(mu/bR0^3)*bR0y
#magnitude of acceleration
AE = sqrt(  AEx^2 + AEy^2 )


## accel wrt moon
m2b = 0.012*ME
# radius of moon body
mmR = sqrt(  (Ex1[1] - bR0x)^2 + (Ey1[1] - bR0y)^2  )
AMR = (mu/6)/mmR^3


# moon component acceleration
AMRx = AMR*(Ex1[1] - bR0x)
AMRy = AMR*(Ey1[1] - bR0y)

mlist = M_o

Xlist = bR0x
Ylist = bR0y

Vxlist = bRv0x
Vylist = bRv0y

Axlist = bRa0x
Aylist = bRa0y

#floor(M_o/M_dot)
for (i in 2:length(Ex1)){

# time of fuel expiration
if(i > floor(M_o/M_dot) ){Thrust = 0}

Xlist[i] = Xlist[i-1] + 0.5*( AEx[i-1] + AMRx[i-1] + Axlist[i-1] + Thrust/(mlist[i-1] - M_dot*ts  ))*ts^2 + (Vxlist[i-1] + Thrust*ts*cos(AV_theta))/(mlist[i-1]-M_dot*ts) 
Ylist[i] = Ylist[i-1] + 0.5*( AEy[i-1] + AMRy[i-1] + Aylist[i-1] + Thrust/(mlist[i-1] - M_dot*ts  ))*ts^2 + (Vylist[i-1] + Thrust*ts*sin(AV_theta))/(mlist[i-1]-M_dot*ts) 

R = sqrt(Xlist[i]^2 +Ylist[i]^2)

mlist[i] = mlist[i-1]-M_dot*ts

AEx[i] = -(mu/R^3)*Xlist[i]
AEy[i] = -(mu/R^3)*Ylist[i]

mmR[i] = sqrt(  (Ex1[i] - Xlist[i])^2 + (Ey1[i] - Ylist[i])^2  )
AMR[i] = (mu/6)/mmR[i]^3


AMRx[i] = AMR[i-1]*(Ex1[i]-Xlist[i])
AMRy[i] = AMR[i-1]*(Ey1[i]-Ylist[i])


Vxlist[i] = (Vxlist[i-1] + Thrust*ts*cos(AV_theta))/(mlist[i-1]-M_dot*ts) 
Vylist[i] = (Vylist[i-1] + Thrust*ts*sin(AV_theta))/(mlist[i-1]-M_dot*ts) 
 
Axlist[i] = AEx[i-1] + AMRx[i-1] + Axlist[i-1] + Thrust/(mlist[i-1] - M_dot*ts)  
Aylist[i] = AEy[i-1] + AMRy[i-1] + Aylist[i-1] + Thrust/(mlist[i-1] - M_dot*ts) 



# Vxlist[i] = Xlist[i]-Xlist[i-1]
# Vylist[i] = Ylist[i]-Ylist[i-1]
# 
# Axlist[i] = Vxlist[i]-Vxlist[i-1]
# Aylist[i] = Vxlist[i]-Vxlist[i-1]


}#ENDOF 

#Xlist
#Ylist

summary(lm(Ylist~Xlist))

plot(Xlist,Ylist, col = rgb(  (1:length(Ex1))/length(Ex1),0,0))

B_LIST = data.frame(Xlist,Ylist)


#library(reticulate)
#py_install("pandas")


data <- data.frame(data)


L = 1:length(data$R)

#linear
summary(lm(data$R~((L))))

# square
summary(  lm(data$R~cos(L)+sqrt(L) ))


plot(1:length(data$R),data$R, col = 'blue')
points(1:length(data$R), data$R -( 0 + 1.563e-02*(1:length(data$R))), col = 'red'  )
points(1:length(data$R), ( 6.692e+03 + 1.563e-02*(1:length(data$R))), col = 'black', lwd = 0.75  )
#points(1:length(data$R), ( 6.407e+03 + 4.611*sqrt(1:length(data$R))), col = 'purple', lwd = 0.5  )
#points(1:length(data$R),  6.407e+03 -6.407e-03*cos(w*L)+ 4.611*sqrt(L), col = 'violet', lwd = 0.5  )

bx1 = 0
bx2 = 5.58e4
bx3 = 1.04e5
bx4 = 1.5e5

by1 = 0
by2 = 1.08e3
by3 = 1.7e3
by4 = 2347

#Bezier Curve

Bx = 0
By = 0
for(i in 1:1.5e5){
I = 1-i/1.5e5
Bx[i] =  I*(I*(I*bx1+i*bx2/span) +i*(I*bx2+i*bx3/span) ) +i*(I*(I*bx2+i*bx3/span) + i*(I*bx3 + i*bx4/span) )
By[i] =  I*(I*(I*by1+i*by2/span) +i*(I*by2+i*by3/span) ) +i*(I*(I*by2+i*by3/span) + i*(I*by3 + i*by4/span) )

}

span = 1.5e5

R_list = data$


bx1 = 0
bx2 = span/2
bx3 = 3*span/4
bx4 = span

by1 = 0
by2 = (R_list[floor(span/2)])
by3 = (R_list[floor(3*span/4)])
by4 = (R_list[span-1])

#Bezier Curve
Bx = 0
By = 0

for (i in 1:span){    
  I = 1-i/span
Bx[i] = (I*(I*(I*bx1+i*bx2/span) +i*(I*bx2+i*bx3/span) ) +i*(I*(I*bx2+i*bx3/span) + i*(I*bx3 + i*bx4/span) ))
By[i] = (I*(I*(I*by1+i*by2/span) +i*(I*by2+i*by3/span) ) +i*(I*(I*by2+i*by3/span) + i*(I*by3 + i*by4/span) ))


}

plot(Bx, By)








