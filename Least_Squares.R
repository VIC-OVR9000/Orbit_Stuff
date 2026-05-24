


INT = 1:3000
######## LS using normal interval
tSIM = t[INT] 
set.seed(3)
xSIM = n1*r0[INT]*cos(w*tSIM) +d + rnorm(length(tSIM),0,2.5)
ySIM = n2*r0[INT]*sin(w*tSIM) + rnorm(length(tSIM),0,2.5)
plot(xSIM,ySIM,col='red', main = "SIM(red), analytical(black)")
points(Ex1[INT],Ey1[INT])
length(xSIM)
length(ySIM)

LIN_SIM = data.frame(xSIM,xSIM^2,d,ySIM,ySIM^2,1)

xREG = glm( (ySIM)^2 ~ xSIM^2+xSIM)

xREG

par(mfrow = c(2,1))
plot(xSIM,ySIM,col = 'red', main = "SIM(red) vs fitted(blue)")
points(xSIM,xREG$fitted.values^0.5,col = 'blue')
#points(Ex1,Ey1)



plot2_res <- plot(INT,xREG$residuals^0.5)

par(mfrow = c(1,1))


plot(Ex1[INT],Ey1[INT], main = "Analytical(black), vs fitted(blue)")
points(xSIM,xREG$fitted.values^0.5,col = 'blue')

B = xREG$coefficients
B

Y_hat_sq = B[1]+B[2]*(n1*r0*cos(w*t) +d + rnorm(length(t),0,2.5))

plot((n1*r0*cos(w*t) +d + rnorm(length(t),0,2.5)),Y_hat_sq^0.5, col = 'blue')
points(Ex1,Ey1, lwd = 0.5)
points(xSIM,xREG$fitted.values^0.5,col = 'red')



### sim data adjusted for full time-span
plot((n1*r0*cos(w*t) +d + rnorm(length(t),0,2.5)),Y_hat_sq^0.5, col = 'blue', xlim = c(200,1150),ylim = c(0,700))
points(Ex1,Ey1, lwd = 0.5)
#points(xSIM,xREG$fitted.values^0.5,col = 'yellow')
points(xSIM,ySIM, col = 'purple')

RR = sqrt(Ex1^2 + Ey1^2) - sqrt( ((n1*r0*cos(w*t) +d + rnorm(length(t),0,2.5)))^2   + Y_hat_sq   )

plot(t, RR )


xSIM_dot = 0
for(i in 1:length(xSIM)){
xSIM_dot[i] = xSIM[i+1]-xSIM[i]
}#ENDOF xSIM dot

xSIM_ddot = 0
for(i in 1:length(xSIM)){
  xSIM_ddot[i] = xSIM_dot[i+1]-xSIM_dot[i]
}#ENDOF xSIM ddot


ySIM_dot = 0
for(i in 1:length(ySIM)){
  ySIM_dot[i] = ySIM[i+1]-xSIM[i]
}#ENDOF xSIM dot

ySIM_ddot = 0
for(i in 1:length(xSIM)){
  ySIM_ddot[i] = ySIM_dot[i+1]-ySIM_dot[i]
}#ENDOF xSIM ddot


rSIM = sqrt(xSIM^2 + ySIM^2)
one = rep(1,length(xSIM))
SIM_diff = data.frame(one, rSIM,xSIM, xSIM_dot,xSIM_ddot,ySIM,ySIM_dot,ySIM_ddot)

xSIMsq = xSIM^2
ySIMsq = ySIM^2
xySIM = xSIM*ySIM

for(i in 1:2000){

rI = sample(seq(-1,1,0.001),2000,replace = TRUE)
seedLIST = 0

seedSETlist = 0

seedA = set.seed(sample(0:1000,1,replace = FALSE))
seedLIST[0] = seedA
A = rI
seedB = set.seed(sample(0:1000,1,replace = FALSE))
seedLIST[1] =seedB
B = rI
seedC = set.seed(sample(0:1000,1,replace = FALSE))
seedLIST[2]= seedC
C = rI
seedD = set.seed(sample(0:1000,1,replace = FALSE))
seedLIST[3] = seedD
D = rI
seedE = set.seed(sample(0:1000,1,replace = FALSE))
seedLIST[4] = seedE
E = rI
seedF_xy = set.seed(sample(0:1000,1,replace = FALSE))
seedLIST[5] = seedF_xy
F_xy = rI

seedSETlist[i] = seedLIST

MSElist = 0

F_XY = A*xSIMsq + B*ySIMsq +C*xySIM +D*xSIM + E*ySIM+F_xy



}#ENDOF allipse for loop




write.csv(SIM_diff, "model_SIM.csv")








