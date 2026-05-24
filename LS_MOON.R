




set.seed(3)
# standard Value
R0 = rnorm(n=length(MOON_TRAK_1$a),mean=60.33*6371.01,sd=10)
R= mean(R0)
R
plot(1:length(MOON_TRAK_1$a),R0)
# extract LAT/LON
LON = MOON_TRAK_1$LON
LAT = MOON_TRAK_1$LAT

# convert to polar coordinates
X = R*cos(LON)
Y = R*sin(LON)*cos(LAT)
Z = R*sin(LON)*sin(LAT)

# extract times(ZULU)
Tzu = MOON_TRAK_1$time_zulu

# Run least squares over coordinates
LX = lm(X~Tzu)
LY = lm(Y~Tzu)
LZ = lm(Z~Tzu)

LX$coefficients[2]
LY$coefficients[2]
LZ$coefficients[2]

summary(lm(X~Tzu))
summary(lm(Y~Tzu))
summary(lm(Z~Tzu))


#magnitude of Rdot from LS of R using 2nd coeff.
Rdot = sqrt(LX$coefficients[2]^2 +LY$coefficients[2]^2+LZ$coefficients[2]^2)
# convert to seconds
Rdot/60




#initialize velocities
Xdot = 0
Ydot = 0
Zdot = 0
# compute velocities
for( i in 1:length(MOON_TRAK_1$time_zulu)){

Xdot[i] = (X[i+1]-X[i])/(Tzu[i+1]-Tzu[i])
Ydot[i] = (Y[i+1]-Y[i])/(Tzu[i+1]-Tzu[i])
Zdot[i] = (Z[i+1]-Z[i])/(Tzu[i+1]-Tzu[i])
}
LVX = lm(Xdot~Tzu)
LVY = lm(Ydot~Tzu)
LVZ = lm(Zdot~Tzu)
summary(lm(Xdot~Tzu))
summary(lm(Ydot~Tzu))
summary(lm(Zdot~Tzu))
LRV = sqrt(LVX$coefficients[1]^2+LVY$coefficients[1]^2+LVZ$coefficients[1]^2)/3600#converted to seconds
LRV


# compare R^ to standard value ~60*Re
R_hat = sqrt(LX$coefficients[1]^2 +LY$coefficients[1]^2+LZ$coefficients[1]^2)/R
R_hat

## prepare Wr = V to compute tangential velocity
# compute omega using least squares results
W = 0
for( i in 1:length(MOON_TRAK_1$time_zulu)){

  # arccos(AB / |A||B|)
  
  Q = (X[i]*X[i+1]+Y[i]*Y[i+1]+Z[i]*Z[i+1])# dot product
  Rq = sqrt(X[i]^2+Y[i]^2+Z[i]^2)*sqrt(X[i+1]^2+Y[i+1]^2+Z[i+1]^2)
  t = Tzu[i+1]-Tzu[i]
  theta = acos(Q/Rq)
  W[i] = theta/t
}

lW = lm(W~Tzu)
summary(lW)
omega = lW$coefficients[1]/3600# converted to seconds
omega

TV = omega*sqrt(LX$coefficients[1]^2 +LY$coefficients[1]^2+LZ$coefficients[1]^2)
TV

TV1 = omega*R## will this always be most accurate?
TV1
plot(Tzu[1:length(Tzu)],W)

SV0 = c(R_hat,TV1)# R in earth radii, Tangential velocity
SV0

MOON_TRAK_1$W = W



# SV = c(LX$coefficients[1],LY$coefficients[1],LZ$coefficients[1],
#        LX$coefficients[2],LY$coefficients[2],LZ$coefficients[2])
# SV


hist(X)
hist(Y)
hist(Z)

hist(Xdot)
hist(Ydot)
hist(Zdot)


plot(Tzu,X)
plot(Tzu,Y)
plot(Tzu,Z)


cov(X,Y)
cov(X,Z)
cov(Y,Z)


lm(R0^2~X^2+Y^2+Z^2)


