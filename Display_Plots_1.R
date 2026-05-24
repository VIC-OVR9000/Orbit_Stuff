


library(tidyverse)
library(ggplot2)
library(cowplot)
library(gganimate)



##Model of Actual Orbit
PLOT_MODEL_ORBIT = function(){
paramLabs = paste0("e0 =", e0, "\n", "P_o =", P_o, "\n"  , "A_o =", A_o,"\n" )
MODEL_0_plot = ggplot(MODEL_0, aes(Ex1,Ey1) ) +
               geom_point()+
              annotate("text", x = 2*Re , y = 0, label = paramLabs )+
              geom_point(x = Re*cos(t),y=Re*sin(t), size = 1)+
              labs(title = "Orbit Model" )
return(MODEL_0_plot)
}#ENDOF function PLOT MODEL ORBIT



PLOT_RDOT_0 = function(){
  
  f_dot = (P_o+Re)*e0*w*sin(w*t)/(1+e0*cos(w*t))^2
  #plot(t,f)
  plot(t,r_dot, lwd = 5)
  points(t,f_dot, col='red', lwd = 1)
  
 
}#ENDOF PLOT_RDOT_0




PLOT_FINAL_ORBIT_MODEL = function(){
  
  paramLabs = paste0("ef =", ef, "\n", "P_f =", P_f, "\n"  , "A_f =", A_f,"\n" )
  MODEL_f_plot = ggplot(MODEL_f, aes(Ex2,Ey2) ) +
    geom_point()+
    annotate("text", x = 2*Re , y = 0, label = paramLabs )+
    geom_point(x = Re*cos(tf),y=Re*sin(tf), size = 1)+
    labs(title = "Orbit Model" )
  return(MODEL_f_plot)
  
}#ENDOF plot final model



PLOT_RDOT_F = function(){
  
  f_dot2 = (P_f+Re)*ef*w1*sin(w1*tf)/(1+ef*cos(w1*tf))^2
  #plot(tf,f2)
  plot(tf,r_dot2, lwd = 5)
  points(tf,f_dot2, col='red', lwd = 1)
  
}#ENDOF PLOT_RDOT_F


PLOT_ELLIPSE_DVslider = function(input){
  
  V_select = input# input
  
  e_n = (Re+P_o)*V_select^2/mu -1
  
  A_n = ( e_n*P_o + 2*e_n*Re+P_o )/( 1 - e_n )
  
  #build new ellipse for display
  a_n = (P_o + A_n +2*Re)/2
  
  c_n = a_n - P_o - Re
  
  b_n = sqrt(a_n^2 - c_n^2  )
  
  T_n = 2*pi*sqrt(a_n^3/mu)
  
  w_n = 2*pi/T_n
  
  t_n = seq(0 , T_n, ts)
  
  plot( a_n*cos(w_n*t_n)+c_n, b_n*sin(w_n*t_n)       )
  points(Ex1,Ey1, col = 'blue')
  
  
  
  
  
  
  plot( a_n*cos(w_n*t_n)+c_n, b_n*sin(w_n*t_n)       )
  points(Ex1,Ey1, col = 'blue')
  
  
  
}#ENDOF f-input DVslider













