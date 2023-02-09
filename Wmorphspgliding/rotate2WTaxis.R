rotate2WTaxis <- function(varying.data,len_data){
   ### ------------ Rotate X,Y,Z Forces and Moments to be Drag, Lift, Sideforce, Yaw, Roll and Pitch
  # -----------Date: 08-Sep-17 ---------------------
  
  ### ---------- 1: CONSTANTS -----------------------
  o1           = (0.550-0.250)*25.4 #mm - measured from solidworks model offset of circular region from back plate
  r1           = 0.95*25.4          #mm - measured from solidworks model radius of circular region
  x0           = 0.0587             #mm - offset of LC electrical center
  y0           = 0.1512             #mm - offset of LC electrical center
  len_0        = sqrt(x0^2+y0^2)    #mm - length of LC offset
  gam          = (180/pi)*atan(y0/x0) + 42; #angle of that offset
  z0           = 36.3232            #mm - offset of LC electrical center
  h1           = 22.2               #mm - length of arm that protrudes from support
  h1_sd        = 0.001;             #from caliper
  h35          = 48.3-(32.4*cos((pi/180)*55)) #horizontal distance from attachment of support to wing rod along 35deg piece
  h35_sd       = sqrt(0.001^2+(0.001*cos((pi/180)*55))^2) #absolute uncertainty of this measurement
  v35          = 29.9+31.1*sin((pi/180)*35)-5.8*sin((pi/180)*55) #vertical distance from attachment of support to wing rod along 35deg piece
  v35_sd       = sqrt(0.001^2+(0.001*sin((pi/180)*55))^2+(0.001*sin((pi/180)*35))^2) #absolute uncertainty of this measurement
  LC_wallinset = 0.5*25.4 - 9.1; #how far the load cell is installed in the wall 0.5"(from drawing)-9.1mm (how far the spacer sticks out) 
  airfoil_t    = (3+(1/16))*25.4; #the thickness of the airfoil glued on wall
  varying.data$theta  = asin((varying.data$I-o1)/r1) 
  
  # ----------------------------------------------------- 2: OFFSET CALCULATIONS ------------------------------------------------

  # -------- HORIZONTAL AND VERTICAL OFFSET CALCULATIONS ----------------
  varying.data$horiz_offset0 = (len_0*cos((pi/180)*(180-gam))  + ((r1+h35+h1+varying.data$A*cos((pi/180)*35))*cos(varying.data$theta)))/1000
  varying.data$vert_offset0  = (len_0*sin((pi/180)*(180-gam)) + (varying.data$A*sin((pi/180)*35)+v35))/1000
  varying.data$arm           = sqrt(varying.data$horiz_offset0^2+varying.data$vert_offset0^2)
  varying.data$armangle      = abs(atan(varying.data$vert_offset0/varying.data$horiz_offset0)) #need absolute due to working in the 4th quadrant
  varying.data$horiz_offset  = varying.data$arm*cos(varying.data$armangle-((pi/180)*varying.data$curr.angle))
  varying.data$vert_offset   = varying.data$arm*sin(varying.data$armangle-((pi/180)*varying.data$curr.angle))
  
  varying.data$armCOG          = sqrt((varying.data$horiz_offset0+varying.data$COG.chord)^2+(varying.data$vert_offset0+varying.data$COG.z)^2)
  varying.data$armangleCOG     = abs(atan((varying.data$vert_offset0+varying.data$COG.z)/(varying.data$horiz_offset0+varying.data$COG.chord))) #need absolute due to working in the 4th quadrant
  varying.data$horiz_offsetCOG = varying.data$armCOG*cos(varying.data$armangleCOG - ((pi/180)*varying.data$curr.angle))
  
  # ------ ABSOLUTE UNCERTAINTY ---------------
  horiz_offset0_sd                = (sqrt(h35_sd^2+h1_sd^2+(0.001*cos((pi/180)*35))^2)*cos(varying.data$theta))/1000
  vert_offset0_sd                 = (sqrt(((0.001*sin((pi/180)*35))^2+v35_sd^2)))/1000
  arm_sd                          = sqrt((varying.data$horiz_offset*horiz_offset0_sd)^2+(varying.data$vert_offset*vert_offset0_sd)^2)/(varying.data$arm) 
  armangle_sd                     = abs(1/((varying.data$vert_offset0/varying.data$horiz_offset0)^2+1))*(varying.data$vert_offset0/varying.data$horiz_offset0)*
                                    sqrt((vert_offset0_sd/varying.data$vert_offset0)^2+(horiz_offset0_sd/varying.data$horiz_offset0)^2)
  varying.data$horiz_offset_sd    = abs(varying.data$horiz_offset)*sqrt((arm_sd/varying.data$arm)^2+(armangle_sd*sin(varying.data$armangle)/cos(varying.data$armangle - ((pi/180)*varying.data$curr.angle)))^2)
  varying.data$vert_offset_sd     = abs(varying.data$vert_offset)*sqrt((arm_sd/varying.data$arm)^2+(armangle_sd*cos(varying.data$armangle)/sin(varying.data$armangle - ((pi/180)*varying.data$curr.angle)))^2)
  
  armCOG_sd                       = sqrt((((varying.data$horiz_offset0+varying.data$COG.chord)*(horiz_offset0_sd+varying.data$COG.chord.std))^2)
                                         +((varying.data$vert_offset0+varying.data$COG.z)*(vert_offset0_sd+varying.data$COG.z.std))^2)/(varying.data$armCOG) 
  armangleCOG_sd                  = abs(1/(((varying.data$vert_offset0+varying.data$COG.z)/(varying.data$horiz_offset0+varying.data$COG.chord))^2+1))*((varying.data$vert_offset0+varying.data$COG.z)/(varying.data$horiz_offset0+varying.data$COG.chord))*
                                    sqrt(((vert_offset0_sd+varying.data$COG.z.std)/(varying.data$vert_offset0+varying.data$COG.z))^2+((horiz_offset0_sd+varying.data$COG.chord.std)/(varying.data$horiz_offset0+varying.data$COG.chord))^2)
  varying.data$horiz_offsetCOG_sd = abs(varying.data$horiz_offsetCOG)*sqrt((armCOG_sd/varying.data$armCOG)^2+(armangleCOG_sd*sin(varying.data$armangleCOG)/cos(varying.data$armangleCOG - ((pi/180)*varying.data$curr.angle)))^2)
  
  # ----------- SIDE OFFSET ------------------
  for (i in 1:len_data){
    if (is.na(varying.data$F[i])){
      varying.data$side_offset[i]    = (z0+varying.data$H[i] + LC_wallinset)/1000
      varying.data$side_offset_sd[i] = (sqrt(2)*0.001)/1000; 
    } else { #both H and wallinset are measuered by caliper
      varying.data$side_offset[i]    = (z0+LC_wallinset+ airfoil_t + varying.data$F[i] + (varying.data$A[i]*cos((pi/180)*35)*sin(varying.data$theta[i])))/1000
      varying.data$side_offset_sd[i] = sqrt(0.001^2+0.001^2+0.001^2+(0.001*cos((pi/180)*35)*sin(varying.data$theta[i]))^2)/1000}}

  # ------ ROTATE TO NEW COORDINATE SYSTEM - AERODYNAMIC LIFT, DRAG -------------
  varying.data$Drag      = - varying.data$Fx*cos((42-varying.data$curr.angle)*pi/180) + varying.data$Fy*sin((42-varying.data$curr.angle)*pi/180) - varying.data$Wing.Weight*9.81*(sin((varying.data$curr.angle)*pi/180));
  varying.data$Lift      = - varying.data$Fx*sin((42-varying.data$curr.angle)*pi/180) - varying.data$Fy*cos((42-varying.data$curr.angle)*pi/180) + varying.data$Wing.Weight*9.81*(1-cos((varying.data$curr.angle)*pi/180));
  #no side force due to the drift error
  
  #       Following Anderson's convention in Intro to Flight these are the opposite direction. Pitch is negative 
  #       as the load cell z axis positive was initially into the wall however Anderson convention has y along the length of the right wing.
  varying.data$Roll0     = -(- varying.data$Mx*cos((42-varying.data$curr.angle)*pi/180) + varying.data$My*sin((42-varying.data$curr.angle)*pi/180) +
                               varying.data$Wing.Weight*9.81*(varying.data$side_offset+varying.data$COG.span)*(1-cos((varying.data$curr.angle)*pi/180)));
  varying.data$Yaw0      = -(- varying.data$Mx*sin((42-varying.data$curr.angle)*pi/180) - varying.data$My*cos((42-varying.data$curr.angle)*pi/180) +
                               varying.data$Wing.Weight*9.81*(varying.data$side_offset+varying.data$COG.span)*(sin((varying.data$curr.angle)*pi/180)));
  varying.data$Pitch0    = -(  varying.data$Mz + varying.data$Wing.Weight*9.81*(varying.data$horiz_offsetCOG-(varying.data$horiz_offset0 + varying.data$COG.chord))); #negative due to the orientation of the load cell 
  
  ##-- RELOCATE ALL MOMENTS TO THE SHOULDER ATTACHMENT LOCATION AND REMOVE EFFECTS OF WEIGHT WITH CENTER OF MASS
 
  varying.data$Pitch     = varying.data$Pitch0 + varying.data$Lift*varying.data$horiz_offset - varying.data$Drag*varying.data$vert_offset 
  
  # ------------- Uncertainty -----------
  #this should be absolute uncertainty for future calculations --- 
  varying.data$Drag_sd        = sqrt(((varying.data$Fx_std*abs(cos((42-varying.data$curr.angle)*pi/180)))^2)+((varying.data$Fy_std*abs(sin((42-varying.data$curr.angle)*pi/180)))^2) + (0.0001*abs(9.81*(sin((varying.data$curr.angle)*pi/180))))^2);
  varying.data$Lift_sd        = sqrt(((varying.data$Fx_std*abs(sin((42-varying.data$curr.angle)*pi/180)))^2)+((varying.data$Fy_std*abs(cos((42-varying.data$curr.angle)*pi/180)))^2) + (0.0001*abs(9.81*(1-cos((varying.data$curr.angle)*pi/180))))^2);
  
  varying.data$Roll0_sd       = sqrt(((varying.data$Mx_std*abs(cos((42-varying.data$curr.angle)*pi/180)))^2)+((varying.data$My_std*abs(sin((42-varying.data$curr.angle)*pi/180)))^2) +
                                      (abs(varying.data$Wing.Weight*9.81*(varying.data$side_offset+varying.data$COG.span)*(1-cos((varying.data$curr.angle)*pi/180)))*
                                      sqrt(0.0001/(varying.data$Wing.Weight*9.81)+(sqrt(varying.data$side_offset_sd^2+varying.data$COG.span.std^2)/(varying.data$side_offset+varying.data$COG.span))))^2);
  
  varying.data$Yaw0_sd        = sqrt(((varying.data$Mx_std*abs(sin((42-varying.data$curr.angle)*pi/180)))^2)+((varying.data$My_std*abs(cos((42-varying.data$curr.angle)*pi/180)))^2)+
                                      (abs(varying.data$Wing.Weight*9.81*(varying.data$side_offset+varying.data$COG.span)*(sin((varying.data$curr.angle)*pi/180)))*
                                      sqrt(0.0001/(varying.data$Wing.Weight*9.81)+(sqrt(varying.data$side_offset_sd^2+varying.data$COG.span.std^2)/(varying.data$side_offset+varying.data$COG.span))))^2);

  varying.data$Pitch0_sd      = varying.data$Mz_std + abs(varying.data$Wing.Weight*9.81*(varying.data$horiz_offsetCOG-varying.data$horiz_offset0 + varying.data$COG.chord))*
                                sqrt((0.0001/(varying.data$Wing.Weight*9.81))^2+
                                       (sqrt(varying.data$horiz_offsetCOG_sd^2+varying.data$horiz_offset_sd^2+varying.data$COG.chord.std^2)/(varying.data$horiz_offsetCOG-varying.data$horiz_offset0 + varying.data$COG.chord))^2)
  
  varying.data$Pitch_sd  = sqrt(varying.data$Pitch0_sd^2 + 
                                  (abs(varying.data$Lift*varying.data$horiz_offset)*sqrt((varying.data$Lift_sd/varying.data$Lift)^2
                                                                                         +(varying.data$horiz_offset_sd/varying.data$horiz_offset)^2))^2 +
                                  (abs(varying.data$Drag*varying.data$vert_offset)*sqrt((varying.data$Drag_sd/varying.data$Drag)^2
                                                                                        +(varying.data$vert_offset_sd/varying.data$vert_offset)^2))^2);
  return(varying.data)
}
