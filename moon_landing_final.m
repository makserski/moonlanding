clc
clear all                            
%%
%Lunar descent stage
%The lunar descent stage of the Apollo 11 mission involves taking the LM (Lunar
%Module) out of a circular orbit into an elliptical orbit, called a DOI (Descent
%Orbit Insertion) is performed. The perigee of the orbit is the closest
%approcah to the lunar surface, where the PDI (Powered Descent Initiation)
%is performed at roughly 15.2km. This is done by performing a retrograde
%maneuver. The velocity of the 
clc;
clear all;
Gconst = 6.67408e-11;               %Gravitational constant (m^3.kg^-1.s^-2)
M_moon = 7.34767309e22;             %Mass of moon (kg)
R_circular_orbit = 1848220;         %Radius of orbit around moon (m)
R_apogee = 1849754;                 %Apogee of the elliptical orbit (m)
R_perigee = 1752340;                %Perigee of the elliptical orbit (m)
R_moon = 1737100;                   %Radius of the moon (m)
Semi_major = 1801047;               %Semi-major axis distance (m)
m0 = 15200;                         %Initial mass of LM (kg)
mf = 7000;                          %Dry mass of LM (kg)
dmdt = 14.76;                       %Mass flow rate (k/s)
Ve = 3051;                          %Exhaust velocity (m/s)
Isp = 311;                          %Specific impulse (s)
g_moon = 1.62;                      %Gravitational field strength (m/s^2)
Thrust1 = 45000;                    %Maximum thrust output (N)
Angle = 90;                         %Angle in degrees
deltaV_DOI = (Gconst*M_moon/R_circular_orbit)^0.5*((2*R_perigee/(R_circular_orbit+R_perigee))^0.5 - 1);    %Hohman transfer equation
disp(abs(deltaV_DOI))               %DeltaV of DOI (m/s)
Burn_time = -(m0 - m0*exp(abs(deltaV_DOI)/Ve))/(dmdt*exp(abs(deltaV_DOI)/Ve));     %Burn time for PDI
disp(abs(Burn_time))                %Burn time for DOI (s)
e = (R_apogee-R_perigee)/(R_apogee+R_perigee);
disp(e)                             %Eccentricity of the orbit
V_apogee = sqrt(Gconst*M_moon*(1+e)/R_perigee);
disp(V_apogee)
m0_PDI = m0 - Burn_time*dmdt;       
disp(m0_PDI)            

h_lander = 0;                       %Initial height
M_lander = 15092;                   %Mass after DOI burn
F1=45040;                           %Lateral thrust
matrix_lander=zeros(1,7);
t0= 0; 
dt_lateral =0.1;                    %Timestep
massflow1=14.76*dt_lateral; 
time1 = 580;                        % Time after lateral burn
g = 1.62;
angle = 0;                      
v_lateral = 1609;                   %Initial lateral velocity
rho=0; 
s_lateral=0;                        %Iteration
hh_lateral=0;                       
totalT=0;                           %Time
for t = t0:dt_lateral:time1         
    a = getacceleration(F1,M_lander,g,angle,v_lateral,rho);     %Calculates acceleration
    h_lander = getheight(h_lander,dt_lateral,v_lateral,angle);  %Calculates the distance
    v_lateral=v_lateral-dt_lateral*a;                           %Calculates the velocity
    M_lander = M_lander - massflow1;                            %Calculates the mass
    matrix_lander = updatetable(s_lateral,a,v_lateral,h_lander,angle,hh_lateral,totalT,matrix_lander,M_lander); %Updates the matrix
    s_lateral=1+s_lateral;                                      %Iteration counter for mattrix 
    totalT=totalT+dt_lateral;       
end


figure (1)

plot(matrix_lander(:,6),matrix_lander(:,1)) %Plot acceleration against time

title('Lateral acceleration against time')

xlabel('Time (s)')

ylabel('Acceleration (m/s^2)')

figure (2)

plot(matrix_lander(:,6),matrix_lander(:,2)) %plot lateral velocity against time

title('Velocity against time')

xlabel('Time (s)')

ylabel('Velocity (m/s)')

figure (3)

plot(matrix_lander(:,6),matrix_lander(:,3)) %Plot displacement against time
title('Downrange distance against time')

xlabel('Time (s)')

ylabel('Distance (m) ')

figure(4)
plot(matrix_lander(:,6),matrix_lander(:,7)) %Plot mass against time
xlabel('Time (s)')
ylabel('Mass (kg)')


h_vertical = 15092;                         %Initial height
g_moon=1.62;                                %Gravitational strength of the moon
F1_vertical=5735;                           %Thrust 1
F2_vertical=14930;                          %Thrust 2
matrix_lander=zeros(1,6);
t0_vertical= 0; 
dt_vertical =0.1;                           %Timestep
massflow2=1.9*dt_vertical;                  
massflow3 = 4.9*dt_vertical
time1_vertical = 150;                       %Time after first stage
time2_vertical= 283;                        %landing time
angle_vertical = 0;                         %Needed to run acceleration func
v_vertical = 0;                             %Inital velocity
rho=0;                                      %Needed for acceleration
s_vertical=0;                               %Iteration
hh_vertical=0;                              %Needed for height plot

for t = t0_vertical:dt_vertical:time1_vertical
    a = getacceleration(F1_vertical,M_lander,g_moon,angle_vertical,v_vertical,rho); %Calculates the acceleration
    h_vertical = getheight(h_vertical,dt_vertical,v_vertical,angle_vertical);       %Calculate the altitude
    v_vertical=v_vertical+dt_vertical*a;                                            %Calculate the velocity
    M_lander = M_lander - massflow2;                                                %Calculate the mass
    matrix_lander = updatetable(s_vertical,a,v_vertical,h_vertical,angle_vertical,hh_vertical,totalT,matrix_lander,M_lander);   %Updates the matrix
    s_vertical=1+s_vertical;                                                        %Iteration counter for matrix 
    totalT=totalT+dt_vertical; 
end

for t = time1_vertical:dt_vertical:time2_vertical
    a = getacceleration(F2_vertical,M_lander,g_moon,angle_vertical,v_vertical,rho); %Calculates the acceleration
    h_vertical = getheight(h_vertical,dt_vertical,v_vertical,angle_vertical);       %Calculate the altitude
    v_vertical=v_vertical+dt_vertical*a;                                            %Calculate the velocity
    M_lander = M_lander - massflow3;                                                %Calculate the mass
    matrix_lander = updatetable(s_vertical,a,v_vertical,h_vertical,angle_vertical,hh_vertical,totalT,matrix_lander,M_lander);   %Updates the matrix
    s_vertical=1+s_vertical;                                                        %Iteration counter for the matrix
    totalT=totalT+dt_vertical;
end

figure (5)

plot(matrix_lander(:,6),matrix_lander(:,1)) %Plot acceleration against time 

title('Vertical acceleration against time')

xlabel('Time (s)')

ylabel('Acceleration (m/s^2)')

figure (6)

plot(matrix_lander(:,6),matrix_lander(:,2)) %Plot velocity against time

title('Vertical velocity against time')

xlabel('Time (s)')

ylabel('Velocity (m/s)')

figure (7)

plot(matrix_lander(:,6),matrix_lander(:,3)) %Plot height against time 
ylim([0 15500])
title('Altitude against time')

xlabel('Time (s)')

ylabel('Altitude (m) ')

figure(8)
plot(matrix_lander(:,6),matrix_lander(:,7)) %Plot mass against time
xlabel('Time (s)')
ylabel('Mass (kg)')