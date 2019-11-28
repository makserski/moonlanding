function [MASS, TIME, VELOCITY] = moon_landing_final()
% clc
% clear all                            


%%
%Lunar descent stage
%The lunar descent stage of the Apollo 11 mission involves taking the LM (Lunar
%Module) out of a circular orbit into an elliptical orbit, called a DOI (Descent
%Orbit Insertion) is performed. The perigee of the orbit is the closest
%approcah to the lunar surface, where the PDI (Powered Descent Initiation)
%is performed at roughly 15.2km. This is done by performing a retrograde
%maneuver. The velocity of the LM at the perigee is calculated in this
%script. The script plots the acceleration, mass, velocity and
%altitude/displacement of the LM at any time.
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
%disp(abs(deltaV_DOI))              %DeltaV of DOI (m/s)
Burn_time = -(m0 - m0*exp(abs(deltaV_DOI)/Ve))/(dmdt*exp(abs(deltaV_DOI)/Ve));     %Burn time for PDI
%disp(abs(Burn_time))               %Burn time for DOI (s)
e = (R_apogee-R_perigee)/(R_apogee+R_perigee);
%disp(e)                            %Eccentricity of the orbit
V_apogee = sqrt(Gconst*M_moon*(1+e)/R_perigee);
%disp(V_apogee)
m0_PDI = m0 - Burn_time*dmdt;       
%disp(m0_PDI)            

h_l = 0;                            %Initial displacement
M_lander = 15092;                   %Inital mass
F1=45040;                           %Thrust1 (N)
matrix_lander=zeros(1,7);
t0= 0; 
dt_lateral = 0.1;                   %Timestep
massflow1=14.76*dt_lateral; 
time1 = 580;                        %Time after lateral burn
g = 1.62;                           %Gravitational strength of the moon
angle = 0;                          %Needed to run acceleration func
v_lateral = 1609;                   %Inital velocity
rho=0;                              %Needed for acceleration
s_lateral=0;                        %Iteration
hh_lateral=0;                       %Needed for the displacement
totalT=0;                           %Time
for t = t0:dt_lateral:time1
    a = getacceleration(F1,M_lander,g,angle,v_lateral,rho); %Calculates the acceleration
    h_l = getheight(h_l,dt_lateral,v_lateral,angle);        %Calculates the altitude
    v_lateral=v_lateral-dt_lateral*a;                       %Calculates the velocity
    M_lander = M_lander - massflow1;                        %Calculates the mass
    matrix_lander = updatetable(s_lateral,a,v_lateral,h_l,angle,hh_lateral,totalT,matrix_lander,M_lander);  %Updates the matrix
    s_lateral=1+s_lateral;                                  %Iteration counter for mattrix 
    totalT=totalT+dt_lateral;                               %Time counter
end

%tiledlayout(3,3);
[X,Y] = peaks(30);

nexttile
%figure (1)
plot(matrix_lander(:,6),matrix_lander(:,1))     %Plot of acceleration against time
title('Lateral acceleration against time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')

nexttile
%figure (2)
plot(matrix_lander(:,6),matrix_lander(:,2))     %Plot of velocity against time
title('Velocity against time')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

nexttile
%figure (3)
plot(matrix_lander(:,6),matrix_lander(:,3))     %Plot of height against time 
title('Downrange distance against time')
xlabel('Time (s)')
ylabel('Distance (m) ')

nexttile
%figure(4)
plot(matrix_lander(:,6),matrix_lander(:,7))
title('Mass against time')
xlabel('Time (s)')
ylabel('Mass (kg)')



h_vertical = 15240;                         %Initial height(m)
g_moon=1.62;                                %Gravity (m/s^2)
F1_vertical=5700;                           %Thrust 2 (N)
F2_vertical=14930;                          %Thrust 3 (N)
matrix_lander=zeros(1,6);
t0_vertical= 0; 
dt_vertical =0.1;                           %Timestep
massflow2=1.9*dt_vertical;                  %Mass flow rate (kg/s)
massflow3 = 4.9*dt_vertical;
time1_vertical = 150; 
time2_vertical= 283;
angle_vertical = 0;                         %Needed to run acceleration function
v_vertical = 0;                             %Inital velocity
rho=0;                                      %Needed for acceleration
s_vertical=0;                               %Iteration
hh_vertical=0;                              %Needed for height plot

for t = t0_vertical:dt_vertical:time1_vertical
    a = getacceleration(F1_vertical,M_lander,g_moon,angle_vertical,v_vertical,rho); %Calculates the acceleration
    h_vertical = getheight(h_vertical,dt_vertical,v_vertical,angle_vertical);       %Calculates the altitude
    v_vertical=v_vertical+dt_vertical*a;                                            %Calculates the velocity
    M_lander = M_lander - massflow2;                                                %Calculates the mass
    matrix_lander = updatetable(s_vertical,a,v_vertical,h_vertical,angle_vertical,hh_vertical,totalT,matrix_lander,M_lander);   %Updates the matrix
    s_vertical=1+s_vertical;                                                        %Iteration counter for matrix 
    totalT=totalT+dt_vertical;                                                      %Time counter
end

for t = time1_vertical:dt_vertical:time2_vertical
    a = getacceleration(F2_vertical,M_lander,g_moon,angle_vertical,v_vertical,rho); %Calculates the acceleration
    h_vertical = getheight(h_vertical,dt_vertical,v_vertical,angle_vertical);       %Calculates the altitude
    v_vertical=v_vertical+dt_vertical*a;                                            %Calculates the velocity
    M_lander = M_lander - massflow3;                                                %Calculatesthe mass
    matrix_lander = updatetable(s_vertical,a,v_vertical,h_vertical,angle_vertical,hh_vertical,totalT,matrix_lander,M_lander);   %Updates the matrix
    s_vertical=1+s_vertical;                                                        %Iteration counter for matrix
    totalT=totalT+dt_vertical;                                                      %Time counter
end

nexttile
%figure (5)
plot(matrix_lander(:,6),matrix_lander(:,1)) %Plot of acceleration against time
title('Vertical acceleration against time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')

nexttile
%figure (6)
plot(matrix_lander(:,6),matrix_lander(:,2)) %Plot of velocity against time
title('Vertical velocity against time')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

nexttile
%figure (7)
plot(matrix_lander(:,6),matrix_lander(:,3)) %Plot of height against time 
ylim([0 15500])
title('Altitude against time')
xlabel('Time (s)')
ylabel('Altitude (m) ')

MASS = matrix_lander(:,7);
TIME = matrix_lander(:,6);
VELOCITY = matrix_lander(:,2);

MASS = transpose(MASS);
TIME = transpose(TIME);
VELOCITY = transpose(VELOCITY);
end
