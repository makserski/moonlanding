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
%maneuver. The velocity of the 
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
Thrust2 = 27000;                    %60% of max thrust
Thrust3 = 4500;                     %10% of max thrust
Angle = 90;                         %Angle in degrees
deltaV_DOI = (Gconst*M_moon/R_circular_orbit)^0.5*((2*R_perigee/(R_circular_orbit+R_perigee))^0.5 - 1);    %Hohman transfer equation
%disp(abs(deltaV_DOI))               %DeltaV of DOI (m/s)
Burn_time = -(m0 - m0*exp(abs(deltaV_DOI)/Ve))/(dmdt*exp(abs(deltaV_DOI)/Ve));     %Burn time for PDI
%disp(abs(Burn_time))                %Burn time for DOI (s)
e = (R_apogee-R_perigee)/(R_apogee+R_perigee);
%disp(e)                             %Eccentricity of the orbit
V_apogee = sqrt(Gconst*M_moon*(1+e)/R_perigee);
%disp(V_apogee)
m0_PDI = m0 - Burn_time*dmdt;       
%disp(m0_PDI)            

h_l = 0; %initial height
M_lander = 15092; %inital mass
F1=45040; %thrust 1
matrix_lander=zeros(1,7);
t0= 0; 
dt_lateral = 0.1; %timestep
massflow1=14.76*dt_lateral; 
time1 = 580; % time after first stage
g = 1.62;
angle = 0; %needed to run acceleration func
v_lateral = 1609; %inital velocity
rho=0; %needed for acceleratiom
s_lateral=0; %iteration
hh_lateral=0; %needed
totalT=0; %time (neede)
for t = t0:dt_lateral:time1
    a = getacceleration(F1,M_lander,g,angle,v_lateral,rho);
    h_l = getheight(h_l,dt_lateral,v_lateral,angle);
    v_lateral=v_lateral-dt_lateral*a;
    M_lander = M_lander - massflow1;
    matrix_lander = updatetable(s_lateral,a,v_lateral,h_l,angle,hh_lateral,totalT,matrix_lander,M_lander);
    s_lateral=1+s_lateral; %iteration counter for mattrix 
    totalT=totalT+dt_lateral; %x axis (stops it plotting in ms) (look at the plot functions)
end

figure (1)
plot(matrix_lander(:,6),matrix_lander(:,1)) %plot time acceleration
title('Lateral acceleration against time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')

figure (2)
plot(matrix_lander(:,6),matrix_lander(:,2)) %plot time velocity
title('Velocity against time')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

figure (3)
plot(matrix_lander(:,6),matrix_lander(:,3)) %plot height time 
title('Downrange distance against time')
xlabel('Time (s)')
ylabel('Distance (m) ')

figure(4)
plot(matrix_lander(:,6),matrix_lander(:,7))
title('Mass against time')
xlabel('Time (s)')
ylabel('Mass (kg)')

% works at F1 =15200 
                                        % F2=26500
%                                         M = 15200 %inital mass
%                                         F1=14100 %thrust 1
%                                         F2=26500 %thrust 2

h_vertical = 15092; %initial height
g_moon=1.62; %gravity
M_lander = 8000; %inital mass after lateral burn
F1_vertical=7210; %thrust 1
F2_vertical=18000; %thrust 2
matrix_lander=zeros(1,6);
t0_vertical= 0; 
dt_vertical =0.1; %timestep
massflow2=1.5*dt_vertical; 
massflow3 = 5*dt_vertical;
time1_vertical = 150; % time after first stage
time2_vertical= 283;%landing time
angle_vertical = 0; %needed to run acceleration func
v_vertical = 0; %inital velocity
rho=0; %needed for acceleratiom
s_vertical=0; %iteration
hh_vertical=0; %needed

for t = t0_vertical:dt_vertical:time1_vertical
    a = getacceleration(F1_vertical,M_lander,g_moon,angle_vertical,v_vertical,rho);
    h_vertical = getheight(h_vertical,dt_vertical,v_vertical,angle_vertical);
    v_vertical=v_vertical+dt_vertical*a;
    M_lander = M_lander - massflow2;
    matrix_lander = updatetable(s_vertical,a,v_vertical,h_vertical,angle_vertical,hh_vertical,totalT,matrix_lander,M_lander);
    s_vertical=1+s_vertical; %iteration counter for mattrix 
    totalT=totalT+dt_vertical; %x axis (stops it plotting in ms) (look at the plot functions)
end

for t = time1_vertical:dt_vertical:time2_vertical
    a = getacceleration(F2_vertical,M_lander,g_moon,angle_vertical,v_vertical,rho);
    h_vertical = getheight(h_vertical,dt_vertical,v_vertical,angle_vertical);
    v_vertical=v_vertical+dt_vertical*a;
    M_lander = M_lander - massflow3;
    matrix_lander = updatetable(s_vertical,a,v_vertical,h_vertical,angle_vertical,hh_vertical,totalT,matrix_lander,M_lander);
    s_vertical=1+s_vertical;
    totalT=totalT+dt_vertical;
end

figure (5)
plot(matrix_lander(:,6),matrix_lander(:,1)) %plot time acceleration
title('Vertical acceleration against time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')

figure (6)
plot(matrix_lander(:,6),matrix_lander(:,2)) %plot time velocity
title('Vertical velocity against time')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

figure (7)
plot(matrix_lander(:,6),matrix_lander(:,3)) %plot height time 
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