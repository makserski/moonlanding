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
R_perigee = 1737200;                %Perigee of the elliptical orbit (m)
R_moon = 1752340;                   %Radius of the moon (m)
Semi_major = 1801047;               %Semi-major axis distance (m)
m0 = 15200;                         %Initial mass of LM (kg)
mf = 7000;                          %Dry mass of LM (kg)
dmdt = 14.76;                       %Mass flow rate (k/s)
Ve = 3051;                          %Exhaust velocity (m/s)
Isp = 311;                          %Specific impulse (s)
g_moon = 1.62;                      %Gravitational field strength (m/s^2)
Thrust1 = 45000;                    %Maximum thrust output (N)
Thrust2 = 0;                    %60% of max thrust
Thrust3 = 4500;                     %10% of max thrust
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
disp(m0_PDI)                        %Mass of LM after DOI burn
%%
%Initial conditions
h = 15240;                          %Starting height (m)
v = abs(V_apogee);                  %Starting velocity (m/s)
tstart_lateral = 0;                 %Start time (s)
dt = 1;                             %Time step (s)
tstop_lateral = 500;                %End time (s)
t_lateral = tstart_lateral:dt:tstop_lateral;
deltav_lander_lateral = Ve.*log(m0_PDI./(m0_PDI - dmdt.*t_lateral));
disp(max(deltav_lander_lateral))        %DeltaV of the LM
    
m_at_t_lateral = m0_PDI.*exp(-deltav_lander_lateral./Ve);
                                    %Mass at time t of the LM
decceleration_lateral = Thrust1./m_at_t_lateral;

lateral_velocity = V_apogee - cumsum(decceleration_lateral);
lateral_velocity(lateral_velocity <0) = 0;
s_lateral = cumsum(lateral_velocity)/1000;


 figure(2)
 plot(t_lateral,m_at_t_lateral)
 xlabel('Time (s)'), ylabel('Mass of LM over time')
 figure(3)
 plot(t_lateral,decceleration_lateral)
 xlabel('Time (s)'), ylabel('Decelleration of LM (m/s^2)')
 figure(4)
 plot(t_lateral,lateral_velocity)
 xlabel('Time (s)'), ylabel('Lateral velocity of LM (m/s)')
 figure(5)
 plot(t_lateral,s_lateral)
 xlabel('Time (s)'), ylabel('Distance travelled along surface (km)')

