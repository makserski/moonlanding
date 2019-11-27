function launch
clc
clear

%%Initial conditions
v=0; %initial velocity
hh=0; %displacement (not really needed)
h=6371000; %initial height
tstart=0; %time 0
dt=0.1; %timestep
rho = 0;
t1=160; %time of stage 1 burn
t2 = 460; % time of stage 2 burn
t3 = 990 ; % time of stage 3 burn
g= 9.81; % g
F1 = 35100000; % stage 1 thrust
F2 = 5141000; %stage 2 thrust 
F3 = 1028200; %stage 3 thrustm
M = 2970000; % initial mass of ship

matrix=zeros(); % create the matrix that holds all the data

% mass flow rate stage 1
massflow1= 12580*dt;
% mass flow rate stage  2
massflow2 = 1137*dt;
  % mass flow rate stage  3
massflow3 = 213*dt;
s=0 ;% variable to help index the matrix
drag = 0;
angle_turn = 1.4*(90/460)*dt; % the angle that the ship turns at
angle = 0;
Mearth = 5.97219*10^24;
orbit_dist = 6569300;
G = 6.673*10^-11;
di = 51.57;
speed_for_orbit = sqrt((Mearth*G)/(orbit_dist)); % speed needed for orbit 
totalT =0;

for t = tstart:dt:t1
    a = getacceleration(F1,M,g,angle,v,rho);
    h = getheight(h,dt,v,angle); %calculating height and dispalcemt 
    rho = density(h);
    hh=hh+dt*v*sind(angle);
    v=v+dt*a;
    M = M - massflow1;
    matrix = updatetable(s,a,v,h,angle,hh,totalT,matrix,M);
    s=1+s;
    totalT=totalT+dt;
    angle = angle + angle_turn;

                            % if the height is the target orbit height
                            % then turn the rocket paralell to the ground

    if h > orbit_dist
        angle = 90;
        angle_turn = 0;
    end

    if cosd(angle) < 0 % stop rotating the rocket once it is paralell
        angle_turn = 0;
    end

    if v > speed_for_orbit % once the speed is the orbit soeed then stop accelerating
        F1 = 0;
    end
end
M = M - 200000; %loss of first stage

for t = t1:dt:t2
    a = getacceleration(F2,M,g,angle,v,rho); % this loop is the same as the firist one bit for the secind stage 
    h = getheight(h,dt,v,angle);
    rho = density(h);
    hh = hh+dt*v*sind(angle);
    v=v+dt*a;
    M = M - massflow2;
    matrix = updatetable(s,a,v,h,angle,hh,totalT,matrix,M);
    totalT=totalT+dt;
    s=1+s;
    angle = angle + angle_turn;

    if h > orbit_dist
        angle = 90;
        angle_turn = 0;
    end

    if cosd(angle) < 0
        angle_turn = 0;
    end

    if v > speed_for_orbit
        F2 = 0;
    end
end

M=M-160000;

for t = t2:dt:t3
    a = getacceleration(F3,M,g,angle,v,rho); % this looop is the same as the firist one bit for the secind stage 
    h = getheight(h,dt,v,angle);
    rho = density(h);
    hh=hh+dt*v*sind(angle);
    v=v+dt*a;
    M = M - massflow3;
    matrix = updatetable(s,a,v,h,angle,hh,totalT,matrix,M);

    s=1+s;
    totalT=totalT+dt;
    
    if h > orbit_dist
        angle = 90;
        angle_turn = 0;
    end

    if cosd(angle) < 0
        angle_turn = 0;
    end

    if v > speed_for_orbit
        F3 = 0;
    end
end

dv = 2*v*sin(di/2); % velocity change to alter incliantion or orbit
fprintf('to change the inclination of the orbit to one around the equator a velocity change of  %d m/s must be made towards the equator\n', dv);
figure (1)
plot(matrix(:,6),matrix(:,1))
title('acceleration')
xlabel('time (s)')
ylabel('acceleration')

figure (2)
plot(matrix(:,6),matrix(:,2))
title('velocity time')
xlabel('time (s)')
ylabel('velocity')

figure (3)
plot(matrix(:,6),matrix(:,3))
title('height time')
xlabel('time (s)')
ylabel('height ')

figure (4)
plot(matrix(:,6),matrix(:,4))
title('angle of the rocket')
xlabel('time (s)')
ylabel('angel ')

figure (5)
plot(matrix(:,6),matrix(:,5))
title('displacement time')
xlabel('time (s)')
ylabel('displacement')

figure(6)
plot(matrix(:,5),matrix(:,3))
xlabel('disp')
ylabel('height')

figure (8)
plot(matrix(:,6),matrix(:,7));
xlabel('Time (s)')
ylabel('Mass (kg)')
end





