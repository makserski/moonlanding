clc
clear
v=0; %initial velocity
hh=0; %displacement (not really needed)
h=0; %initial height
tstart=0; %time 0
dt=0.1; %timestep
t1=160; %time of stage 1 burn
t2 = 460; % time of stage 2 burn
g= 9.81; % g
F1 = 34000000; % stage 1 thrust
F2 = 4900000; %stage 2 thrust 
M = 2970000; % initial mass of ship
matrix=zeros(4600,4); % create the matrix that holds all the data
massflow1 = (2100000/160)*dt; % mass flow rate stage 1
massflow2 = (444000/360)*dt; % mass flow rate stage  2
s=0 % variable to help index the matrix
drag = 0.25*pi*5^2*0.05*v^2/2
angle_turn = 1.3*(90/460)*dt % the angle that the ship turns at
angle = 0
Mearth = 5.97219*10^24
Radius_orbit = 185000
G = 6.673*10^-11
speed_for_orbit = sqrt((Mearth*G)/Radius_orbit) % speed needed for orbit 
orbit_dist = 185000
for t = tstart:dt:t1
    a = getacceleration(F1,M,g,angle);
    h= h + dt*v*cosd(angle); %calculating height and dispalcemt 
    hh=hh+dt*v*sind(angle);
    v=v+dt*a;
    M = M - massflow1;
    matrix(1+s,1) = a; % these just update the matrix
    matrix(1+s,2) = v;
    matrix(1+s,3) = h;
    matrix(1+s,4) = angle;
    matrix(1+s,5) = hh;
    s=1+s;
    angle = angle + angle_turn;
                            % if the height is the target orbit height
                            % then turn the rocket paralell to the ground
    if h > orbit_dist; 
        angle = 90
        angle_turn = 0
    end
    if cosd(angle) < 0 % stop rotating the rocket once it is paralell
        angle_turn = 0
    end
    if v > speed_for_orbit % once the speed is the orbit soeed then stop accelerating
        F1 = 0
    end
end
for t = t1:dt:t2
     a = getacceleration(F2,M,g,angle); % this looop is the same as the firist one bit for the secind stage 
    h= h + dt*v*cosd(angle);
    hh=hh+dt*v*sind(angle);
    v=v+dt*a;
    M = M - massflow2;
    matrix(1+s,1) = a;
    matrix(1+s,2) = v;
    matrix(1+s,3) = h;
    matrix(1+s,4) = angle;
    matrix(1+s,5) = hh;
    s=1+s;
    angle = angle + angle_turn;
    if h > orbit_dist;
        angle = 90
        angle_turn = 0
    end
    if cosd(angle) < 0
        angle_turn = 0
    end
    if v > speed_for_orbit
        F2 = 0
    end
end
figure (1)
plot(matrix(:,1))
title('acceleration')
xlabel('time (ms)')
ylabel('acceleration')
figure (2)
plot(matrix(:,2))
title('velocity time')
xlabel('time (ms)')
ylabel('velocity')
figure (3)
plot(matrix(:,3))
title('height time')
xlabel('time (ms)')
ylabel('height ')
figure (4)
plot(matrix(:,4))
title('angle of the rocket')
xlabel('time (ms)')
ylabel('angel ')
figure (5)
plot(matrix(:,5))
title('displacement time')
xlabel('time (ms)')
ylabel('displacement')
