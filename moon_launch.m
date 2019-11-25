clc
clear

%%Initial conditions

v=0; %initial velocity

hh=0; %displacement (not really needed)

h=0; %initial height

tstart=0; %time 0

dt=0.1; %timestep

t1=650; %time of stage 1 burn

g= 1.62; % g

F = 16000; % stage 1 thrust

M = 4700; % initial mass of ship

matrix=zeros(4600,4); % create the matrix that holds all the data

massflow1 = 2250/t1
s=0 % variable to help index the matrix

drag = 0

angle_turn = 0.98*(90/435)*dt % the angle that the ship turns at

angle = 0

Mmoon = 7.342*10^(22)

orbit_dist = 90000

G = 6.673*10^-11

rho = 0

speed_for_orbit = sqrt((Mmoon*G)/(orbit_dist+1737400)) % speed needed for orbit 

for t = tstart:dt:t1

    a = getacceleration(F,M,g,angle,v,rho);

    h = getheight(h,dt,v,angle); %calculating height and dispalcemt 

    hh=hh+dt*v*sind(angle);

    v=v+dt*a;

   % M = M - massflow1;

    matrix(1+s,1) = a; % these just update the matrix

    matrix(1+s,2) = v;

    matrix(1+s,3) = h;

    matrix(1+s,4) = angle;

    matrix(1+s,5) = hh;

    s=1+s;

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

        F = 0;

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
