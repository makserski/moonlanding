function moon_launch
clc
clear

%%Initial conditions
totalT=0;
v=0; %initial velocity

hh=0; %displacement (not really needed)

h=1737400; %initial height

tstart=0; %time 0

dt=0.1; %timestep

t1=650; %time of stage 1 burn

g= 1.62; % g

F = 16000; % stage 1 thrust

M = 4700; % initial mass of ship

matrixmoon=zeros(); % create the matrix that holds all the data

massflow1 = dt*2250/615 ;
s=0; % variable to help index the matrix

drag = 0;

angle_turn = 1.5*(90/650)*dt; % the angle that the ship turns at

angle = 0;

Mmoon = 7.342*10^(22);

orbit_dist = 1848220;

G = 6.673*10^-11;

rho = 0;

speed_for_orbit = sqrt((Mmoon*G)/(orbit_dist)); % speed needed for orbit 

for t = tstart:dt:t1

    a = getacceleration(F,M,g,angle,v,rho);

    h = getheight(h,dt,v,angle); %calculating height and dispalcemt 

    hh=hh+dt*v*sind(angle);

    v=v+dt*a;

    M = M - massflow1;
    matrixmoon = updatetable(s,a,v,h,angle,hh,totalT,matrixmoon);

    totalT=totalT+dt;
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

plot(matrixmoon(:,6),matrixmoon(:,1))

title('acceleration')

xlabel('time (s)')

ylabel('acceleration')

figure (2)

plot(matrixmoon(:,6),matrixmoon(:,2))

title('velocity time')

xlabel('time (s)')

ylabel('velocity')

figure (3)

plot(matrixmoon(:,6),matrixmoon(:,3))

title('height time')

xlabel('time (s)')

ylabel('height ')

figure (4)

plot(matrixmoon(:,6),matrixmoon(:,4))

title('angle of the rocket')

xlabel('time (s)')

ylabel('angel ')

figure (5)

plot(matrixmoon(:,6),matrixmoon(:,5))

title('displacement time')

xlabel('time (s)')

ylabel('displacement')
end
