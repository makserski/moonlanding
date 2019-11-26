clc
clear                                   % works at F1 =15200 
                                        % F2=26500

h = 15000 %initial height
g=1.62 %gravity
M = 15200 %inital mass
F1=14100 %thrust 1
F2=26500 %thrust 2
matrix=zeros(1,6)
t0= 0 
dt =0.1 %timestep
massflow=14.76*dt 
time1 = 150 % time after first stage
time2= 300 %landing time
angle = 0 %needed to run acceleration func
v = 0 %inital velocity
rho=0 %needed for acceleratiom
s=0 %iteration
hh=0 %needed
totalT=0 %time (neede)
for t = t0:dt:time1
    a = getacceleration(F1,M,g,angle,v,rho);
    h = getheight(h,dt,v,angle);
    v=v+dt*a;
    M = M - massflow;
    matrix = updatetable(s,a,v,h,angle,hh,totalT,matrix);
    s=1+s; %iteration counter for mattrix 
    totalT=totalT+dt; %x axis (stops it plotting in ms) (look at the plot functions)
end

for t = time1:dt:time2
    a = getacceleration(F2,M,g,angle,v,rho);
    h = getheight(h,dt,v,angle);
    v=v+dt*a;
    M = M - massflow;
    matrix = updatetable(s,a,v,h,angle,hh,totalT,matrix);
    s=1+s;
    totalT=totalT+dt;
end

figure (1)

plot(matrix(:,6),matrix(:,1)) %plot time acceleration

title('acceleration')

xlabel('time (ms)')

ylabel('acceleration')

figure (2)

plot(matrix(:,6),matrix(:,2)) %plot time velocity

title('velocity time')

xlabel('time (ms)')

ylabel('velocity')

figure (3)

plot(matrix(:,6),matrix(:,3)) %plot height time 

title('height time')

xlabel('time (s)')

ylabel('height ')
