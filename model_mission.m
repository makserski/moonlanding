%This code is a basic representation of the Apollo 11 mission, showing
%important stages. It shows the orbit of the Earth and the Moon, and the
%Translunar injection and return journey. This model is not to scale,
%since it is just meant to be a rough representation of the mission.
clc 
clear all

[x,y,z] = sphere;
figure(1)
surf(3.67.*x,3.67.*y,3.67.*z)           %This produces the model sphere for the Earth
axis equal                              %This command keeps the axes separation equal
hold on                                 %Allows further plots to be added

surf(x+100,y,z)                         %This plots the model sphere for the moon
th = linspace(0,2*pi) ;                 %This plots an array for the angle
R = (10) ;                              %This produces the radius of the orbit around the Earth
x = R*cos(th-pi/2) ;                    %This plots the circular orbit around the Earth
y = R*sin(th-pi/2) ;
comet(x,y) ;                            %Comet plot of the orbit

xc = [0:1:100];                         %This produces an array for the Translunar injection plot
yc = 0.15.*xc -10;                      %This produces the straight line for the Translunar injection
comet(xc,yc)                            %Comet plot of the Translunar injection

t = linspace(0,2*pi) ;                  %This plots an array for the angle
R1 = (1)*5 ;                            %This is thh radius of the Lunar orbit
xb = R1*cos(t+(3*pi/2)) + 100 ;         %This plots the circular orbit around the Moon
yb = -R1*sin(t+(3*pi/2)) ;
comet(xb,yb) ;                          %Comet plot of the Lunar orbit

xd = [100:-1:0];                         %This produces an array for the return journey plot
yd = -0.15.*xd + 10;                     %This produces the straight line for the return journey
comet(xd,yd)                            %Comet plot of the return journey

hold off