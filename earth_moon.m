clc
clear all
[x,y,z] = sphere;
figure(1)
surf(6371e3.*x,6371e3.*y,36371e3.*z)
axis equal
hold on 

% surf(x+100,y,z) % centered at (0,1,-3)
     th = linspace(0,2*pi) ; 
 R = (6569.3e3) ;
 x = R*cos(th+ 3*pi/4) ; 
 y = R*sin(th+ 3*pi/4) ; 
 comet(x,y) ;
% t = linspace(0,2*pi,200) ; 
% 
% xc = [0:1:100] 
% yc = 0.15.*xc -10
% comet(xc,yc)
% 
% R1 = (1)*5 ;
% xb = R1*cos(t+(3*pi/2)) + 100 ; 
% yb = -R1*sin(t+(3*pi/2)) ; 
% comet(xb,yb) ;
% 
% % compute points corresponding to axis-oriented ellipse
% r1 = (1.1)*5;
% r2 = (0.9)*5;
% xt = r1 * cos(t);
% yt = -r2 * sin(t);
% % aply rotation by angle theta
%  %angle = pi/2;
%  %cot = cos(angle); sit = sin(angle);
%  %xa = xt * cot - yt * sit + 100;
%  %ya = xt * sit - yt * cot;
%  %draw the curbe
% %comet(xa, ya);
% xd = [100:-1:0] 
% yd = -0.15.*xd + 10
% comet(xd,yd)

 EARTH_ORBIT = 6569.3e3; % Parking orbit of the Earth
    MOON_ORBIT = 1848e3; % Orbit of the Moon
    EARTH_MOON_DIST = 384400e3; % Distance from the Earth to the Moon
    X_COORD(1) = -4645e3; % Starting point at x=0
    Y_COORD(1) = 4645e3; % Starting point at y=0
    %G = 6.673e-11;
    
    % Initialise the time array
    dt = 1;
    TIME = 0:dt:10300;
    
    % Calculating the geometrical properties of the elliptical orbit during
    % the Hohmann transfer to the Moon
    %
    % here
    
    % Properties of the elliptical orbit during the Hohman transfer from
    % Earth's parking orbit to the large orbit which coincides with Moon's
    % orbit:
    SEMI_MAJOR_AXIS = 194560e3; % (a)
    SEMI_MINOR_AXIS = 50129.61e3; % (b)
    ECCENTRICITY = 0.966237; % Calculated from geometry (e)
    ASCENDING_NODE = 1.12; % Longitude of ascending node (OMEGA)
    ARG_PERICENTRE = 1.23; % Argument of pericentre (omega)
    INCLINATION = 0; % Angle between the equator and the orbit (i)
    % True anomaly (v) also needs to be found but its value changes when the
    % location of the rocket changes so its values will be found within the
    % for-loop below
    
    dE = pi/10300;
    E = 0:dE:pi; % This assumes a constant angular velocity about the centre point of the big circle

    %MASS_EARTH = 5.9724e24;
    
    for i = 1:length(TIME)

        TRUE_ANOMALY(i) = 2*atan(sqrt((1+ECCENTRICITY)/(1-ECCENTRICITY)).*tan((E(i))/2));
        r(i) = (SEMI_MAJOR_AXIS*(1+ECCENTRICITY^2))./(1+ECCENTRICITY.*cos(TRUE_ANOMALY(i)));
        
        X_COORD(i+1) = r(i).*(cos(TRUE_ANOMALY(i) + ARG_PERICENTRE)*cos(ASCENDING_NODE) - ...
        sin(TRUE_ANOMALY(i) + ARG_PERICENTRE)*sin(ASCENDING_NODE)*cos(INCLINATION));
        
        Y_COORD(i+1) = r(i).*(cos(TRUE_ANOMALY(i) + ARG_PERICENTRE)*sin(ASCENDING_NODE) + ...
        sin(TRUE_ANOMALY(i) + ARG_PERICENTRE)*cos(ASCENDING_NODE)*cos(INCLINATION));
        
        %add comments here
    
        %VELOCITY(i) = sqrt((MASS_EARTH*G).*((2/r(i))-(1/SEMI_MAJOR_AXIS)));
        
    end
    
    INITIAL_ERROR_X = X_COORD(2)-X_COORD(1);
    
    for i = 1:length(X_COORD)-1
        X_COORD(i+1) = X_COORD(i+1) - INITIAL_ERROR_X;
    end
    
    INITIAL_ERROR_Y = Y_COORD(2)-Y_COORD(1);
    
    for i = 1:length(Y_COORD)-1
        Y_COORD(i+1) = Y_COORD(i+1) - INITIAL_ERROR_Y;
    end
    
    comet(X_COORD,Y_COORD);
hold off