%The following code shows the translunar injection. It shows how vast the
%distance to the moon from Earth really is, since it is to scale. 

clc
clear all

[x,y,z] = sphere;               %This plot creates a unit sphere centered at (0,0,0)
SC = 6371e3                     %SC is the scaling factor to make a realistic sized earth
figure(1)
surf(SC.*x,SC.*y,SC.*z)         %This produces a surface plot of the sphere
axis equal                      %This sets the axes to be equally separated
hold on                         %This stops the plot from finishing to allow a second plot to be overlayed
th = linspace(0,2*pi) ;         %This produces an array of angle sizes for the circular orbit
R = (6569.3e3) ;                %Radius of the earth
x = R*cos(th+ 3*pi/4) ;         %Polar coordinate plot to produce an orbit plot
y = R*sin(th+ 3*pi/4) ;
comet(x,y) ;                    %Comet plot shows the progression of the translunar injection


EARTH_ORBIT = 6569.3e3;         % Parking orbit of the Earth
MOON_ORBIT = 1848e3;            % Orbit of the Moon
EARTH_MOON_DIST = 384400e3;     % Distance from the Earth to the Moon
X_COORD(1) = -4645e3;           % Starting point at x=0
Y_COORD(1) = 4645e3;            % Starting point at y=0


% Initialise the time array

dt = 1;
TIME = 0:dt:10300;

% Calculating the geometrical properties of the elliptical orbit during
% the Hohmann transfer to the Moon


% Properties of the elliptical orbit during the Hohman transfer from
% Earth's parking orbit to the large orbit which coincides with Moon's
% orbit:

SEMI_MAJOR_AXIS = 194560e3;     % (a)
SEMI_MINOR_AXIS = 50129.61e3;   % (b)
ECCENTRICITY = 0.966237;        % Calculated from geometry (e)
ASCENDING_NODE = 1.12;          % Longitude of ascending node (OMEGA)
ARG_PERICENTRE = 1.23;          % Argument of pericentre (omega)
INCLINATION = 0;                % Angle between the equator and the orbit (i)

% True anomaly (v) also needs to be found but its value changes when the
% location of the rocket changes so its values will be found within the
% for-loop below

dE = pi/10300;
E = 0:dE:pi;                    % This assumes a constant angular velocity 
                                % about the centre point of the big circle

for i = 1:length(TIME)
    
    TRUE_ANOMALY(i) = 2*atan(sqrt((1+ECCENTRICITY)/(1-ECCENTRICITY)).*tan((E(i))/2));
    r(i) = (SEMI_MAJOR_AXIS*(1+ECCENTRICITY^2))./(1+ECCENTRICITY.*cos(TRUE_ANOMALY(i)));
    
    X_COORD(i+1) = r(i).*(cos(TRUE_ANOMALY(i) + ARG_PERICENTRE)*cos(ASCENDING_NODE) - ...
        sin(TRUE_ANOMALY(i) + ARG_PERICENTRE)*sin(ASCENDING_NODE)*cos(INCLINATION));
    
    Y_COORD(i+1) = r(i).*(cos(TRUE_ANOMALY(i) + ARG_PERICENTRE)*sin(ASCENDING_NODE) + ...
        sin(TRUE_ANOMALY(i) + ARG_PERICENTRE)*cos(ASCENDING_NODE)*cos(INCLINATION));
    
    %This code produces a set of X and Y coordinates, which can be plotted
    %on the figure. This code produces a Hohmann transfer, which is a large
    %elliptical orbit, bringing to where the Command Module will go into a
    %parking orbit around the moon.
    

    
end

%The following code removes an initial error that creates a straight line which
%throws the injection plot off. Removing it leaves the correct ellipse.

INITIAL_ERROR_X = X_COORD(2)-X_COORD(1);

for i = 1:length(X_COORD)-1
    X_COORD(i+1) = X_COORD(i+1) - INITIAL_ERROR_X;
end

INITIAL_ERROR_Y = Y_COORD(2)-Y_COORD(1);

for i = 1:length(Y_COORD)-1
    Y_COORD(i+1) = Y_COORD(i+1) - INITIAL_ERROR_Y;
end

comet(X_COORD,Y_COORD);
hold off                        %This allows all of the plots to finally be plotted