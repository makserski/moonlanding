%function [X_COORD, Y_COORD] = translunar_injection()

    EARTH_ORBIT = 6569.3e3; % Parking orbit of the Earth
    MOON_ORBIT = 1848e3; % Orbit of the Moon
    EARTH_MOON_DIST = 384400e3; % Distance from the Earth to the Moon
    X_COORD(1) = 0; % Starting point at x=0
    Y_COORD(1) = 0; % Starting point at y=0
    G = 6.673e-11;
    
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

    MASS_EARTH = 5.9724e24;
    
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
    %figure(1)
    %comet(X_COORD,Y_COORD)
    
    %VELOCITY(1) = 7.7889e+03;
    for i = 1:length(TIME)-1
        R(i) = sqrt((Y_COORD(i).^2)+((X_COORD(i)-EARTH_ORBIT).^2));
        VELOCITY(i+1) = sqrt((MASS_EARTH*G).*((2/R(i))-(1/SEMI_MAJOR_AXIS)));
    end
    figure(2)
    plot(TIME,VELOCITY)
    % Draw x and y which represent the path of the rocket in 2D
%end
    va = sqrt(((G*MASS_EARTH)/SEMI_MAJOR_AXIS)*((1-ECCENTRICITY)/(1+ECCENTRICITY)))
    vp = sqrt(((G*MASS_EARTH)/SEMI_MAJOR_AXIS)*((1+ECCENTRICITY)/(1-ECCENTRICITY)))