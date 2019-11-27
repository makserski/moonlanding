function [X_COORD, Y_COORD, VELOCITY, TIME] = translunar_injection(z)

    %clear all
    
    EARTH_ORBIT = 6569.3e3; % Parking orbit of the Earth
    MOON_ORBIT = 1848e3; % Orbit of the Moon
    EARTH_MOON_DIST = 384400e3; % Distance from the Earth to the Moon
    X_COORD(1) = 0; % Starting point at x=0
    Y_COORD(1) = 0; % Starting point at y=0
    G = 6.673e-11; % Gravitational constant
    MASS_EARTH = 5.9724e24; % Mass of the Earth
    MU = G*MASS_EARTH; % Standard gravitational parameter
    
    % Initialise the time array
    dt = 20;
    TIME = 0:dt:262763;
    
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
    
    dE = pi/length(TIME);
    E = 0:dE:pi; % This assumes a constant angular velocity about the centre point of the big circle

    for i = 1:length(TIME)
        
        if i == length(TIME) % Making sure the number of array elements does not exceed length(TIME)
            TRUE_ANOMALY(i) = 2*atan(sqrt((1+ECCENTRICITY)/(1-ECCENTRICITY))*tan((E(i))/2));
            r(i) = (SEMI_MAJOR_AXIS*(1+ECCENTRICITY^2))/(1+ECCENTRICITY*cos(TRUE_ANOMALY(i)));
        
        else % The normal case
            TRUE_ANOMALY(i) = 2*atan(sqrt((1+ECCENTRICITY)/(1-ECCENTRICITY))*tan((E(i))/2));
            r(i) = (SEMI_MAJOR_AXIS*(1+ECCENTRICITY^2))/(1+ECCENTRICITY*cos(TRUE_ANOMALY(i)));

            X_COORD(i+1) = r(i)*(cos(TRUE_ANOMALY(i) + ARG_PERICENTRE)*cos(ASCENDING_NODE) - ...
            sin(TRUE_ANOMALY(i) + ARG_PERICENTRE)*sin(ASCENDING_NODE)*cos(INCLINATION));

            Y_COORD(i+1) = r(i)*(cos(TRUE_ANOMALY(i) + ARG_PERICENTRE)*sin(ASCENDING_NODE) + ...
            sin(TRUE_ANOMALY(i) + ARG_PERICENTRE)*cos(ASCENDING_NODE)*cos(INCLINATION));
        end
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
    
    comet(z,X_COORD,Y_COORD)
    
    % Constants and arrays required for the velocity calculation
    n = sqrt(MU/(SEMI_MAJOR_AXIS^3));
    l1 = cos(ASCENDING_NODE)*cos(ARG_PERICENTRE) - sin(ASCENDING_NODE)*sin(ARG_PERICENTRE)*cos(INCLINATION);
    l2 = - cos(ASCENDING_NODE)*sin(ARG_PERICENTRE) - sin(ASCENDING_NODE)*cos(ARG_PERICENTRE)*cos(INCLINATION);
    m1 = sin(ASCENDING_NODE)*cos(ARG_PERICENTRE) + cos(ASCENDING_NODE)*sin(ARG_PERICENTRE)*cos(INCLINATION);
    m2 = - sin(ASCENDING_NODE)*sin(ARG_PERICENTRE) - cos(ASCENDING_NODE)*cos(ARG_PERICENTRE)*cos(INCLINATION);
    n1 = sin(ARG_PERICENTRE)*sin(INCLINATION);
    n2 = cos(ARG_PERICENTRE)*sin(INCLINATION);
    dNU = pi/(length(TIME)-1); 
    NU = 0:dNU:pi;
    
    %VELOCITY(1) = 7.7889e+03;
    % Calculating the velocity of the spacecraft
    for i = 1:length(TIME)
        R(i) = sqrt((Y_COORD(i)^2)+((X_COORD(i)-EARTH_ORBIT)^2));
%         VELOCITY(i) = sqrt((MASS_EARTH*G)*((2/R(i))-(1/SEMI_MAJOR_AXIS)));
        X_VELOCITY(i) = ((n*SEMI_MAJOR_AXIS)/r(i))*(SEMI_MINOR_AXIS*l2*cos(NU(i)) - SEMI_MAJOR_AXIS*l1*sin(NU(i)));
        Y_VELOCITY(i) = ((n*SEMI_MAJOR_AXIS)/r(i))*(SEMI_MINOR_AXIS*m2*cos(NU(i)) - SEMI_MAJOR_AXIS*m1*sin(NU(i)));
        VELOCITY(i) = sqrt(X_VELOCITY(i)^2+Y_VELOCITY(i)^2);
    end
%     figure(2)
%     plot(TIME,X_VELOCITY)
    % Draw x and y which represent the path of the rocket in 2D
    
    va = sqrt(((G*MASS_EARTH)/SEMI_MAJOR_AXIS)*((1-ECCENTRICITY)/(1+ECCENTRICITY)));
    vp = sqrt(((G*MASS_EARTH)/SEMI_MAJOR_AXIS)*((1+ECCENTRICITY)/(1-ECCENTRICITY)));
end