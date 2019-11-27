function a = getacceleration(F,M,g,angle,v,rho)
    drag = (0.25*pi*5^2*rho*v^2);
    a = F/M - g*cosd(angle);
end 