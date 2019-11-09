function a = getacceleration(F,M,g,angle);
 a = F/M - g*cosd(angle);
end 