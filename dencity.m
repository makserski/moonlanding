function rho = dencity(h)
if h < 11019.13
    temp=59-0.00356*h;
    preassure=2116*((temp+459.7)/518.6)^5.256;
elseif h < 25098
    temp = -70;
    preassure = 473.1*exp(1.73-0.000048*h);
elseif h > 25098
    temp = -205.05 + 0.00164*h;
    preassure = 51.97*((temp+459.7)/389.98)^-11.388;
end
rho = preassure/(1718*(temp+459.7));
end

      
