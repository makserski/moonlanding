    th = linspace(0,2*pi) ; 
R = (6564.03e3) ;
x = R*cos(th-pi/2+pi/32) ; 
y = -R*sin(th-pi/2+pi/32) ;
comet(x,y) ;
