% F2t & F2c calculator
% F in 0 degree after transformation
F = [-1.98*10^4; 6.6074*10^3; 0];
t = -30
c = cosd(t)
s = sind(t)
T = [c s 0; -s c 0; 0 0 1];
F_angle = T*F
