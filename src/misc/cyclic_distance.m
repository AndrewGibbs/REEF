function r = cyclic_distance(x,y)
% get angular distance on circle

% first restrict real part to [0,2\pi)
    x = mod(real(x),2*pi) + imag(x);
    y = mod(real(y),2*pi) + imag(y);
% now consider every case
    r0 = abs(x-y);
    r1 = abs(x-y-2*pi);
    r2 = abs(x-y+2*pi);
    r = min(min(r0,r1),r2);
end