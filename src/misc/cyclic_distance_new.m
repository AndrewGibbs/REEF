function r = cyclic_distance(x,y)
% get angular distance on circle

% first restrict real part to [0,2\pi)
    x = mod(real(x),2*pi) + imag(x);
    y = mod(real(y),2*pi) + imag(y);
% now consider every case

    % get signed displacements
    v1 = x-y;
    v2 = x-y-2*pi;
    v3 = x-y+2*pi;

    % get distances
    r0 = abs(v1);
    r1 = abs(v2);
    r2 = abs(v3);

    % work out which way around the circle is the closest
    min_dist_inds = ones(length(x));
    min_dist_inds(r1<r0 & r1<=r2) = 2;
    min_dist_inds(r2<r0 & r1>=r2) = 3;

    % allocate min displacements
    r(min_dist_inds==1) = v1(min_dist_inds==1);
    r(min_dist_inds==2) = v2(min_dist_inds==2);
    r(min_dist_inds==3) = v3(min_dist_inds==3);

%     r = min(min(r0,r1),r2);
end