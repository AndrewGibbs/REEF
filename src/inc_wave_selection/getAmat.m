function A = getAmat(alpha_in,N)
%computes the matrix A for a regular N-gon. Should extend this for
%irregular polygons too, replacing N by a vector of interior angles
    M_large = length(alpha_in);
    M = get_M_for_Ngon(N);
    A = zeros(M,M_large);
    if N <= 2
        interior_angle = 0;
    else
        interior_angle = pi*(N-2)/N;
    end
    
    %interior_angle = 2*pi-twoPhi;
    twoPhi = 2*pi-interior_angle;
        
    for m=1:M_large
        rotation = 0;
        for corner=0:(N-1)
            alpha_rot = mod(alpha_in(m) - rotation, 2*pi);
            if (alpha_rot <=pi) || (alpha_rot >= mod(pi-twoPhi,2*pi))
%                 (mod(corner*interior_angle+pi,2*pi) <= alpha_rot) ...
%                     && (alpha_rot<mod((corner+1)*interior_angle+pi,2*pi))
%                 2+2;
%             else
                for n=1:(M/N)
                    A(n+corner*M/N,m) = sin(n*pi/twoPhi*(alpha_rot - pi + twoPhi));
                end
            end
            rotation = mod(rotation + (pi - interior_angle), 2*pi);
        end
    end
end

