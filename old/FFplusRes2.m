function Dout = FFplusRes2(theta,D,B,p,Z,alpha_in,alpha_out)
% proof of concept for the residue-based approach to contour integrals

    for n=1:length(theta)
        [z1(n), z2(n)] = find_closest_two(theta(n));
    end

    hat=@(theta,alpha) cos(p*repmat(theta,1,length(alpha)))-(-1)^p*cos(p*(repmat(alpha,length(theta),1)));
    
    Dout = zeros(size(theta));
    M = length(alpha_in);
    
    a = @(chi) -1/(p*sin(p*chi));
    
    n = 0;
    for t=(theta.')
        n = n + 1;
        for m=1:M
            Dout(n) = Dout(n) + ...
                            B(m)*(...
                                    hat(t,alpha_in(m))*D{m}(t)*(z1(n)-t)*(z2(n)-t)...
                                    +hat(z1(n),alpha_in(m))*D{m}(z1(n))*hat(t,alpha_out)*(z2(n)-t)*a(z1(n))...
                                    +hat(z2(n),alpha_in(m))*D{m}(z2(n))*hat(t,alpha_out)*(z1(n)-t)*a(z2(n))...
                                  )...
                                /(hat(t,alpha_out)*(z1(n)-t)*(z2(n)-t));
        end
    end

    function [z1, z2] = find_closest_two(t)
        z1 = NaN;
        z2 = NaN;
        d1 = inf;
        d2 = inf;
        for z = Z
            d = abs(exp(1i*t)-exp(1i*z));
            if d< d1
                %demote z1 and replace z2
                d2 = d1;
                z2 = z1;
                %replace z1
                z1 = z;
                d1 = d;
            elseif d< d2
                %replace z2
                d2 = d;
                z2 = z;
            end
        end
    end
end