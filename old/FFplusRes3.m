function Dout = FFplusRes3(theta,D,B,p,Z,alpha_in,alpha_out)
% proof of concept for the residue-based approach to contour integrals
% next thing to change - group the divisions by the Laurent coefficients

    for n=1:length(theta)
        [z1(n), z2(n)] = find_closest_two(theta(n));
    end

    hat=@(theta,alpha) cos(p*repmat(theta,1,length(alpha)))-(-1)^p*cos(p*(repmat(alpha,length(theta),1)));
    
    Dout = zeros(size(theta));
    M = length(alpha_in);
    
    a = @(chi) -1/(p*sin(p*chi)); % =-1/(p*sin(alpha_out))
    %ainv = @(chi) -(p*sin(p*chi));
    ainv = - (p*sin(alpha_out));
    
    n = 0;
    for t=(theta.')
        n = n + 1;
        Dsum = 0;
%         a_1 = a(z1(n));s
%         a_2 = a(z2(n));
        a_1 = -1/(p*sin(alpha_out));
        a_2 = -1/(p*sin(-alpha_out));
        for m=1:M
            Dsum = Dsum +  B(m)*(...
                                    hat(t,alpha_in(m))*D{m}(t)*(z1(n)-t)*(z2(n)-t)*ainv...
                                    +hat(z1(n),alpha_in(m))*D{m}(z1(n))*hat(t,alpha_out)*(z2(n)-t)...
                                    -hat(z2(n),alpha_in(m))*D{m}(z2(n))*hat(t,alpha_out)*(z1(n)-t)...
                                  );
                                    %may have messed up the signs above and below here
                                
        end
        Dout(n) = Dout(n) + Dsum/(ainv*hat(t,alpha_out)*(z1(n)-t)*(z2(n)-t));
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
        Z = sort([z1 z2]);
        z1 = Z(1);
        z2 = Z(2);
    end
end