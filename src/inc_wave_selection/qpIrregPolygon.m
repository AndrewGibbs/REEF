function [ q_best,p_best ] = qpIrregPolygon( angle, pMax, qMax )
%approximate the internal angle by a rational number
    error_best=inf;
    for q=1:qMax
        for p=(q+1):pMax %internal angle, so p>q
            qp_error=abs(angle-pi*q/p);
            if qp_error<error_best
               q_best=q; p_best=p; 
               error_best=qp_error;
            end
        end
    end
    
    %incase error is somehow less in scaled up q/p, find lowest common
    %factors:
    GCD=gcd(p_best,q_best);
    while GCD>1
        p_best=p_best/GCD;
        q_best=q_best/GCD;
        GCD=gcd(p_best,q_best);
    end
                
    if error_best>10^(-10)
        warning('Error in rational angle approximation by q/p is quite large, consider choosing larger values of qmax and pmax');
    end
end