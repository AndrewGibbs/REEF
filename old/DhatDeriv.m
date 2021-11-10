function vals = DhatDeriv(theta, alphas, FF, p, derivs)
% computes the n'th derivative of \hat{D}(\theta,\alphas) w.r.t \theta.
% returns a matrix of the n'th derivative, rows of theta columns of alpha

    vals = zeros(length(theta),length(alphas));
    %following sum is based on Leibnitz's generalised product rule:
    for dn_=0:derivs
        vals = vals + nchoosek(derivs,dn_).*hatDeriv(theta,dn_).*(FFDeriv(theta,derivs-dn_));
    end
      
    % derivative of hat function
    function vals = hatDeriv(theta,dn)
        if dn== 0
            vals = cos(p*repmat(theta,1,length(alphas)))-(-1)^p*cos(p*(repmat(alphas.',length(theta),1)));
        else
            vals = p^dn*cos(p*(dn*pi/2 + repmat(theta,1,length(alphas))));
        end
    end

    % derivative of far-field pattern:
    function vals = FFDeriv(theta,dn)
        vals = zeros(length(theta),length(alphas));
        for m = 1:length(alphas)
            vals(:,m) = FF{m,dn+1}(theta);
        end 
    end
end

