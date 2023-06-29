function Dhat_dn = get_Dhat_derivs(theta,alpha_,dn,p,FF,FF_derivs)

    dhat=@(theta,alpha,dn) real((1i*p)^dn*exp(1i*p*repmat(theta,1,length(alpha))))-(dn==0)*(-1)^p*cos(p*(repmat(alpha,length(theta),1)));
  
    %compute generalised Leibnitz product rule over all hats and far fields
    Dhat_dn=0;
    for n_=0:dn
        if dn-n_ == 0
            Dhat_dn = Dhat_dn  +  nchoosek(dn,n_).* dhat(theta,alpha_,n_).*FF(theta);
        else
            Dhat_dn = Dhat_dn  +  nchoosek(dn,n_).* dhat(theta,alpha_,n_).*FF_derivs{dn-n_}(theta);
        end
    end

end