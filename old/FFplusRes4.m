function Dout = FFplusRes4(theta, D, B, p, Z, alpha_in, alpha_out)
    Z = sort(Z);
    h = 1e-3;
    qppw = 100;
    kwave = 1;
    
    ResIndices1 = zeros(length(alpha_out), length(theta));
    ResIndices2 = zeros(length(alpha_out), length(theta));
    
    for n = 1:length(alpha_out)
        [res_index1, res_index2] = split_obs_angles(theta, Z);
        ResIndices1(n,:) = res_index1;
        ResIndices2(n,:) = res_index2;
    end
    ResVals1 = Z(ResIndices1);
    ResVals2 = Z(ResIndices2);

    hat=@(theta,alpha) cos(p*repmat(theta,1,length(alpha)))-(-1)^p*cos(p*(repmat(alpha,length(theta),1)));
    
    % store the residue values so they are not recomputed loads of times
    Dhat_res_vals = zeros(length(alpha_out),length(Z));
    for m =1:length(alpha_in)
        Dhat_res_vals(m,:) = hat(Z.',alpha_in(m)).*D{m}(Z);
    end
    %these distances will be wrong over the 2pi-0 jump
    ResDist1 = (ResVals1.'-theta);
    ResDist2 = (ResVals2.'-theta);
    
    a = @(chi) -1./(p*sin(p*chi));
    
    top = zeros(length(alpha_out), length(theta));
    bottom = (hat(theta,alpha_out).*ResDist1.*ResDist2).';
    for m=1:length(alpha_in)
        top = top + B(m)*(...
                        ResDist1.*ResDist2.*hat(theta,alpha_in(m)).*D{m}(theta) ...
                        + a(ResVals1.').*ResDist2.*(Dhat_res_vals(m,ResIndices1)).' ...
                        + a(ResVals2.').*ResDist1.*(Dhat_res_vals(m,ResIndices2)).'...
                    ).';
    end
    
    Dout = top./bottom;
    
    % now replace the dodgy bits which may be wrong due to rounding errors
    % need to make the roundL and roundR vectors, and account for both
    % cases
    
    % make the limits for the endpoints of the regions which we integrate
    % around to avoid rounding errors
    
    if abs(Z(2)-Z(1))<2*h %case 2a
        roundL = Z(1:2:end) - 2*h;
        roundR = Z(2:2:end) + 2*h;
    elseif abs(Z(3)-Z(2))<2*h %case 2b
        roundL = Z(2:2:end) - 2*h;
        roundL = [0 roundL Z(end)];
        roundR = Z(3:2:end) + 2*h;
        roundR = [Z(1) roundR 2*pi];
    else %case 1 (most common)
        roundL = Z - 2*h;
        roundR = Z + 2*h;
    end
    
    for n=1:length(alpha_out)
        for q =1:length(roundL)
            theta_dodgy_inds = (roundL(q)-h<theta & theta<roundR(q)+h);
            theta_dodgy = theta(theta_dodgy_inds);
            [z,w] = Cauchy_box_quad(theta_dodgy,roundL(q)-2*h,roundR(q)+2*h,h,qppw,kwave);
            replacement = zeros(length(theta_dodgy),1);
            for m=1:length(alpha_in)
                replacement = replacement + B(m)*w.'*(hat(z,alpha_in(m)).*D{m}(z)./hat(z,alpha_out(n)));
            end
            Dout(n,theta_dodgy_inds) = replacement;
        end
    end
end