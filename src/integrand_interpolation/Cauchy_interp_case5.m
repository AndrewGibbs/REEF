function I = Cauchy_interp_case5(theta, theta_0, B, Dhatin, Dhatin_derivs, ...
                                hatout, p, h, Nquad, dtheta)


    a = min([theta; theta_0])-h;
    b = max([theta; theta_0])+h;

    [z,w] = Cauchy_box_quad(theta, a, b, h, Nquad, 1, true);

%     P = zeros(size(w));

    %% consider interpolating polynomial p(z;\theta,\theta_0),
    % which interpolates 'sum b Dhat' at theta and theta_0.
    % here this is represented by a matrix, as theta and z are vectors.
    P1 = zeros(length(theta),1);
    P2 = 0;
    for m=1:length(B)
        P1 = P1 + B(m)*Dhatin{m}(theta);
        P2 = P2 + B(m)*Dhatin{m}(theta_0);
    end
    P = ((z-theta_0).'.*P1 - (z.'-theta)*P2)./(theta-theta_0);

    I = (w.'.*P)*(1./hatout(z));

    %% special cases - where difference between residues is negligable
    theta2theta_0 = cyclic_distance(theta,theta_0);
    inds_5i = theta2theta_0<dtheta;
    
    % need to use second derivative (e.g. L'hopital's)
    if sum(inds_5i)>0
%         dhatout = @(theta) -p*sin(p*theta);
%         I(inds_5i) = 0;
%         for m=1:length(B)
%             I(inds_5i) = I(inds_5i) + B(m)*Dhatin_derivs{m}{1}(theta(inds_5i));
%         end
%         I(inds_5i) = I(inds_5i)./dhatout(theta(inds_5i));
        %create linear function which interpolates value and derivative at
        %theta_0
        P1_ = 0; % gradient
        for m=1:length(B)
            P1_ = P1_ + B(m)*Dhatin_derivs{m}{1}(theta_0);
        end
        P2_ = P2-P1_*theta_0; % intercept

        [z_,w_] = Cauchy_box_quad(theta_0, theta_0-h, theta_0+h, h, Nquad, 1, true);
        P_ = P1_*z + P2_;
        I(inds_5i) = (w_.'.*P_.')*(1./hatout(z_));
    end
end