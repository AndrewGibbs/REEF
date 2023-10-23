function I = Cauchy_interp_case7(theta, theta_0, theta_0_, B_coeffs, Dhatin, ...
                            Dhatin_derivs, hatout, p, h, Nquad, dtheta, sum_bm_Dhat_at_theta)

    a = min([theta; theta_0; theta_0_])-h;
    b = max([theta; theta_0; theta_0_])+h;

    [z,w] = Cauchy_box_quad(theta, a, b, h, Nquad, 1, true);
%     z = linspace(a,b,1000).';

%     P = zeros(size(w));

    %% consider interpolating polynomial p(z;\theta,\theta_0,\theta_0'),
    % which interpolates 'sum b Dhat' at theta and theta_0.
    % here this is represented by a matrix, as theta and z are vectors.

    % want this to be barycentric style, so construct the weights first:

    % potential for speed up here - repeated values
    lambda_1 = 1./((theta-theta_0).*(theta-theta_0_));
    lambda_2 = 1./((theta_0-theta).*(theta_0-theta_0_));
    lambda_3 = 1./((theta_0_-theta).*(theta_0_-theta_0));
    % further speedup possible - could multiply each by its 1/(z-?) bit
    % here

    % construct the integrand-inerpolating polynomial, with three terms
%     P1 = zeros(length(theta),1);
    P1 = sum_bm_Dhat_at_theta;
    P2 = 0;
    P3 = 0;
    for m=1:length(B_coeffs)
%         P1 = P1 + B_coeffs(m)*Dhatin{m}(theta);
        % possible speedup option - have these already been computed?
        P2 = P2 + B_coeffs(m)*Dhatin{m}(theta_0);
        P3 = P3 + B_coeffs(m)*Dhatin{m}(theta_0_);
    end


    %changed the below as it wasn't vectorised
%     P = ((P1.'.*lambda_1)./(z.'-theta) + P2.*lambda_2./(z.'-theta_0) + P3.*lambda_3./(z.'-theta_0_))...
%         ./(lambda_1./(z.'-theta) + lambda_2./(z.'-theta_0) + lambda_3./(z.'-theta_0_));

    P = ((P1.'.*lambda_1)./(z.'-theta) + P2.*lambda_2./(z.'-theta_0) + P3.*lambda_3./(z.'-theta_0_))...
        ./(lambda_1./(z.'-theta) + lambda_2./(z.'-theta_0) + lambda_3./(z.'-theta_0_));

%     n = 80;
%     plot(z,P(n,:),theta(n),P1(n),'x',theta_0,P2,'o',theta_0_,P3,'*');

    I = (w.'.*P)*(1./hatout(z));

    %% special cases - where difference between residues is negligable
    theta2theta_0 = cyclic_distance(theta,theta_0);
    theta2theta_0_ = cyclic_distance(theta,theta_0_);
    theta_02theta_0_ = cyclic_distance(theta_0,theta_0_);

    % collect indices corresponding to special cases
    inds_7i = (theta2theta_0<dtheta) & (theta_02theta_0_>=dtheta);
    inds_7i_ = (theta2theta_0_<dtheta) & (theta_02theta_0_>=dtheta);
    inds_7ii = (theta2theta_0>=dtheta) & (theta_02theta_0_<dtheta);
    inds_7iii = (theta2theta_0<dtheta) & (theta_02theta_0_<dtheta);

    % case one - where \theta=\theta_0\neq\theta_0',
    % need to use derivative (e.g. L'hopital's)
    if sum(inds_7i)>0
%         dhatout = @(theta) p*sin(p*theta);
%         I(inds_7i) = 0;
%         for m=1:length(B_coeffs)
%             I(inds_7i) = I(inds_7i) + B_coeffs*Dhatin_derivs{m}{1}(theta(inds_7i));
%         end
%         I(inds_7i) = (inds_7i)./dhatout(theta(inds_7i));

        % don't need a loop over theta, because they are all approximately
        % \theta_0 :-)
        interp_pts = [theta_0; theta_0_];
        interp_vals = zeros(size(interp_pts));
        dinterp_val_at_theta_0 = 0;
        for m=1:length(B_coeffs)
            interp_vals = interp_vals + B_coeffs(m)*Dhatin{m}(interp_pts);
            dinterp_val_at_theta_0 = dinterp_val_at_theta_0 + B_coeffs(m)*Dhatin_derivs{m}{1}(theta_0);
        end
        interp_coeffs = interp_cubic(interp_pts(1), interp_pts(2), interp_vals(1), interp_vals(2), dinterp_val_at_theta_0);
        interp_poly = polyval(flipud(interp_coeffs),z);
        I(inds_7i) = (w(:,inds_7i).')*(interp_poly./hatout(z));
    end

    % case one' - where \theta=\theta_0'\neq\theta_0,
    %(can't happen in elegant theory but can in inelegant practice, as we are dealing with lots of theta at once)
    % need to use derivative (e.g. L'hopital's)
    if sum(inds_7i_)>0
%         dhatout = @(theta) p*sin(p*theta);
%         I(inds_7i) = 0;
%         for m=1:length(B_coeffs)
%             I(inds_7i) = I(inds_7i) + B_coeffs*Dhatin_derivs{m}{1}(theta(inds_7i));
%         end
%         I(inds_7i) = (inds_7i)./dhatout(theta(inds_7i));

        % don't need a loop over theta, because they are all approximately
        % \theta_0 :-)
%         interp_pts = [theta_0; theta_0_];
        interp_pts = [theta_0_; theta_0];
        interp_vals = zeros(size(interp_pts));
        dinterp_val_at_theta_0_ = 0;
        for m=1:length(B_coeffs)
            interp_vals = interp_vals + B_coeffs(m)*Dhatin{m}(interp_pts);
            dinterp_val_at_theta_0_ = dinterp_val_at_theta_0_ + B_coeffs(m)*Dhatin_derivs{m}{1}(theta_0_);
        end
        interp_coeffs = interp_cubic(interp_pts(1), interp_pts(2), interp_vals(1), interp_vals(2), dinterp_val_at_theta_0_);
        interp_poly = polyval(flipud(interp_coeffs),z);
        I(inds_7i_) = (w(:,inds_7i_).')*(interp_poly./hatout(z));
    end

    % case two - double pole residue at \theta_0=\theta_0'\neq\theta
    % need to use a different interpolation approach, which also
    % interpolates first derivative of integrand at \theta_0
    if sum(inds_7ii)>0
        %
        dinterp_val_at_theta_0 = 0;
        interp_val_at_theta_0 = 0;
        for m=1:length(B_coeffs)
            dinterp_val_at_theta_0 = dinterp_val_at_theta_0 + B_coeffs(m)*Dhatin_derivs{m}{1}(theta_0);
            interp_val_at_theta_0 = interp_val_at_theta_0 + B_coeffs(m)*Dhatin{m}(theta_0);
        end

        for n=inds_7ii
            interp_val_at_theta = 0;
            for m=1:length(B_coeffs)
                interp_val_at_theta = interp_val_at_theta + B_coeffs(m)*Dhatin{m}(theta(n));
            end
            interp_coeffs = interp_cubic(theta_0, theta(n), ...
                                    interp_val_at_theta_0, interp_val_at_theta,...
                                    dinterp_val_at_theta_0);
            interp_poly = polyval(flipud(interp_coeffs),z);
            I(n) = (w(:,n).')*(interp_poly./hatout(z));
        end
    end

    % case three - where \theta=\theta_0=\theta_0',
    % need to use second derivative (e.g. L'hopital's)
    if sum(inds_7iii)>0
        ddhatout = @(theta) -p^2*cos(p*theta);
        I(inds_7iii) = 0;
        for m=1:length(B_coeffs)
            I(inds_7iii) = I(inds_7iii) + B_coeffs(m)*Dhatin_derivs{m}{2}(theta(inds_7iii));
        end
        I(inds_7iii) = I(inds_7iii)./ddhatout(theta(inds_7iii));
    end

end