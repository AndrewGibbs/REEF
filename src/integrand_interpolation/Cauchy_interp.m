function I = Cauchy_interp(theta_set, theta_0, Dhatin, hatout, h, deriv, Nquad)
    small_shift = 10^5*eps;

%     if isempty(deriv)
%         % approximate derivative if none provided
%         dDhatin = @(x) (Dhatin(x+small_shift/2)-Dhatin(x-small_shift/2))/(small_shift);
%     else
%         dDhatin = deriv{1};
%     end
% 
%     if length(deriv)<=1
%         % approximate second derivative if none provided
%         ddDhatin = @(x) (dDhatin(x+small_shift/2)-dDhatin(x-small_shift/2))/(small_shift);
%     else
%         ddDhatin = deriv{2};
%     end
    dDhatin = deriv{1};
    ddDhatin = deriv{2};

    theta_set = theta_set(:);

    % get focii of ellipse (sort of)

    a = min([theta_set; theta_0])-h;
    b = max([theta_set; theta_0])+h;

    if (length(theta_0)==2) && (theta_0(1)==theta_0(2))
        theta_0(2) = theta_0(2) + small_shift;
    end

    count = 0;
    I = zeros(size(theta_set));
%     I_check = zeros(size(theta_set));

    [z,w] = Cauchy_box_quad(theta_set, a, b, h, Nquad, 1, false);
    % can be made more efficient and looped over
        
    for theta = theta_set'
            count = count + 1;
        if abs(theta-theta_0(1))<small_shift
            interp_pts = theta_0;
            interp_vals = Dhatin(interp_pts);
            interp_coeffs = interp_cubic(interp_pts(1),interp_pts(2),interp_vals(1),interp_vals(2),dDhatin(theta_0(1)));
            interp_poly = polyval(flipud(interp_coeffs),z);
            I(count) = (w(:,count).')*(interp_poly./hatout(z));

        else
    %         if abs(theta-theta_0(1))<small_shift
    %             theta = theta - small_shift;
    %         elseif length(theta_0)==2 && theta==theta_0(2)
    %             theta = theta - small_shift;
    %         end
            interp_pts = [theta; theta_0];
            interp_vals = Dhatin(interp_pts);
            interp_poly = barylag([interp_pts interp_vals],z);
            I(count) = (w(:,count).')*(interp_poly./hatout(z));
    %         I_check(count) = (w(:,count).')*(Dhatin(z)./hatout(z));

        end
    end
    
end