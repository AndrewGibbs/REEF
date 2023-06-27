function I = Cauchy_interp(theta_set, theta_0, Dhatin, hatout, h)
if nargin == 6
    %use the divide by one trick
    div1 = true;
end
    qppw = 20;
    small_shift = 10^5*eps;
    theta_set = theta_set(:);

    % get focii of ellipse (sort of)

    a = min([theta_set; theta_0])-h;
    b = max([theta_set; theta_0])+h;

    if (length(theta_0)==2) && (theta_0(1)==theta_0(2))
        theta_0(2) = theta_0(2) + small_shift;
    end

    count = 0;
    I = zeros(size(theta_set));
    I_check = zeros(size(theta_set));

    [z,w] = Cauchy_box_quad(theta_set, a, b, h, qppw, 1, false);
    % can be made more efficient and looped over
        
    for theta = theta_set'
            count = count + 1;
        if abs(theta-theta_0(1))<small_shift
            interp_pts = theta_0;
            interp_vals = Dhatin(interp_pts);
            deriv_approx = (Dhatin(theta_0(1)+small_shift/2)-Dhatin(theta_0(1)-small_shift/2))/(small_shift);
            interp_coeffs = interp_cubic(interp_pts(1),interp_pts(2),interp_vals(1),interp_vals(2),deriv_approx);
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