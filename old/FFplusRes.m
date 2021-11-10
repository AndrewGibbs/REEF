function [Dout, Rout] = FFplusRes(theta,D,B,p,Z,alpha_in,alpha_out)
% proof of concept for the residue-based approach to contour integrals

    hat=@(theta,alpha) cos(p*repmat(theta,1,length(alpha)))-(-1)^p*cos(p*(repmat(alpha,length(theta),1)));
    h = @(z,theta) p*sin(p*z)*(z-theta);
    %-hat(z,alpha_out) - p*sin(p*z)*(z-theta);
    
    H = 1;
    for z = Z
        H = H.*h(z,theta);
    end
    Dout = 0;
    Rout = 0;
%     for m=1:length(alpha_in)
%         R = 0;
%         for z=Z
%             w_z = 1;
%             for z_=Z
%                 if z ~= z_
%                     w_z = w_z.*h(z_,theta);
%                 end
%             end
%             R = R + w_z.*D{m}(z).*hat(z,alpha_in(m)).*hat(theta,alpha_out);
%         end
%         Dout = Dout + B(m).*(R + H.*D{m}(theta).*hat(theta,alpha_in(m)))./(H.*hat(theta,alpha_out));
%         Rout = Rout + B(m).*(R)./(H.*hat(theta,alpha_out));
%     end
    for m=1:length(alpha_in)
        R = 0;
        for z=Z
            w_z = 1;
            for z_=Z
                if z ~= z_
                    w_z = w_z.*h(z_,theta);
                end
            end
            R = R + w_z.*D{m}(z).*hat(z,alpha_in(m)).*hat(theta,alpha_out);
        end
        Dout = Dout + B(m).*(R + H.*D{m}(theta).*hat(theta,alpha_in(m)));
        Rout = Rout + B(m).*(R);
    end
    Dout = Dout./(H.*hat(theta,alpha_out));
    Rout = Rout./(H.*hat(theta,alpha_out));

end