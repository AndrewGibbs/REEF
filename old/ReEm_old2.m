classdef ReEm < handle
    %class which stores far-field data,
    %creates new far-field data without solving a new discretised problem,
    %using Cauchy's integral theorem
    
    properties
        FF_in
        alpha_in
        kwave
        qppw = 50; %optimised for pmax = 4, 1e-3 accurate far-fields
        p
        M
        RotShift = 0
        
        % thresholds/tolerances
        round_error_width = 1e-3/2
        rank_tol = 1e-5
        
        
        %frequently used variables:
        hat
        %matrix for coefficients
        ffDhat
        
        %analysis parameters:
        total_quad_pts
        imag_dist
        alpha_cond
        full_rank % does the collocation matrix have full rank?
        
    end
    
    methods
        function self = ReEm(FF_in,alpha_in,kwave,p,varargin)
            alpha_in = alpha_in(:);
            if length(FF_in) ~= length(alpha_in)
                error('Number of far-fields (first arg) must match number of incident angles');
            end
            self.kwave = kwave;
            self.p = p;
            self.M = length(alpha_in);
            
            self.imag_dist = 1/kwave;
            
            for n=1:length(varargin)
               if strcmp(varargin{n},'imag dist')
                   self.imag_dist = varargin{n+1};
               elseif strcmp(varargin{n},'qppw')
                   self.qppw = varargin{n+1};
               elseif strcmp(varargin{n},'M')
                   self.M = varargin{n+1};
                   if self.M>length(alpha_in)
                      error('There are too few incident angles for this embedding formula');
                   end
               end
            end
            
            self.hat=@(theta,alpha) cos(p*repmat(theta,1,length(alpha)))-(-1)^p*cos(p*(repmat(alpha,length(theta),1)));
            % .  hat=@(theta,alpha) cos(p*repmat(theta,1,length(alpha)))-(-1)^p*cos(p*(repmat(alpha,length(theta),1)));
            
            % Part One - The collocation matrix can be singular. 
            % sTaking more canonical incident angles gives a better chance
            % of finding a submatrix which is invertible.
            
            ffDhat_temp = zeros(length(alpha_in));
            for m = 1:length(alpha_in)
                ffDhat_temp(:,m) = self.hat(alpha_in,alpha_in(m)).*FF_in{m}(alpha_in);
            end
            
            % the next several lines are based on mathworks file exchange code at:
            %https://uk.mathworks.com/matlabcentral/fileexchange/77437-extract-linearly-independent-subset-of-matrix-columns
           tol = self.rank_tol;
           if max(max(abs(ffDhat_temp)))<tol %check if all entries are almost zero
               r = 0;
               E = 1:self.M;
           else
               [~, R, E] = qr(ffDhat_temp,0); 
               if ~isvector(R)
                diagr = abs(diag(R));
               else
                diagr = R(1);   
               end
               %Rank estimation
               r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation, should be <=M
               uniqe_m=sort(E(1:self.M));
           end
            
            if r < self.M
               warning(sprintf('%dx%d collocation matrix has rank %d - results will be highly inaccurate - please consider different/more indicent angles',self.M,self.M,r));
               self.full_rank = false;
               self.alpha_in = alpha_in;
               self.ffDhat = ffDhat_temp(1:self.M,1:self.M);
               uniqe_m_trunc = 1:self.M;
            else
               uniqe_m_trunc = uniqe_m(1:self.M);
               self.alpha_in = alpha_in(uniqe_m_trunc);
               self.full_rank = true;
               self.ffDhat = ffDhat_temp(uniqe_m_trunc,uniqe_m_trunc);
            end
            self.alpha_cond = cond(self.ffDhat);
            for m=1:self.M
                self.FF_in{m} = FF_in{uniqe_m_trunc(m)};
            end
        end
        
        function B = getBmatrix(self,alpha_out)
            alpha_out = alpha_out(:);
            %initialise matrices
            RHS = zeros(self.M,length(alpha_out));
            for m = 1:self.M
                RHS(m,:) = self.hat(alpha_out,self.alpha_in(m).').*self.FF_in{m}(alpha_out);
            end
            u_3p4 = (-1)^(self.p+1)*RHS;
            B = self.ffDhat\(u_3p4);
        end
        
        function Dout = getFarField(self,theta_out,alpha_out)
            %things need to be the right shape
            theta_out = theta_out(:);
            alpha_out = alpha_out(:).';
            
            %get the residues
            [Z_,~] = getBadAngles(alpha_out,self.p);
            Z_ = sort(Z_')';
            %extend the Z vector a little to catch residues just outside of
            %[0,2pi]
            Z = [Z_(:,(end-1):end)-2*pi Z_ Z_(:,1:2)+2*pi];
            %get the coefficients b_m(\alpha)
            B = self.getBmatrix(alpha_out);
            h = self.round_error_width;
            
            ResIndices1 = zeros(length(alpha_out), length(theta_out));
            ResIndices2 = zeros(length(alpha_out), length(theta_out));

            for n = 1:length(alpha_out)
                [res_index1, res_index2] = split_obs_angles(theta_out, Z(n,:));
                ResIndices1(n,:) = res_index1;
                ResIndices2(n,:) = res_index2;
            end
            
            for n=1:length(alpha_out)
                ResVals1(n,:) = Z(n,ResIndices1(n,:));
                ResVals2(n,:) = Z(n,ResIndices2(n,:));
            end
            
            % store the residue values so they are not recomputed loads of times
            
            Zlong = reshape(Z,numel(Z),1);
            Res1IndsLong = reshape(ResIndices1,numel(ResIndices1),1);
            Res2IndsLong = reshape(ResIndices2,numel(ResIndices2),1);
            ResVals1long = reshape(ResVals1,numel(ResVals1),1);
            ResVals2long = reshape(ResVals2,numel(ResVals2),1);
            
            for m =1:length(self.alpha_in)
                Dhat_res_vals{m} = zeros(size(Z));
                Dhat_res_vals_long = self.hat(Zlong,self.alpha_in(m)).*self.FF_in{m}(Zlong);
                Dhat_res_vals{m} = reshape(Dhat_res_vals_long,size(Z));
                
                %Dhat_res_vals_short = self.hat(Z(n,:), self.alpha_in(m)).*self.FF_in{m}(Z(n,:));
                for n=1:length(alpha_out)
%                     Dhat_res_vals1_long(n,:) = Dhat_res_vals{m}(n,ResIndices1(n,:)).';
                    Dhat_res_vals1{m}(n,:) = Dhat_res_vals{m}(n,ResIndices1(n,:)).';
                    Dhat_res_vals2{m}(n,:) = Dhat_res_vals{m}(n,ResIndices2(n,:)).';
                end
                %Dhat_res_vals1_long = Dhat_res_vals_short(ResIndices1(n,:)).';
%                 Dhat_res_vals1_long_old = self.hat(ResVals1long, self.alpha_in(m)).*self.FF_in{m}(ResVals1long);
%                 Dhat_res_vals1{m} = reshape(Dhat_res_vals1_long_old, size(ResVals1)); %5x1000
                
                
%                 for n=1:length(alpha_out)
%                     %Dhat_res_vals2_long(n,:) = Dhat_res_vals{m}(n,ResIndices2(n,:)).';
%                 end
%                 Dhat_res_vals2_short = self.hat(Z(n,:), self.alpha_in(m)).*self.FF_in{m}(Z(n,:));
                %Dhat_res_vals2_long = Dhat_res_vals_short(ResIndices2(n,:)).';
%                 Dhat_res_vals2_long_old = self.hat(ResVals2long, self.alpha_in(m)).*self.FF_in{m}(ResVals2long);
%                 Dhat_res_vals2{m} = reshape(Dhat_res_vals2_long_old, size(ResVals2));
            end
            %these distances will be wrong over the 2pi-0 jump
            ResDist1 = (ResVals1.'-theta_out);
            ResDist2 = (ResVals2.'-theta_out);

            Laurent_coeff = @(chi) -1./(self.p*sin(self.p*chi));

            top = zeros(length(alpha_out), length(theta_out));
            bottom = (self.hat(theta_out,alpha_out).*ResDist1.*ResDist2).';
            for m=1:length(self.alpha_in)
                top = top + (B(m,:).*(...
                                ResDist1.*ResDist2.*self.hat(theta_out,self.alpha_in(m)).*self.FF_in{m}(theta_out) ...
                                + Laurent_coeff(ResVals1.').*ResDist2.*(Dhat_res_vals1{m}).'.*self.hat(theta_out,alpha_out) ...
                                + Laurent_coeff(ResVals2.').*ResDist1.*(Dhat_res_vals2{m}).'.*self.hat(theta_out,alpha_out)...
                            )).';
            end
%                                 + a(ResVals1.').*ResDist2.*(Dhat_res_vals{m}(ResIndices1)).'.*self.hat(theta_out,alpha_out) ...
%                                 + a(ResVals2.').*ResDist1.*(Dhat_res_vals{m}(ResIndices2)).'.*self.hat(theta_out,alpha_out)...
%  
            Dout = top./bottom;

            % make the limits for the endpoints of the regions which we integrate
            % around to avoid rounding errors
            
            type = zeros(length(alpha_out));
            for n=1:length(alpha_out)
                if abs(Z(n,2)-Z(n,1))<2*h %case 2a
                    roundL{n} = Z(n,1:2:end);
                    roundR{n} = Z(n,2:2:end);
                    type(n) = 2.1;
                elseif abs(Z(n,3)-Z(n,2))<2*h %case 2b
                    roundL{n} = Z(n,2:2:end);
                    roundL{n} = [0 roundL{n} Z(n,end)];
                    roundR{n} = Z(n,3:2:end);
                    roundR{n} = [Z(n,1) roundR{n} 2*pi];
                    type(n) = 2.2;
                else %case 1 (most common (but still uncommon))
%                     roundL{n} = Z;
%                     roundR{n} = Z;
                    roundL{n} = Z(n,:);
                    roundR{n} = Z(n,:);
                    type(n) = 1;
                end
            end

            % now compute the values in the danger zone using a Cauchy
            % integral, which is far more stable and can avoid the rounding
            % errors:
            
            for n=1:length(alpha_out)
                %first replace the values where the residues are too close
                %together. Case two in my notes.
                if round(type(n)) == 2
                    num_clumps = numel(Z)/length(alpha_out)/2;
                    for q = 1:num_clumps
                        a = roundL{n}(q);
                        b = roundR{n}(q);
                        theta_dodgy_inds = (a-h<theta_out & theta_out<b+h);
                        if sum(theta_dodgy_inds)>0
                            [z,w] = Cauchy_box_quad(theta_out(theta_dodgy_inds), a-2*h, b+2*h, h, self.qppw, self.kwave, false);
                            theta_dodgy = theta_out(theta_dodgy_inds);
                            replacement = zeros(length(theta_dodgy),1);
                            for m=1:length(self.alpha_in)
                                replacement = replacement ...
                                    + B(m,n).*(w.'*(self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z)./self.hat(z,alpha_out(n)))).'...
                                    + B(m,n).*(self.hat(theta_dodgy,self.alpha_in(m)).*self.FF_in{m}(theta_dodgy)./(self.hat(theta_dodgy,alpha_out(n))));
                            end
                            Dout(n,theta_dodgy_inds) = replacement;
                        end
                    end
          
                end
                
                for q =1:length(roundL{n})
                    theta_dodgy_inds = (roundL{n}(q)-h<theta_out & theta_out<roundR{n}(q)+h);
                    if sum(theta_dodgy_inds)>0
                        theta_dodgy = theta_out(theta_dodgy_inds);
                        [z,w] = Cauchy_box_quad(theta_dodgy,roundL{n}(q)-2*h,roundR{n}(q)+2*h,h,self.qppw,self.kwave);
                        replacement = zeros(length(theta_dodgy),1);
                        for m=1:length(self.alpha_in)
                            replacement = replacement + B(m,n)*w.'*(self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z)./self.hat(z,alpha_out(n)));
                            if type(n) == 1
                                replacement = replacement + B(m,n)*(...
                                    Laurent_coeff(ResVals2(n,theta_dodgy_inds).').*(Dhat_res_vals{m}(ResIndices2(n,theta_dodgy_inds))).'...
                                    );
                            end
                        end
                        Dout(n,theta_dodgy_inds) = replacement;
                    end
                end
            end
        end
        
        function FFslider(self,theta_test)
            theta_test = theta_test(:);
            fig = uifigure('Position',[100 100 600 600]);
            ax = uiaxes(fig,'Position',[100 175 400 300]);
%             ax.legend('real','imaginary');
            fig.Name = 'Far-field pattern';
%             title('inc field slider');

            sld = uislider(fig,...
                'Position',[150 150 300 3],...
                'ValueChangedFcn',@(sld,event) updateGauge(event,ax));
            sld.MajorTickLabels = {'0','pi/2','pi','3pi/2','2pi'};
            sld.MajorTicks = [0 pi/2 pi 3*pi/2 2*pi];
            sld.Limits = [0 2*pi];

            function updateGauge(event,ax)
                Evals = self.getFarField(theta_test,event.Value);
                plot(ax,theta_test,real(Evals),theta_test,imag(Evals));
                ax.YLim = ([-2.5*self.kwave 2.5*self.kwave]);
                ax.XLim = ([0 2*pi]);
                ax.XTick = ([0 pi/2 pi 3*pi/2 2*pi]);
                ax.XTickLabel = ({'0','pi/2','pi','3pi/2','2pi'});
            end

        end
    end
end

