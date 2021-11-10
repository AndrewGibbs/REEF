classdef Embed
    %class which stores far-field data,
    %creates new far-field data without solving a new discretised problem,
    %using Cauchy's integral theorem
    
    properties
        FF_in
        alpha_in
        kwave
        qppw = 25; %optimised for pmax = 4, 1e-3 accurate far-fields
        p
        M
        RotShift = 0
        cap_ends;
        %frequently used variables:
        hat
        unstable_period
        %integration nodes
        z_rail
        z_Lcap
        z_Rcap
        %weights
        w_rail
        w_Lcap
        w_Rcap
        %values
        D_rail
        D_Lcap
        D_Rcap
        %matrix for coefficients
        ffDhat
        %option for dividing by one, extra computation but should improve
        %accuracy of contour integrals:
        divide_by_one = true
        
        %analysis parameters:
        total_quad_pts
        imag_dist
        alpha_cond
        full_rank % does the collocation matrix have full rank?
        
    end
    
    methods
        function self = Embed(FF_in,alpha_in,kwave,p,varargin)
            self.FF_in = FF_in;
            alpha_in = alpha_in(:);
            if length(FF_in) ~= length(alpha_in)
                error('Number of far-fields (first arg) must match number of incident angles');
            end
            self.kwave = kwave;
            self.p = p;
            self.M = length(alpha_in);
            self.unstable_period = 2*pi/p;
            self.cap_ends = [self.unstable_period self.unstable_period/2];
            
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
            
            %construct the 'rail' component of the complex contour:
            [ x, w ] = gauss_quad_wave_split3(0, 2*pi, self.qppw, kwave,  2*pi );
            self.z_rail = [x+1i*self.imag_dist; x-1i*self.imag_dist];
            self.w_rail = [-w; w]/(2i*pi);
            self.D_rail = get_D_matrix(self.M,FF_in,self.z_rail,self.w_rail,self.hat,alpha_in);
            clear x w;
            
            %construct the 'left caps' of the complex contour:
            Ldist = -self.cap_ends;
            for j = [1 2]
                [ x, w_x ] = gauss_quad_wave_split3(Ldist(j), 0, self.qppw, kwave,  abs(Ldist(j)) );
                [ y, w_y ] = gauss_quad_wave_split3(Ldist(j)-1i*self.imag_dist, Ldist(j)+1i*self.imag_dist, self.qppw, kwave,  2*self.imag_dist );
                self.z_Lcap{j} = [x+1i*self.imag_dist; y; x-1i*self.imag_dist];
                self.w_Lcap{j} = [-w_x; -w_y; w_x]/(2i*pi);
                self.D_Lcap{j} = get_D_matrix(self.M,FF_in,self.z_Lcap{j},self.w_Lcap{j},self.hat,alpha_in);
                clear x y w_x w_y;
            end
            
            %construct the 'right caps' of the complex contour:
            Rdist = self.cap_ends;
            for j = [1 2]
                [ x, w_x ] = gauss_quad_wave_split3(2*pi, 2*pi+Rdist(j), self.qppw, kwave,  abs(Rdist(j)) );
                [ y, w_y ] = gauss_quad_wave_split3(2*pi+Rdist(j)-1i*self.imag_dist, 2*pi+Rdist(j)+1i*self.imag_dist, self.qppw, kwave,  2*self.imag_dist );
                self.z_Rcap{j} = [x+1i*self.imag_dist; y; x-1i*self.imag_dist];
                self.w_Rcap{j} = [-w_x; w_y; w_x]/(2i*pi);
                self.D_Rcap{j} = get_D_matrix(self.M,FF_in,self.z_Rcap{j},self.w_Rcap{j},self.hat,alpha_in);
                clear x y w_x w_y;
            end
            
            
            ffDhat_temp = zeros(length(alpha_in));
            for m = 1:length(alpha_in)
                ffDhat_temp(:,m) = self.hat(alpha_in,alpha_in(m)).*FF_in{m}(alpha_in);
            end
            
            % The collocation matrix can be singular. We haven't been able
            % to develop theory on how to avoid this entirely, however
            % taking more canonical incident angles gives a better chance
            % of finding a submatrix which is invertible.
            
            % the next several lines are based on mathworks file exchange code at:
            %https://uk.mathworks.com/matlabcentral/fileexchange/77437-extract-linearly-independent-subset-of-matrix-columns
            
           tol = 1e-5;
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
            
            %store the number of quad points for analysis
            self.total_quad_pts = length([self.w_rail; self.w_Lcap{1}; self.w_Rcap{1}]);
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
        
        function D = getFarField(self,theta_out,alpha_out,Cauchy_frac_input)
            %optional extra input for fast plotting with slider:
            if nargin == 3
                Cauchy_frac_in = false;
            else
                Cauchy_frac_in = true;
            end
            
            alpha_out = alpha_out(:).';
            D = zeros(length(theta_out),length(alpha_out));
            sub_alpha_inds = -(-1)^self.p*cos(self.p*alpha_out)>=0;
            sub_alpha_inds_split = {sub_alpha_inds,~sub_alpha_inds};
            for j=[1 2]
                alpha_out_split = alpha_out(sub_alpha_inds_split{j});
                if ~isempty(alpha_out_split)
                    B = self.getBmatrix(alpha_out_split);
                    z = [self.z_rail;
                        self.z_Lcap{j};
                        self.z_Rcap{j}];
                    Dhat_in = [self.D_rail; self.D_Lcap{j}; self.D_Rcap{j};];
                    if Cauchy_frac_in
                        Cauchy_frac = Cauchy_frac_input{j};
                    else
                        Cauchy_frac = 1./(z.'-theta_out);
                    end
                    hat_frac = 1./(self.hat(z,alpha_out_split));
                    if self.divide_by_one
                        div_by_one = 1./(Cauchy_frac*[self.w_rail; self.w_Lcap{j}; self.w_Rcap{j}]);
                    else
                        div_by_one = 1;
                    end
                    D(:,sub_alpha_inds_split{j}) = Cauchy_frac*((Dhat_in*B).*hat_frac) .* div_by_one;
                end
            end
            
        end
        
        function FFslider(self,theta_test)
            %precompute this stuff to save flops each time the slider moves
            for j=[1 2]
                z = [self.z_rail;
                    self.z_Lcap{j};
                    self.z_Rcap{j}];
                Cauchy_frac{j} = 1./(z.'-theta_test);
            end
            fig = uifigure('Position',[100 100 600 600]);
            ax = uiaxes(fig,'Position',[100 175 400 300]);
            %ax.legend('real','imaginary');
            fig.Name = 'Far-field pattern';
            %title('inc field slider');

            sld = uislider(fig,...
                'Position',[150 150 300 3],...
                'ValueChangedFcn',@(sld,event) updateGauge(event,ax));
            sld.MajorTickLabels = {'0','pi/2','pi','3pi/2','2pi'};
            sld.MajorTicks = [0 pi/2 pi 3*pi/2 2*pi];
            sld.Limits = [0 2*pi];

            function updateGauge(event,ax)
                Evals = self.getFarField(theta_test,event.Value,Cauchy_frac);
                %plot(ax,theta_test,real(Evals),theta_test,imag(Evals));
                plot(ax,theta_test,real(Evals),theta_test,imag(Evals));
                ax.YLim = ([-2.5*self.kwave 2.5*self.kwave]);
                ax.XLim = ([0 2*pi]);
                ax.XTick = ([0 pi/2 pi 3*pi/2 2*pi]);
                ax.XTickLabel = ({'0','pi/2','pi','3pi/2','2pi'});
            end

        end
    end
end

