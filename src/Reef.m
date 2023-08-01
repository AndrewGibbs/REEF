classdef Reef < handle
    %class which stores far-field data,
    %creates new far-field data without solving a new discretised problem,
    %using Cauchy's integral theorem
    
    properties
        FF_in
        Dhat_in
        alpha_in
        kwave
        qppw = 20; %optimised using screen_testRES_test_quad.m (was 50)
        p
        M
        RotShift = 0
        
        % thresholds/tolerances
        %round_error_width_per_can_inc = 0.0015 %optimised using screen_testRES_test_quad.m (was 1e-3/2)
        round_error_width = 0.001 % 0.001 width at which rounding errors become catestrophic
        clas_error_width = 0.15 % 0.1 width at which classical embedding formula becomes inaccurate
        rank_tol = 1e-5
        reconstruction_strategy = 1
        % 1 - least squares via Moore-Penrose inverse
        % 2 - L1 minimisation via basis pursuit
        % 3 - (Approx) maximum (determinant) subvolume, via greedy algorithm
        
        %frequently used variables:
        hat
        %matrix for coefficients
        sampling_matrix

        % inverse/pseudo-inverse
        inv_sampling_matrix
        pinv_tol = 1e-8

        % parameter for basis pursuit (unused)
        rho = 10^9
        
        %analysis parameters:
        alpha_cond
        full_rank % does the collocation matrix have full rank?

        %should the residues be interpolated in the Cauchy integral?
        res_interp = true

        % derivs
        num_derivs = 0
        FF_in_derivs
        Dhat_in_derivs = {}
        deriv_approx_width = 1E-12;
        
    end
    
    methods
        function self = Reef(FF_in_plus_derivs,alpha_in,kwave,p,varargin)
            % note that FF_in and self.FF_in are not the same thing, one is
            % a sub-cell array of the other
            alpha_in = alpha_in(:);
            [length_FF_in, width_FF_in] = size(FF_in_plus_derivs);
            if length_FF_in == 1 || width_FF_in == 1
                self.num_derivs = 0;
                length_FF_in = max(length_FF_in,width_FF_in);
            else
                self.num_derivs = width_FF_in - 1;
            end

            if length_FF_in ~= length(alpha_in)
                error('Number of far-fields (cell columns of first arg) must match number of incident angles');
            end

            % move zeroth derivatives into convinient cell array, if
            % required
            if self.num_derivs>0
                for n=1:length_FF_in
                    self.FF_in{n} = FF_in_plus_derivs{n,1};
%                     for j=1:self.num_derivs
%                         self.FF_in_derivs{n}{j} = FF_in_plus_derivs{n,j+1};
%                     end
                end
            else
                self.FF_in = FF_in_plus_derivs;
                self.FF_in_derivs = {}; % empty cell array, no derivs provided
            end

            for j=1:2
                if self.num_derivs>=j
                    for n=1:length_FF_in
                        self.FF_in_derivs{n}{j} = FF_in_plus_derivs{n,j+1};
                    end
                else % approximate derivative, incase it is needed
                    for n=1:length_FF_in
                        self.FF_in_derivs{n}{j} = ...
                            @(x) (FF_in_plus_derivs{n,j}(x+self.deriv_approx_width/2) - FF_in_plus_derivs{n,j}(x-self.deriv_approx_width/2))...
                                    /self.deriv_approx_width;
                    end
                end
            end

            self.kwave = kwave;
            self.p = p;
            self.M = length(alpha_in);
            
            % set this to be not too small that it interferes with the
            % rounding error limits, and small enough that it accounts for
            % the oscillations in kwave and p
            self.round_error_width = self.round_error_width;
            self.clas_error_width = min(pi/p,max(self.clas_error_width,6*self.round_error_width));
            keep_incs = [];     
            manual_M = false;
            
            for n=1:length(varargin)
                if strcmp(varargin{n},'qppw')
                   self.qppw = varargin{n+1};
               elseif strcmp(varargin{n},'M')
                   self.M = varargin{n+1};
                   manual_M = true;
                   manual_M_val = self.M;
                   if self.M>length(alpha_in)
                      error('There are too few incident angles for this embedding formula');
                   end
               elseif strcmp(varargin{n},'keep')
                   keep_incs = varargin{n+1};
                elseif strcmp(varargin{n},'strategy')
                    self.reconstruction_strategy = varargin{n+1};
                elseif strcmp(varargin{n},'res interp')
                    self.res_interp = varargin{n+1};
                    if self.num_derivs<2
                        warning("Far-field derivative(s) may be approximated");
                    end
                elseif strcmp(varargin{n},'delta')
                    self.pinv_tol = varargin{n+1};
                elseif strcmp(varargin{n},'h')
                    self.round_error_width = varargin{n+1};
                elseif strcmp(varargin{n},'H')
                    self.clas_error_width = varargin{n+1};
                end
            end

            % check for a reasonable reconstruction strategy index
            if ~ismember(self.reconstruction_strategy,1:3)
                error('strategy must be in {1,2,3}');
            end

            if self.reconstruction_strategy <3
                self.M = length(alpha_in);
            end

            self.hat=@(theta,alpha) cos(p*repmat(theta,1,length(alpha)))-(-1)^p*cos(p*(repmat(alpha,length(theta),1)));

            % Part One - The sampling matrix can be singular. 
            % Taking more canonical incident angles gives a better chance
            % of finding a submatrix which is invertible.

            sampling_matrix = get_sampling_matrix(alpha_in,self.FF_in,self.hat);
%             if manual_M
%                 if rank(sampling_matrix,self.rank_tol) < manual_M_val
%                     error("Numerical rank too low, add more inc waves");
% %                 elseif rank(sampling_matrix,self.rank_tol) > manual_M_val
% %                     error("Numerical rank too high, reduce rank_tol");
%                 end
%             end

            switch self.reconstruction_strategy
                case 1 % option 1: use least squares to minimise coefficients 2-norm
                    [self.alpha_in, self.sampling_matrix, self.FF_in, self.inv_sampling_matrix] ...
                        = get_LS_MP(alpha_in, self.FF_in, self.hat, sampling_matrix, self.pinv_tol);
%                 case 3 % option 3: use basis pursuit to minimise coefficients 1-norm
%                     % nothing much to do now
%                     error("This approach doesn't work - due to absence of basis persuit for complex vectors");
%                     self.alpha_in = alpha_in;
%                     self.FF_in = FF_in;
%                     self.sampling_matrix = sampling_matrix;

                case 2 % option 2: maximise determinant of submatrix
                    self.M = manual_M_val;
                    if isempty(keep_incs)
                        [self.alpha_in, self.sampling_matrix, self.FF_in, self.FF_in_derivs, self.full_rank] ...
                            = choose_submatrix(alpha_in, self.FF_in, self.FF_in_derivs, sampling_matrix, self.M, self.rank_tol);
                    else
                        [self.alpha_in, self.sampling_matrix, self.FF_in, self.FF_in_derivs] ...
                            = choose_submatrix_keeps(alpha_in, self.FF_in, self.FF_in_derivs, sampling_matrix, self.M, self.rank_tol,keeps);
                    end
                    self.inv_sampling_matrix = inv(self.sampling_matrix);
            end


            self.alpha_cond = cond(self.sampling_matrix);

            % make some convinient functions for later:
            for m=1:self.M
                self.Dhat_in{m} = @(theta) self.hat(theta,self.alpha_in(m)).*self.FF_in{m}(theta);
                for j=1:2 % first two derivatives is all that is needed
                    self.Dhat_in_derivs{m}{j} = @(theta) get_Dhat_derivs(theta,self.alpha_in(m),j,p,self.FF_in{m},self.FF_in_derivs{m});
                end
            end
        end
        
        function [B,res] = getBmatrix(self,alpha_out)
            alpha_out = alpha_out(:);
            %initialise matrices
            RHS = zeros(self.M,length(alpha_out));
            for m = 1:self.M
                RHS(m,:) = self.hat(alpha_out,self.alpha_in(m).').*self.FF_in{m}(alpha_out);
            end

            u_3p4 = (-1)^(self.p+1)*RHS;

            if ismember(self.reconstruction_strategy,[1 2])
                % multiply by (pseudo) inverse
                B = self.inv_sampling_matrix * u_3p4;
%             else
%                 B = self.sampling_matrix\(u_3p4);
            else % Basis pursuit - strategy 2
                B = zeros(length(alpha_out),length(self.alpha_in));
                % convex optimisation problem needs to be solved one column
                % at a time
                for j = 1:length(alpha_out)
                    B(j,:) = Basis_Pursuit(self.sampling_matrix, u_3p4, self.rho);
                end
            end
%             disp(norm(B)^2);
            res = norm(self.sampling_matrix*B-u_3p4);
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
            H = self.clas_error_width;
            
            [ResDist1,ResDist2,ResVals1,ResVals2,Dhat_res_vals1,Dhat_res_vals2,Dhat_res_vals] = ...
                get_residue_data(theta_out,alpha_out,Z,self.alpha_in,self.FF_in,self.hat);
          

            [res_pairs,~,cases,case4_theta_inds,case5_theta_inds,case6_theta_inds,case7_theta_inds]...
            = choose_cases(theta_out,alpha_out,Z,h,H);
            
             % create tiled copies of b_m for indexing purposes
             for m=1:length(self.alpha_in)
                 B_m_tiled{m} = repmat(B(m,:),length(theta_out),1); 
             end
            
            Laurent_coeff = @(chi) -1./(self.p*sin(self.p*chi));

            % construct the classical embedding formula (case one)
            top = zeros(length(theta_out),length(alpha_out));
            full_hat = (self.hat(theta_out,alpha_out));
            bottom = full_hat;
            for m=1:length(self.alpha_in)
                top = top + B(m,:).*(self.hat(theta_out,self.alpha_in(m)).*self.FF_in{m}(theta_out));
            end
            sum_bm_Dhat_at_theta = top;
            Dout = top./bottom;
            
            % now add a correction if required (cases 2-4)
            
            % case 2
            top = zeros(length(theta_out),length(alpha_out));
            bottom = full_hat.*ResDist1;
            for m=1:length(self.alpha_in)
                top(cases==2) = top(cases==2) + B_m_tiled{m}(cases==2).*(...
                                + Laurent_coeff(ResVals1(cases==2)).*(Dhat_res_vals1{m}(cases==2)).*full_hat(cases==2) ...
                                );
            end
            Dout(cases==2) = Dout(cases==2) + top(cases==2)./bottom(cases==2);
            
            % case 3 (can def be lumped together with case 2 one day)
            top = zeros(length(theta_out),length(alpha_out));
            bottom = (full_hat.*ResDist1.*ResDist2);
            for m=1:length(self.alpha_in)
                top(cases==3) = top(cases==3) + B_m_tiled{m}(cases==3).*(...
                                + Laurent_coeff(ResVals1(cases==3)).*ResDist2(cases==3).*(Dhat_res_vals1{m}(cases==3)).*full_hat(cases==3) ...
                                + Laurent_coeff(ResVals2(cases==3)).*ResDist1(cases==3).*(Dhat_res_vals2{m}(cases==3)).*full_hat(cases==3) ...
                                );
            end
            Dout(cases==3) = Dout(cases==3) + top(cases==3)./bottom(cases==3);
            
       %cases 4-7 all require some kind of Cauchy integral.
       
            for n=1:length(alpha_out)
                for res_pair_index=1:length(res_pairs{n})
                    % case 4 Cauchy integral around residues, not theta_out
                    if ~isempty(case4_theta_inds{n,res_pair_index})
                        theta_0 = res_pairs{n}{res_pair_index}(1);
                        theta_0_ = res_pairs{n}{res_pair_index}(2);
                        a = min([theta_0 theta_0_])-h;
                        b = max([theta_0 theta_0_])+h;
                        [z,w] = Cauchy_box_quad(theta_out(case4_theta_inds{n,res_pair_index}),a,b,h,self.qppw,self.kwave,false);
                        
%                         if false
%                             for m=1:length(self.alpha_in)
%                                 if self.res_interp
%                                     residues = [theta_0; theta_0_];
%                                     p_at_resisudes = self.hat(residues,self.alpha_in(m)).*self.FF_in{m}(residues);
%                                     denom = @(z) self.hat(z,alpha_out(n));
%                                     Dout(case4_theta_inds{n,res_pair_index},n) = Dout(case4_theta_inds{n,res_pair_index},n) + ...
%                                         B(m,n)*Cauchy_interp2(theta_out(case4_theta_inds{n,res_pair_index}), residues, p_at_resisudes, denom, h, self.qppw);
%                                 else
%                                     Dout(case4_theta_inds{n,res_pair_index},n) = Dout(case4_theta_inds{n,res_pair_index},n) + ...
%                                         (B(m,n)*w.'*(self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z)./self.hat(z,alpha_out(n))));
%                                 end
%                             end
%                         else

                            if self.res_interp
                                
                                % slightly streamlined version - sum up
                                % b_m*Dhat(theta_0,alpha_m)
                                % at residue, and interpolate THAT
                                residues = [theta_0; theta_0_];
                                p_at_residues = [0; 0];
                                for m=1:length(self.alpha_in)
                                    p_at_residues = p_at_residues+B(m,n)*self.Dhat_in{m}(residues);
                                end
                                denom = @(z) self.hat(z,alpha_out(n));

                                coalesced_residues = (abs(residues(1)-residues(2))<self.deriv_approx_width);

                                if coalesced_residues % case when \theta_0=\theta_0'
                                    dp_at_residues = 0;
                                    for m=1:length(self.alpha_in)
                                        dp_at_residues = dp_at_residues+B(m,n)*self.Dhat_in_derivs{m}{1}(residues(1));
                                    end
                                    Dout(case4_theta_inds{n,res_pair_index},n) =  Dout(case4_theta_inds{n,res_pair_index},n)+...
                                            Cauchy_interp_case4_coalescing_residues(theta_out(case4_theta_inds{n,res_pair_index}), residues(1), ...
                                            p_at_residues(1), dp_at_residues, denom, h, self.qppw);
                                else % more common case
                                    Dout(case4_theta_inds{n,res_pair_index},n) =  Dout(case4_theta_inds{n,res_pair_index},n)+...
                                            Cauchy_interp_case4(theta_out(case4_theta_inds{n,res_pair_index}), residues, ...
                                            p_at_residues, denom, h, self.qppw);
                                end
        
                            else
                                for m=1:length(self.alpha_in)
                                    Dout(case4_theta_inds{n,res_pair_index},n) = Dout(case4_theta_inds{n,res_pair_index},n) + ...
                                        (B(m,n)*w.'*(self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z)./self.hat(z,alpha_out(n))));
                                end
                            end
%                         end
                    end
                    
                    % Cases 5-7 depend on theta too. Need to cluster the
                    % theta together, for efficient contour integration.

                        %case 5 integral
                        for q=1:2
                            res_index_loc = mod(q,2)+1;
                            res_nearest_index_loc = mod(q+1,2)+1;
                            theta_0 = res_pairs{n}{res_pair_index}(res_index_loc);
                            theta_0_ = res_pairs{n}{res_pair_index}(res_nearest_index_loc);
                            case5U6_theta_inds = [case5_theta_inds{n,res_pair_index,res_index_loc} case6_theta_inds{n,res_pair_index,res_index_loc}];
                            if ~isempty(case5U6_theta_inds)
                                Dout(case5U6_theta_inds,n) = 0;
                                if self.res_interp % takes place outside of m loop
                                    Dout(case5U6_theta_inds,n) = Dout(case5U6_theta_inds,n) +...
                                        Cauchy_interp_case5(theta_out(case5U6_theta_inds), theta_0, B(:,n), ...
                                        self.Dhat_in, self.Dhat_in_derivs, @(z) self.hat(z,alpha_out(n)),...
                                        self.p, h, self.qppw, self.deriv_approx_width);
                                else % evaluate oscillatory contour
                                    a = min([theta_out(case5U6_theta_inds).' theta_0])-h;
                                    b = max([theta_out(case5U6_theta_inds).' theta_0])+h;
                                    [z,w] = Cauchy_box_quad(theta_out(case5U6_theta_inds),a,b,h,self.qppw,self.kwave,true);
                                    for m=1:length(self.alpha_in)
    
    %                                     if self.res_interp
    %                                         Dout(case5U6_theta_inds,n) = Dout(case5U6_theta_inds,n) + ...
    %                                             B(m,n)*Cauchy_interp(theta_out(case5U6_theta_inds), theta_0,...
    %                                                                     @(z) self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z),...
    %                                                                     @(z) self.hat(z,alpha_out(n)),h,...
    %                                                                     self.Dhat_in_derivs{m},...
    %                                                                     self.qppw);
    %                                         
    %                                     else
                                            Dout(case5U6_theta_inds,n) = Dout(case5U6_theta_inds,n) +...
                                                B(m,n)*w.'*(self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z)./self.hat(z,alpha_out(n)));
    %                                     end
                                    end
                                end

                                %case 6 resdiue correction
                                if ~isempty(case6_theta_inds{n,res_pair_index,res_index_loc})
                                    switch res_index_loc
                                        case 1
                                            Dhat_res_vals = Dhat_res_vals1;
                                        case 2
                                            Dhat_res_vals = Dhat_res_vals2;
                                    end
                                    for m=1:length(self.alpha_in)
%                                         Dout(case6_theta_inds{n,res_pair_index,res_index_loc},n) =  Dout(case6_theta_inds{n,res_pair_index,res_index_loc},n)+...
%                                             B(m,n)*(Laurent_coeff(theta_0_).*(Dhat_res_vals{m}(n,res_nearest_index_loc)))./(theta_0_-theta_out(case6_theta_inds{n,res_pair_index,res_index_loc}));
                                        Dout(case6_theta_inds{n,res_pair_index,res_index_loc},n) =  Dout(case6_theta_inds{n,res_pair_index,res_index_loc},n)+...
                                            B(m,n)*(Laurent_coeff(theta_0_).*(Dhat_res_vals{m}(case6_theta_inds{n,res_pair_index,res_index_loc},n)))...
                                            ./(theta_0_-theta_out(case6_theta_inds{n,res_pair_index,res_index_loc}));
                                    end
                                end
                            end
                        end
                    
                    %case 7 - basically the same as case 4, except theta is inside
                    %the integral, so can divide by one
                    if ~isempty(case7_theta_inds{n,res_pair_index})
                        theta_0 = res_pairs{n}{res_pair_index}(1);
                        theta_0_ = res_pairs{n}{res_pair_index}(2);
                        Dout(case7_theta_inds{n,res_pair_index},n) = 0;
                            if self.res_interp
                                Dout(case7_theta_inds{n,res_pair_index},n) =...
                                    Cauchy_interp_case7(theta_out(case7_theta_inds{n,res_pair_index}), theta_0, theta_0_, B(:,n),...
                                        self.Dhat_in, self.Dhat_in_derivs, @(z) self.hat(z,alpha_out(n)), self.p, ...
                                        h, self.qppw, self.deriv_approx_width, sum_bm_Dhat_at_theta(case7_theta_inds{n,res_pair_index}));
                                
                                %TO DO: account for special cases (see notes)

%                                 for m=1:length(self.alpha_in)
%                                     Dout(case7_theta_inds{n,res_pair_index},n) = Dout(case7_theta_inds{n,res_pair_index},n) + ...
%                                         B(m,n)*Cauchy_interp(theta_out(case7_theta_inds{n,res_pair_index}), [theta_0; theta_0_],...
%                                                                         @(z) self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z),...
%                                                                         @(z) self.hat(z,alpha_out(n)),h,self.Dhat_in_derivs{m},...
%                                                                         self.qppw);
%                                 end
                            else
                                a = min([theta_out(case7_theta_inds{n,res_pair_index}).' theta_0 theta_0_])-h;
                                b = max([theta_out(case7_theta_inds{n,res_pair_index}).' theta_0 theta_0_])+h;
                                [z,w] = Cauchy_box_quad(theta_out(case7_theta_inds{n,res_pair_index}),a,b,h,self.qppw,self.kwave,true);
                                for m=1:length(self.alpha_in)
                                    Dout(case7_theta_inds{n,res_pair_index},n) = Dout(case7_theta_inds{n,res_pair_index},n) + ...
                                            B(m,n)*w.'*(self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z)./self.hat(z,alpha_out(n)));
                                end
                            end
%                         end
%                         I = Cauchy_interp_case7(theta_out(case7_theta_inds{n,res_pair_index}), theta_0, theta_0_, B(:,n),...
%                             self.Dhat_in, @(z) self.hat(z,alpha_out(n)), h, self.qppw);
                    end
                end
            end

        end

        function Dout = getFarFieldNaive(self,theta_out,alpha_out)
            %things need to be the right shape
            theta_out = theta_out(:);
            alpha_out = alpha_out(:).';

            % construct the classical embedding formula
            B = self.getBmatrix(alpha_out);
            top = zeros(length(theta_out),length(alpha_out));
            full_hat = (self.hat(theta_out,alpha_out));
            bottom = full_hat;
            for m=1:length(self.alpha_in)
                top = top + B(m,:).*(self.hat(theta_out,self.alpha_in(m)).*self.FF_in{m}(theta_out));
            end
            Dout = top./bottom;

        end
        
        function FFslider(self,theta_test)
            theta_test = theta_test(:);
            fig = uifigure('Position',[100 100 600 600]);
            ax = uiaxes(fig,'Position',[100 175 400 300]);
            fig.Name = 'Far-field pattern';

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

