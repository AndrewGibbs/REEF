classdef ReEm < handle
    %class which stores far-field data,
    %creates new far-field data without solving a new discretised problem,
    %using Cauchy's integral theorem
    
    properties
        FF_in
        alpha_in
        kwave
        qppw = 20; %optimised using screen_testRES_test_quad.m (was 50)
        p
        M
        RotShift = 0
        
        % thresholds/tolerances
        round_error_width = 0.001 %optimised using screen_testRES_test_quad.m (was 1e-3/2)
        wavelengths_from_residue = 1.5
        clas_error_width; % width at which classical embedding formula becomes inaccurate
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
            wavelength = 2*pi/kwave; 
            
            % set this to be not too small that it interferes with the
            % rounding error limits, and small enough that it accounts for
            % the oscillations in kwave and p
            self.clas_error_width = min(pi/p,max(self.wavelengths_from_residue*wavelength,6*self.round_error_width));
            
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
            H = self.clas_error_width;
            
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
            
            for m =1:length(self.alpha_in)
                Dhat_res_vals{m} = zeros(size(Z));
                Dhat_res_vals_long = self.hat(Zlong,self.alpha_in(m)).*self.FF_in{m}(Zlong);
                Dhat_res_vals{m} = reshape(Dhat_res_vals_long,size(Z));
                
                for n=1:length(alpha_out)
                    Dhat_res_vals1{m}(:,n) = Dhat_res_vals{m}(n,ResIndices1(n,:));
                    Dhat_res_vals2{m}(:,n) = Dhat_res_vals{m}(n,ResIndices2(n,:));
                end
            end
            %these distances will be wrong over the 2pi-0 jump
            ResDist1 = (ResVals1.'-theta_out);
            ResDist2 = (ResVals2.'-theta_out);
            
            %these two were defined the wrong way around before
            ResVals1 = ResVals1.';
            ResVals2 = ResVals2.';
            
            % now partition input data into all seven cases
            obs_to_res1_dists = abs(ResDist1);
            res1_to_res2_dists = abs(ResDist1-ResDist2);
            
            cases = ones(length(theta_out),length(alpha_out));
%             cases(1.5*h<obs_to_res1_dists & obs_to_res1_dists<=H & H<res1_to_res2_dists) = 2;
%             cases(1.5*h<obs_to_res1_dists & obs_to_res1_dists<=H & 3*h<res1_to_res2_dists & res1_to_res2_dists<=H) = 3;
            %everything below depends on alpha_out, so can't be fully
            %vectorised, so there is little point in the indexing steps for
            %4-7

            one_to_length_theta = 1:length(theta_out);
            for n=1:length(alpha_out)
                [res_pairs{n}, res_pairs_inds{n}] = pair_residues(Z(n,:));
                for j=1:length(res_pairs{n})
                    
                    case4_theta_inds{n,j} = [];
                    case7_theta_inds{n,j} = [];
                    ress = [res_pairs{n}{j}(1) res_pairs{n}{j}(2)];
                    res_dist = abs(ress(1)-ress(2));
                    
                    
                    case5_theta_inds{n,j,1} = [];
                    case6_theta_inds{n,j,1} = [];
                    case5_theta_inds{n,j,2} = [];
                    case6_theta_inds{n,j,2} = [];
                    
                    if  res_dist <= 3*h
                        theta_pair_dist = min(abs(theta_out-ress(1)),abs(theta_out-ress(2)));
                        case4_theta_inds{n,j} = one_to_length_theta(1.5*h < theta_pair_dist & theta_pair_dist <= H);
                        case7_theta_inds{n,j} = one_to_length_theta(theta_pair_dist <= 1.5*h);
                    elseif 3*h < res_dist && res_dist <= H
                        for q = 1:2 % loop over both residues in pair
                            theta_dist = abs(theta_out-ress(q));
                            % case 3:
                            cases(h<theta_dist & theta_dist<=H,n) = 3;
                            %case 6:
                            case6_theta_inds{n,j,q} = one_to_length_theta(theta_dist <= h);
                        end
                    elseif H < res_dist
                        for q = 1:2 % loop over both residues in pair
                            theta_dist = abs(theta_out-ress(q));
                            % case 2:
                            cases(h<theta_dist & theta_dist<=H,n) = 2;
                            %case 5:
                            case5_theta_inds{n,j,q} = one_to_length_theta(theta_dist <= h);
                        end
                    end
                end
            end
            
             % create tiled copied of b_m for indexing purposes
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
            
       %cases 4-7 all require some kind of Cauchy integral
       
            for n=1:length(alpha_out)
                for res_pair_index=1:length(res_pairs{n})
                    if ~isempty(case4_theta_inds{n,res_pair_index})
                        theta_0 = res_pairs{n}{res_pair_index}(1);
                        theta_0_ = res_pairs{n}{res_pair_index}(2);
                        a = min([theta_0 theta_0_])-h;
                        b = max([theta_0 theta_0_])+h;
                        [z,w] = Cauchy_box_quad(theta_out(case4_theta_inds{n,res_pair_index}),a,b,h,self.qppw,self.kwave,false);
                        for m=1:length(self.alpha_in)
                            Dout(case4_theta_inds{n,res_pair_index},n) = Dout(case4_theta_inds{n,res_pair_index},n) + ...
                                    (B(m,n)*w.'*(self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z)./self.hat(z,alpha_out(n))));
                        end
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
                                a = min([theta_out(case5U6_theta_inds).' theta_0])-h;
                                b = max([theta_out(case5U6_theta_inds).' theta_0])+h;
                                [z,w] = Cauchy_box_quad(theta_out(case5U6_theta_inds),a,b,h,self.qppw,self.kwave,true);
                                Dout(case5U6_theta_inds,n) = 0;
                                for m=1:length(self.alpha_in)
                                    Dout(case5U6_theta_inds,n) = Dout(case5U6_theta_inds,n) +...
                                        B(m,n)*w.'*(self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z)./self.hat(z,alpha_out(n)));

                                    %case 6 resdiue correction
                                    if ~isempty(case6_theta_inds{n,res_pair_index})
                                        Dout(case6_theta_inds{n,res_pair_index,res_index_loc},n) =  Dout(case6_theta_inds{n,res_pair_index,res_index_loc},n)+...
                                            B(m,n)*(Laurent_coeff(theta_0_).*(Dhat_res_vals{m}(n,res_nearest_index_loc)))./(theta_0_-theta_out(case6_theta_inds{n,res_pair_index,res_index_loc}));
                                    end
                                end
                            end
                        end
%                     end
                    
                    %case 7 - basically the same as case 4, except theta is inside
                    %the integral, so can divide by one
                    if ~isempty(case7_theta_inds{n,res_pair_index})
                        theta_0 = res_pairs{n}{res_pair_index}(1);
                        theta_0_ = res_pairs{n}{res_pair_index}(2);
                        a = min([theta_out(case7_theta_inds{n,res_pair_index}).' theta_0 theta_0_])-h;
                        b = max([theta_out(case7_theta_inds{n,res_pair_index}).' theta_0 theta_0_])+h;
                        [z,w] = Cauchy_box_quad(theta_out(case7_theta_inds{n,res_pair_index}),a,b,h,self.qppw,self.kwave,true);
                        Dout(case7_theta_inds{n,res_pair_index},n) = 0;
                        for m=1:length(self.alpha_in)
                            Dout(case7_theta_inds{n,res_pair_index},n) = Dout(case7_theta_inds{n,res_pair_index},n) + ...
                                    B(m,n)*w.'*(self.hat(z,self.alpha_in(m)).*self.FF_in{m}(z)./self.hat(z,alpha_out(n)));
                        end
                    end
                end
            end

            % make the limits for the endpoints of the regions which we integrate
            % around to avoid rounding errors

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

