function [ResDist1,ResDist2,ResVals1,ResVals2,Dhat_res_vals1,Dhat_res_vals2,Dhat_res_vals] = get_residue_data(theta_out,alpha_out,Z,alpha_in,FF_in,hat)

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
            
            for m =1:length(alpha_in)
                Dhat_res_vals{m} = zeros(size(Z));
                Dhat_res_vals_long = hat(Zlong,alpha_in(m)).*FF_in{m}(Zlong);
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
end

