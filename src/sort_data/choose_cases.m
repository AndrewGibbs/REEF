function [res_pairs,res_pairs_inds,cases,case4_theta_inds,case5_theta_inds,case6_theta_inds,case7_theta_inds]...
            = choose_cases(theta_out,alpha_out,Z,h,H)
        
        cases = ones(length(theta_out),length(alpha_out));
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

                if  res_dist < 3*h
                    theta_pair_dist = min(abs(theta_out-ress(1)),abs(theta_out-ress(2)));
                    case4_theta_inds{n,j} = one_to_length_theta(1.5*h <= theta_pair_dist & theta_pair_dist < H);
                    case7_theta_inds{n,j} = one_to_length_theta(theta_pair_dist < 1.5*h);
                elseif 3*h <= res_dist && res_dist < H
                    for q = 1:2 % loop over both residues in pair
                        theta_dist = abs(theta_out-ress(q));
                        % case 3:
                        cases(h<=theta_dist & theta_dist<H,n) = 3;
                        %case 6:
                        case6_theta_inds{n,j,q} = one_to_length_theta(theta_dist < h);
                    end
                elseif H <= res_dist
                    for q = 1:2 % loop over both residues in pair
                        theta_dist = abs(theta_out-ress(q));
                        % case 2:
                        cases(h<=theta_dist & theta_dist<H,n) = 2;
                        %case 5:
                        case5_theta_inds{n,j,q} = one_to_length_theta(theta_dist < h);
                    end
                end
            end
        end
end

