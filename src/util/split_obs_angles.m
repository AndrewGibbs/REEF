function [closest_res_indices, next_closest_res_indices, closest_res, nearby_nearby_res] ...
            = split_obs_angles(theta, T)
    theta = theta(:).';
    T = T(:).';
    % get residue which is closest to observation point
    R = abs(theta - T.');
    [~, closest_res_indices] = min(R);
    closest_res = T(closest_res_indices);
    
    % now get residue which is closest to that residue
    R2 = abs(T - T.') + 2*pi*eye(length(T));
    %R2 = exp(1i*T(closest_res_indices)) - exp(1i*T.');
    [~, closest_res_to_res] = min(R2);
    next_closest_res_indices = closest_res_to_res(closest_res_indices);
    nearby_nearby_res = T(next_closest_res_indices);
    
%     function r = circle_dist(x,y)
%        r = min(min(abs(x-y), abs(x+2*pi-y)), abs(x-2*pi-y));
%     end
end