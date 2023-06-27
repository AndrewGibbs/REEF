function ffDhat_temp = get_sampling_matrix(alpha_in,FF_in,hat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    ffDhat_temp = zeros(length(alpha_in));
    for m = 1:length(alpha_in)
        ffDhat_temp(:,m) = hat(alpha_in,alpha_in(m)).*FF_in{m}(alpha_in);
    end
end