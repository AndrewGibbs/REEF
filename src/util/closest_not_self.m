function [Z_val,Z_index] = closest_not_self(z,Z)
    R = abs(z-Z);
    if sum(R==0) >= 2
        Z_val = z;
        Z_index = find(Z==z,1,'first');
    else
        R(R==0)=2*pi;
        [~,Z_index] = min(R);
        Z_val = Z(Z_index);
    end
end

