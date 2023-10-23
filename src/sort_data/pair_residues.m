function [Z_,Z_inds] = pair_residues(Zn)
%clumps residues into pairs
    count = 1;
    num_pairs = (length(Zn)/2);
    if abs(Zn(1)-Zn(2)) < abs(Zn(3)-Zn(2))
        for j=0:(num_pairs-1)
            z_1 = Zn(2*j+1);
            z_2 = Zn(2*j+2);
            z_inds = [2*j+1 2*j+2];
            store_pair();
        end
    else
        for j=0:(num_pairs-2)
            z_1 = Zn(2*j+2);
            z_2 = Zn(2*j+3);
            z_inds = [2*j+2 2*j+3];
            store_pair();
        end
    end
    
    function store_pair()
        if max(z_1,z_2)>=0 && min(z_1,z_2)<=2*pi*(1+eps)
            Z_{count} = [z_1 z_2];
            Z_inds{count} = z_inds;
            count = count + 1;
        end 
    end
end

