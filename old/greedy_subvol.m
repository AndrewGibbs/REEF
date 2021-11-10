function S = greedy_subvol(A)
    [M,N] = size(A);
    S = [];
    
    Atemp = A;
    for m=1:M
        
        %get max norm
        max_norm = 0;
        for n=1:N
            norm_v = norm(Atemp(:,n));
            if norm_v > max_norm
                max_n = n;
                max_norm = norm_v;
            end
        end
        
        % join v to S
        S = [S A(:,max_n)];
        v_max = Atemp(:,max_n);
        
        %now subtract projection from A
        for n=1:N
            Atemp(:,n) = Atemp(:,n) - proj(v_max,Atemp(:,n));
        end
    end
    
    function I = proj(u,v)
        I = u*((u.'*v)/(u.'*u));
    end
end

