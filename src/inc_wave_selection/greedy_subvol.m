function [S, column_indices] = greedy_subvol(A)
    %vectorised version of greedy_subvol
    [M,~] = size(A);
    S = zeros(M,M);
    column_indices = zeros(1,M);
    
    Atemp = A;
    for m=1:M
        
        [max_norm, max_n] = max(sqrt(sum(Atemp.^2)));
        column_indices(m) = max_n;
        % fill in m'th row of S
        S(:,m) = A(:,max_n);
        v_max = Atemp(:,max_n);
        
        Atemp = Atemp - v_max*v_max.'*Atemp/(max_norm)^2;
    end
end

