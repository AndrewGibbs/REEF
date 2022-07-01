function [S, column_indices] = greedy_subvol_mod(A,keeps)
    %vectorised version of greedy_subvol, will keep the first N columns as
    %is and choose near optimal remaining columns
    [M,~] = size(A);
    S = zeros(M,M);
    column_indices = zeros(1,M);
    
    Atemp = A;
    %just ignore the other columns for now
    for m=keeps
        [max_norm, max_n_] = max(sqrt(sum(Atemp(keeps,keeps).^2)));
        max_n = keeps(max_n_);
        column_indices(m) = max_n;
        % fill in m'th row of S
        S(:,m) = A(:,max_n);
        v_max = Atemp(:,max_n);
        
        Atemp = Atemp - v_max*v_max.'*Atemp/(max_norm)^2;
    end
    
    for m=1:M-length(keeps)
%         if ~ismember(m,keeps)
            [max_norm, max_n] = max(sqrt(sum(Atemp.^2)));
%         else
%             max_n = m;
%             max_norm = max(sqrt(sum(Atemp(max_n,:).^2)));
%         end
        column_indices(m) = max_n;
        % fill in m'th row of S
        S(:,m) = A(:,max_n);
        v_max = Atemp(:,max_n);
        
        Atemp = Atemp - v_max*v_max.'*Atemp/(max_norm)^2;
    end
end

