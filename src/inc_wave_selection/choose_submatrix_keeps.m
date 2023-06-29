function [alpha_in,ffDhat,FF_in,FF_derivs] = choose_submatrix_keeps(alpha_in,FF_in,FF_derivs,ffDhat_temp,hat,M,tol,keeps)
    
    % now choose submatrix, subject to the constraint of keeping ros/columns indexed by 'keeps'
    [~, column_indices] = greedy_subvol_mod(ffDhat_temp,keeps);
    
    % if the condition number is bad, try again, dropping one of the keeps,
    % (if there are any keeps specified)
    if cond(ffDhat_temp(FF_incolumn_indices,column_indices))>tol && ~isempty(keeps)
        warning(sprintf('%dx%d collocation matrix has condition number %.1e - trying again, with fewer specified inc angles',M,M,cond(FF_in)));
        keeps_drop_one = keeps(column_indices(1:(length(keeps)-1)));
        [alpha_in,ffDhat,FF_in,FF_derivs] = ...
            choose_submatrix_keeps(alpha_in,FF_in,FF_derivs,ffDhat_temp,hat,M,tol,keeps_drop_one);
    else
        % reduce the matrix
        alpha_in = alpha_in(column_indices);
        FF_in = FF_in(column_indices);
        FF_derivs = FF_derivs(column_indices);
        ffDhat = ffDhat_temp(column_indices,column_indices);
    end
    %final check / warning
    if cond(ffDhat)>tol
        warning(sprintf('%dx%d collocation matrix has condition number %.1e - results will be highly inaccurate - please consider different/more indicent angles',M,M,cond(FF_in)));
    end
    
end

