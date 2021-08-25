function canon_angles = get_canon_angles_for_Ngon(N,num_guesses)
%uses GTD to select canonical incident angles which are less
%likely to lead to an ill-conditioned sampling matrix
if nargin == 1
    num_guesses = 1000;
end
    %quickly try all of these inc angles, and select the best ones:
    alpha_guesses = linspace(0, 2*pi, num_guesses);
    alpha_guesses = alpha_guesses(1:(end-1));
    % get matrix of coefficients of local expansion at each corner:
    A_wide = getAmat(alpha_guesses, N);
    % now choose submatrix which has (close to) maximum determinant:
    [~, column_indices] = greedy_subvol(A_wide);
    % select incident angles corresponding to maximised determinant:
    canon_angles = alpha_guesses(column_indices);
    
end

