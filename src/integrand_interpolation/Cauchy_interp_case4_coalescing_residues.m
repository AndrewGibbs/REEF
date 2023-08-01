function I = Cauchy_interp_case4_coalescing_residues(theta_set, residue, p_at_resisude, dp_at_residues, denom, h, Nquad)

    theta_set = theta_set(:);

    % get focii of ellipse (sort of)
    a = (residue)-h;
    b = (residue)+h;

    [z,w] = Cauchy_box_quad(theta_set, a, b, h, Nquad, 1, false);
    % can be made more efficient and looped over
    P1 = dp_at_residues; % gradient
    P2 = p_at_resisude-P1*residue; % intercept
    P = P1*z + P2;
    I = (w.')*(P./denom(z));
end