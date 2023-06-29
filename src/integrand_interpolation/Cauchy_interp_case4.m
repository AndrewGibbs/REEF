function I = Cauchy_interp2(theta_set, residues, p_at_resisudes, denom, h, Nquad)
    % computes cauchy integral for the case when theta is not inside
    % contour
    theta_set = theta_set(:);

    % get focii of ellipse (sort of)
    a = min(residues)-h;
    b = max(residues)+h;

    [z,w] = Cauchy_box_quad(theta_set, a, b, h, Nquad, 1, false);
    % can be made more efficient and looped over
    interp_poly = barylag([residues p_at_resisudes],z);
    I = (w.')*(interp_poly./denom(z));
    
end