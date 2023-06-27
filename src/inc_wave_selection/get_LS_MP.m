function [alpha_in, ffDhat, FF_in, pinverse] = get_LS_MP(alpha_in,FF_in,hat,ffDhat,tol)

    ffDhat = get_sampling_matrix(alpha_in,FF_in,hat);
    pinverse = pinv(ffDhat,tol);
    
end

