function [alpha_in, ffDhat,FF_in,full_rank] = choose_submatrix(alpha_in,FF_in,hat,M,tol)
ffDhat_temp = zeros(length(alpha_in));
    for m = 1:length(alpha_in)
        ffDhat_temp(:,m) = hat(alpha_in,alpha_in(m)).*FF_in{m}(alpha_in);
    end

    % the next several lines are based on mathworks file exchange code at:
    %https://uk.mathworks.com/matlabcentral/fileexchange/77437-extract-linearly-independent-subset-of-matrix-columns
%    tol = self.rank_tol;
   if max(max(abs(ffDhat_temp)))<tol %check if all entries are almost zero
       r = 0;
       E = 1:M;
   else
       [~, R, E] = qr(ffDhat_temp,0); 
       if ~isvector(R)
        diagr = abs(diag(R));
       else
        diagr = R(1);   
       end
       %Rank estimation
       r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation, should be <=M
       uniqe_m=sort(E(1:M));
   end

    if r < M
       warning(sprintf('%dx%d collocation matrix has rank %d - results will be highly inaccurate - please consider different/more indicent angles',M,M,r));
       full_rank = false;
       ffDhat = ffDhat_temp(1:M,1:M);
       uniqe_m_trunc = 1:M;
    else
       uniqe_m_trunc = uniqe_m(1:M);
       alpha_in = alpha_in(uniqe_m_trunc);
       full_rank = true;
       ffDhat = ffDhat_temp(uniqe_m_trunc,uniqe_m_trunc);
    end
    
    for m=1:M
        FF_in{m} = FF_in{uniqe_m_trunc(m)};
    end
end

