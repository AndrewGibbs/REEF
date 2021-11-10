function Dmat = get_D_matrix(M,D,x,w,hat,alpha_in)
    Dmat = zeros(length(w),M);
    for m=1:M
        Dmat(:,m) = hat(x,alpha_in(m)).*D{m}(x).*w;
    end
end

