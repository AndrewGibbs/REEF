function [M,V,p,q,name] = get_embedding_params(N)
    [X,~,name] = eg_bank(N);
    M = X.M;
    p = X.p;
    q = X.q;
    V = X.V;
end