function z_out = tweakBadAngles(z,thresh)

    [zs,I] = sort(z);
    
    for n=1:(length(z)-1)
        if abs(zs(n)-zs(n+1))<thresh
           zs(n) = zs(n) - thresh/2;
           zs(n+1) = zs(n+1) - thresh/2;
        end
    end
    if abs(zs(1)-zs(end)-2*pi)<thresh
        zs(1) = zs(1) + thresh/2;
        zs(end) = zs(end) + thresh/2;
    end
    
    %resort any weird order there was before
    z_out(I) = zs;

end