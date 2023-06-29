function a = interp_cubic(x0,x1,y0,y1,dy0)
    p = [y0;y1;dy0];
    V = [1 x0 x0^2;
        1 x1 x1^2;
        0 1 2*x0];
    a = V\p;
end