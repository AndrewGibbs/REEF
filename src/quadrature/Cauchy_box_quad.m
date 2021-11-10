function [z,W] = Cauchy_box_quad(x,a,b,h,qppw,kwave,div1)
%returns quadrature for a Cauchy integral, to evaluate a function at the
%points x
if nargin == 6
    %use the divide by one trick
    div1 = true;
end
    x = x(:);
    [ zz_top, w_top ] = gauss_quad_wave_split3(b+1i*h, a+1i*h, qppw, kwave);
    [ zz_left, w_left ] = gauss_quad_wave_split3(a+1i*h, a-1i*h, qppw, kwave);
    [ zz_bottom, w_bottom ] = gauss_quad_wave_split3(a-1i*h, b-1i*h, qppw, kwave);
    [ zz_right, w_right ] = gauss_quad_wave_split3(b-1i*h, b+1i*h, qppw, kwave);
    
    z = [zz_top; zz_left; zz_bottom; zz_right];
    w = [w_top; w_left; w_bottom; w_right]/(2i*pi);
    W = w./(z-x.');
    if div1
        W = W./sum(W); % divide by one
    end
end