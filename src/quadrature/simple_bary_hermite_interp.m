function p = simple_bary_hermite_interp(x0,x1,f0,f1,x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    x_m_x0_squared = (x-x0)^2;
    x_m_x1_squared = (x-x1)^2;
    x1_m_x0_squared = (x1-x0)^2;

    p_top = (f0./x_m_x0_squared + f1./x_m_x1_squared)/x1_m_x0_squared;
    p_bottom = (1./x_m_x0_squared + 1./x_m_x1_squared)/x1_m_x0_squared;

    p = p_top./p_bottom;
end