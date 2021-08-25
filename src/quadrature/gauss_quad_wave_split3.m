function [ x, w ] = gauss_quad_wave_split3(a, b, qppw, k,  ab_width )
    if nargin==4
        ab_width=abs(b-a);
    end
    
    if ab_width == 0
        x = [];
        w = [];
        return;
    end
    
    wavelength = 2*pi/k;
    
    wavelengthsInWidth = ceil(ab_width/wavelength);
    
    sub_width = ab_width/wavelengthsInWidth;
    
    x=[]; w=[];
    
    dz = (b-a)/ab_width;
    
    for q=1:wavelengthsInWidth
        [x_,w_]=gauss_quad_complex(a+(q-1)*sub_width*dz,a+q*sub_width*dz,qppw);
        x=[x; x_];
        w=[w; w_*dz];
    end
end