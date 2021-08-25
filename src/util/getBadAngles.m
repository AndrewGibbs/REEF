function [ Z1, Z2 ] = getBadAngles(alpha_out,p,RotShift)
    if nargin ==2
        RotShift = 0;
    end
    n_count=1;
    if mod(p,2)==0 %p is even
        for n=0:(p-1)
           Z1(:,n_count)=alpha_out + 2*pi*n/p; 
           n_count=n_count+1;
           Z1(:,n_count)=2*pi -alpha_out + 2*pi*n/p - 2*RotShift; 
           n_count=n_count+1;
        end
    else %p is odd
        for n=0:(p-1)
           Z1(:,n_count)=alpha_out+(2*n+1)*pi/p; 
           n_count=n_count+1;
           Z1(:,n_count)=2*pi - alpha_out+(2*n+1)*pi/p - 2*RotShift; 
           n_count=n_count+1;
        end
    end
    Z2=repmat((0:pi/p:(2*pi)),length(alpha_out),1)-RotShift;

    Z1=(mod(Z1,2*pi));
    Z2=(mod(Z2,2*pi));
end