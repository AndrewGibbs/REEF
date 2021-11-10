a = 1;
%p = {@(x) (-1i-1)*lim+x, @(x) (-1i+1)*lim+x*1i , @(x) (1i-1)*lim+x, @(x) (-1+1i)*lim-1i*x};
% p{2} = @(x) (-1i+1)*lim+x*1i;
% p{3} = @(x) (1i+1)*lim-x;
% p{4} = @(x) (-1+1i)*lim-1i*x;

I(1) = integral(@(x) cos(x)./(x-a),-1i,-1i+2*pi);%integral(@(x) cos(p{1}(x))./(p{1}(x)),0,2*lim);
I(2) = integral(@(x) cos(x)./(x-a),-1i+2*pi,1i+2*pi);%integral(@(x) cos(p{2}(x))./(p{2}(x)),0,2*lim);
I(3) = integral(@(x) cos(x)./(x-a),1i+2*pi,1i);%pintegral(@(x) cos(p{3}(x))./(p{3}(x)),0,2*lim);
I(4) = integral(@(x) cos(x)./(x-a),1i,-1i);%integral(@(x) cos(p{4}(x))./(p{4}(x)),0,2*lim);
I