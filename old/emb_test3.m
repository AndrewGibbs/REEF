a = 1;
%p = {@(x) (-1i-1)*lim+x, @(x) (-1i+1)*lim+x*1i , @(x) (1i-1)*lim+x, @(x) (-1+1i)*lim-1i*x};
% p{2} = @(x) (-1i+1)*lim+x*1i;
% p{3} = @(x) (1i+1)*lim-x;
% p{4} = @(x) (-1+1i)*lim-1i*x;
embed_std = @(theta) (B(1)*Dhat1(theta.') + B(2)*Dhat2(theta.'))./(hat(theta,alphaTest));

I(1) = integral(@(x) embed_std(x)./(x-a),1i+2*pi+pi/2,1i-pi/2,'ArrayValued',true);%top;
I(2) = integral(@(x) embed_std(x)./(x-a),-1i-pi/2,-1i+2*pi+pi/2,'ArrayValued',true);%bottom
I(3) = integral(@(x) embed_std(x)./(x-a),1i-pi/2,-1i-pi/2,'ArrayValued',true);%left;
I(4) = integral(@(x) embed_std(x)./(x-a),-1i+2*pi+pi/2,1i+2*pi+pi/2,'ArrayValued',true);%right

I = I.';
sum(I)/(2i*pi)

return;
Cauchy_embed_FF_vec = @(theta) [w_x*sum(-(B(1)*top_integral_ff{1}+B(2)*top_integral_ff{2})./(((top_line-theta)).*(hat(top_line,alphaTest))));
                        w_x*sum((B(1)*bot_integral_ff{1}+B(2)*bot_integral_ff{2})./((bot_line-theta).*(hat(bot_line,alphaTest))));
                        w_y*1i*sum(-(B(1)*left_integral_ff{1}+B(2)*left_integral_ff{2})./((left_line-theta).*(hat(left_line,alphaTest))));
                        w_y*1i*sum((B(1)*right_integral_ff{1}+B(2)*right_integral_ff{2})./((right_line-theta).*(hat(right_line,alphaTest))))];
J = Cauchy_embed_FF_vec(a);
sum(J)/(2i*pi)

actual_farfield(a)

Cauchy_embed_FF_vec2 = @(theta) [(-w_x*(1./(top_line.'-theta))*((B(1)*top_integral_ff{1}+B(2)*top_integral_ff{2})./(hat(top_line,alphaTest)))).';
                        (w_x*(1./(bot_line.'-theta))*((B(1)*bot_integral_ff{1}+B(2)*bot_integral_ff{2})./(hat(bot_line,alphaTest)))).';
                        (-w_y*1i*(1./(left_line.'-theta))*((B(1)*left_integral_ff{1}+B(2)*left_integral_ff{2})./(hat(left_line,alphaTest)))).';
                        (w_y*1i*(1./(right_line.'-theta))*((B(1)*right_integral_ff{1}+B(2)*right_integral_ff{2})./(hat(right_line,alphaTest)))).'];
                    
J2 = @(theta) sum(Cauchy_embed_FF_vec2(theta))/(2i*pi);