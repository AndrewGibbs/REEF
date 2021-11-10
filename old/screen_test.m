clear all;
clc;
addpath '../../Dropbox/HNABEMLAB';
addPathsHNA;
kwave=1;
d = [1 -1]./sqrt(2); %direction as a vector
qppw = 500;
imag_dist = .5/kwave;
Nquad_x = qppw*kwave;
Nquad_y = ceil(qppw*imag_dist/(2*pi/kwave));

[v_N, GOA, Gamma] = getSolScreen(kwave,d);

%get all the quadrature points:
theta_quad = linspace(0,2*pi,Nquad_x);
w_x = 2*pi/(Nquad_x-1);
theta_quad = theta_quad(1:(end-1))+w_x/2;

vert_line = 1i*linspace(-imag_dist, imag_dist, Nquad_y);
w_y = 2*imag_dist/(Nquad_y-1);
vert_line = vert_line(1:(end-1))+1i*w_y/2;

top_line = theta_quad+imag_dist*1i;
bot_line = theta_quad-imag_dist*1i;
left_line = vert_line;
right_line = 2*pi + vert_line;

top_integral_ff = getFarFieldScreen(top_line, v_N, GOA, Gamma, kwave).';
bot_integral_ff = getFarFieldScreen(bot_line, v_N, GOA, Gamma, kwave).';
left_integral_ff = getFarFieldScreen(left_line, v_N, GOA, Gamma, kwave).';
right_integral_ff = getFarFieldScreen(right_line, v_N, GOA, Gamma, kwave).';

% Construct Cauchy integral version of the far-field pattern:
Cauchy_FarField = @(theta_) arrayfun(@(theta) (w_x*(sum(-top_integral_ff./(top_line-theta)) ...
                         + sum(bot_integral_ff./(bot_line-theta)))...
                        +w_y*1i*(sum(-left_integral_ff./(left_line-theta)) ...
                         +sum(right_integral_ff./(right_line-theta)))...
                        )/(2i*pi),theta_);
                    
Ifun = @(theta) [w_x*sum(-top_integral_ff./(top_line-theta));
                    w_x*sum(bot_integral_ff./(bot_line-theta));
                    1i*w_y*sum(-left_integral_ff./(left_line-theta));
                    1i*w_y*sum(right_integral_ff./(right_line-theta))];
                     
actual_farfield = getFarFieldScreen(theta_quad, v_N, GOA, Gamma, kwave);
figure(1);
plot(theta_quad,Cauchy_FarField(theta_quad.'),theta_quad,actual_farfield,'.');
figure(2);
plot(theta_quad,imag(Cauchy_FarField(theta_quad.')),theta_quad,imag(actual_farfield),'.');