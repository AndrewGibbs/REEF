%% For this example to work, you must have downloaded HNABEMLAB and added it
% to the matlab search path, see:
% https://github.com/AndrewGibbs/HNABEMLAB
clc;
clear all;

addPathsHNA;
addPathsReef;

%% This code demonstrates Residue Enhanced Embedding Formulae 
% on a regular polygon. A screen is considered two-sided.
kwave=15; % wavenumber
num_sides = 4; % sides of polygon
num_tests = 1000;

%% set approximation parameters:
HNApolydeg = 5;
hybrid_basis = true;
[M,V,p] = get_embedding_params(2);
OS = 1.5; % take a few extra incident angles to be safe

%% get canonical far-field patterns for M_ incident angles
M_ = ceil(M*OS);
alpha_in = linspace(0,2*pi,M_+1);
alpha_in = alpha_in(1:(end-1));

alphaTest = linspace(0,2*pi,num_tests+1);
alphaTest = alphaTest(1:(end-1));
d = @(a) [cos(a+pi) sin(a+pi)];

D = HNAwrapper(V,kwave,alpha_in,HNApolydeg,hybrid_basis);

%% Use REEF to compute far field for large number of incident angles
obs_test = linspace(0,2*pi,num_tests).';
E = Reef(D,alpha_in,kwave,p,'M',M); %create the embedding object
Eout = E.getFarField(obs_test,inc_test);

% Plot the output
imagesc(obs_test,inc_test, abs(Eout));
xlabel('Inc angle $\alpha$','Interpreter','latex');
ylabel('Obs angle $\theta$','Interpreter','latex');
title('Far-field cross-section $D(\theta,\alpha)$','Interpreter','latex');