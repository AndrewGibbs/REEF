%% For this example to work, you must have downloaded HNABEMLAB
% and added it to the matlab search path, see:
% https://github.com/AndrewGibbs/HNABEMLAB
clc;
clear all;
close all;

addpath("../");
addPathsHNA;
addPathsReef;

%% This code demonstrates Residue Enhanced Embedding Formulae 
% on a regular polygon. A screen is considered two-sided.
kwave=15; % wavenumber
num_sides = 4; % sides of polygon
num_tests = 1000;

%% set approximation parameters:
HNApolydeg = 5; %HNABEMLAB polnomial degree
hybrid_basis = true; % hybrid or standard BEM basis
OS = 1.5; % take a few extra incident angles to be safe
[M,V,p] = get_embedding_params(num_sides);

%% get canonical far-field patterns for M_ incident angles
M_ = ceil(M*OS);
alpha_in = linspace(0,2*pi,M_+1);
alpha_in = alpha_in(1:(end-1));

d = @(a) [cos(a+pi) sin(a+pi)];

D = HNAwrapper(V,kwave,alpha_in,HNApolydeg,hybrid_basis);

%% Use REEF to compute far field for large number of incident angles
E = Reef(D,alpha_in,kwave,p,'M',M); %create the embedding object

inc_test = linspace(0,2*pi,num_tests+1);
inc_test = inc_test(1:(end-1));
obs_test = linspace(0,2*pi,num_tests).';

Eout = E.getFarField(obs_test,inc_test); % get far-field data

% Plot the output
imagesc(obs_test,inc_test, abs(Eout));
xlabel('Inc angle $\alpha$','Interpreter','latex');
ylabel('Obs angle $\theta$','Interpreter','latex');
title('Far-field cross-section $|D(\theta,\alpha)|$','Interpreter','latex');
set(gca,'FontSize',20);