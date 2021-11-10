function MPS_FarField_polygon()

%-------------------------------------------
% example of computing the far field using
% MPSPACK.
%
% Copyright (C) 2014 Stuart C. Hawkins
%-------------------------------------------

clear all;

%-------------------------------------------
% set main parameters
%-------------------------------------------

% wavenumber
kwave = 10;
n = 50;
tau = 5e-2;
m = 2*n;
cornern = n;
direction = 0;

% setup structure with the parameters for the MFS method
opts= struct('eta',kwave,'fast',2,'nmultiplier',2,'tau',tau);

% setup structure with parameters for the corner basis
% functions
nuopts=struct('type','s','cornermultipliers',[0 0 1 0 0],'rescale_rad',1);

% set number of vertices
V = 6;

% set radius of polygon
r0 = 1;

% set radius of bounding circle
r1 = 2;

% compute polar angles of the radial edges of the artifical
% domains
t = 2*pi*((0:V-1)-0.5)/V;

% compute corners of artificial domains that lie on the
% bounding circle
r = r1 * exp(1i*t);

% compute the corners of the polygon 
a = r0 * exp(1i*2*pi*(0:V)/V);

% compute the midpoints of the faces
c = 0.5*(a(2:end)+a(1:end-1));

% set the radius for the MFS basis
R = 0.8*r1;

% create straight line segments that are used to construct the 
% first artifical domain 
s = segment.polyseglist(m,[r(2) c(1) a(1) c(V) r(1)]);

% add the boundary that lies on the bounding circle
s = [s(1:3) segment(3*m,[0 r1 t(1) t(2)])];

% create a temporary copy that we can rotate to make the other
% artificial domains
tmp = [s(1:3) segment(3*m,[0 r1 t(1) t(2)])];

% rotate the temporary copy to make the other artificial
% domains
for k=1:V-1
    s = [s rotate(tmp,2*pi*k/V)];
end

% collect segments that make up all the  artificial boundaries
sart = s(interlace(1+4*(0:V-1),4+4*(0:V-1)));

% collect segments that make up the artificial outer boundary
sext = s([4+4*(V-1:-1:0)]);

% collect segments that make the polygon
polygon = s(interlace(2+4*(V-1:-1:0),3+4*(V-1:-1:0)));

% setup domains for the artifical domains
for k=1:V
    ii = [1 2 3 -3 4] + 4*(k-1);
    ii = mod(ii-1,4*V)+1;
    d(k)=domain(s(ii),[1 1 1 -1 1]);
end

% setup a domain for the polygon
dpoly = domain([],[],polygon,1);

% setup the exterior domain
ext = domain([], [], sext, -1);

% set transmission BC on the artificial boundaries
sart.setmatch([-kwave kwave],[1 -1]);

% add the corner basis functions
for k=1:V
    d(k).addcornerbases(cornern,nuopts);
end

% setup MFS basis for the exterior
Z=@(t) R*exp(2i*pi*t);
Zp=@(t) 2i*pi*R*exp(2i*pi*t);            
ext.addmfsbasis({Z,Zp},n,opts);

% initialise the scattering problem
scatteringObject = scattering(ext,d);
scatteringObject.setoverallwavenumber(kwave);
scatteringObject.setincidentwave(direction);
scatteringObject.fillquadwei();
scatteringObject.rhs(:,k) = scatteringObject.fillrighthandside();
scatteringObject.solvecoeffs;
%coeffs{1} = scatteringObject.co(:,1);

%-------------------------------------------
% plot farfield
%-------------------------------------------

%clear obj
p = scatteringObject;

% setup structure describing the mesh... 
% dx is theta spacing, I think
obj = struct('dx', 0.01);

getFarField(self,points,index, dn)
end
