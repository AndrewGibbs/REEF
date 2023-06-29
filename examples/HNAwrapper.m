function [FF, v_N, GOA,Gamma] = HNAwrapper(V,kwave,alphas,pMax,hybrid,derivs)
%     if nargin == 4
%         hybrid = false;
%     end
    if nargin == 5
        derivs = 0;
    end
    Vdims = size(V);
    if Vdims(1)>2
        Gamma=ConvexPolygon(V);    
        % construct the combined layer potential operator ---------------------------
        A = combinedLayer(kwave,Gamma);
    else
        Gamma=Screen(V);
        A = singleLayer(kwave,Gamma);
    end

    num_incs = length(alphas);
    %inident plane wave -------------------------------------------------------
    
    for n_=1:num_incs
        d = [cos(alphas(n_)+pi) sin(alphas(n_)+pi)];
        uinc{n_} = planeWave(kwave,d);
    end

    %make an HNA basis on Gamma -----------------------------------------------
    cL = 2; %layers of grading per polynomial degree
    sigmaGrad=0.15; %grading ratio
    nLayers = cL*(pMax+1)-1; %number of layers of grading
    OverSample = 1.25; %choose amount to oversample by
    hMax = 2*pi/(2*kwave);
%     qppw = 15;
    % construct the HNA basis (single mesh):
    if hybrid
        V = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
    else
        V = hpStandardBasis(Gamma, pMax, hMax, nLayers, sigmaGrad);
    end

    %solve
    [v_N,GOA] = ColHNA(A, V, uinc, Gamma,'oversample', OverSample,'symmetry');
    
    %compute derivatives if requested to do so
    if derivs>0
        for j = 0:derivs
            if hybrid
                for n=1:num_incs
                    FF{n,j+1} = @(theta) FarField(Gamma, v_N{n}, kwave, theta) ...
                        + FarField(Gamma, GOA{n}, kwave, theta);
                end
            else
                for n=1:num_incs
                    FF{n,j+1} = @(theta) FarField(Gamma, v_N{n}, kwave, theta);
                end
            end
        end
    else
        if hybrid
            for n=1:num_incs
                FF{n} = @(theta) FarField(Gamma, v_N{n}, kwave, theta) ...
                    + FarField(Gamma, GOA{n}, kwave, theta);
            end
        else
            for n=1:num_incs
                FF{n} = @(theta) FarField(Gamma, v_N{n}, kwave, theta);
            end
        end
    end
end

