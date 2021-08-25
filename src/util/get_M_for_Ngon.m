function M = get_M_for_Ngon(N)
%computes M, the number of reqired incident angles, for a regular N-sided
%polygon
    if N==1 %ambiguous meaning - but is interpreted as screen here
        M = 2;
    elseif mod(N,2)==0
        M = N^2/2;
    else
        M = N*(N+1);
    end
end