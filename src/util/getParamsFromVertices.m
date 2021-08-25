function [p,M,RotShift,q] = getParamsFromVertices(vertices)
    thresh1 = 1e-8;
    pMax = 50;
    qMax = 100;
    [is_regular, intAngles]=angleTest(vertices, thresh1);
    N = length(intAngles);
    RotShift = acos((vertices(1,:)- vertices(2,:))*[1;0]/norm(vertices(2,:)- vertices(1,:)));
    if is_regular
        p=N; q=N+2;
        %now divide top and bottom to make sure p and q are least
        %integers possibe (coprime)
        GCD=gcd(p,q);
        while GCD>1
            p=p/GCD;
            q=q/GCD;
            GCD=gcd(p,q);
        end

        M=(q-1)*N;

    else %irregular
        for j=1:N
           [qIntl(j),pIntl(j)]=qpIrregPolygon( intAngles(j), pMax, qMax );
        end         
        pIntlAll=LCM2(pIntl);
        qIntAll=round(min(pIntl./qIntl)/pIntlAll);
        %now re-arrange to get external angles:
        p=pIntlAll;
        qExtAll=2*pIntlAll -qIntAll; %this is probably quite irrelevant

        %now divide top and bottom to make sure p and q are least
        %integers possibe (coprime)
        GCD=gcd(p,qExtAll);
        while GCD>1
            p=p/GCD;
            qExtAll=qExtAll/GCD;
            GCD=gcd(p,qExtAll);
        end

        for j=1:N
            q_j(j)=round((2*pi-intAngles(j))*p/pi);
            m_j(j)=q_j(j)-1;
        end
        M=sum(m_j);
        q = q_j;
    end
end

