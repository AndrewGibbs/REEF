% function tf=angleTest(V,T)
%     [N,~]=size(V);
%     A=zeros(1,N);
%     for n=1:(N-1)
%         A(n)=[ V(n,1) V(n,2) ]*[ V(n+1,1); V(n+1,2) ];
%     end
%     A(N)=[ V(N,1) V(N,2)]*[ V(1,1); V(1,2) ];
%     if min(abs(mean(A)-A))<abs(T)
%         tf=true;
%     else
%         tf=false;
%     end
% end
function [tf, A]=angleTest(V,T)
    [N,~]=size(V);
    A=zeros(1,N);
    side_mod=@(n,N) mod(n-1,N)+1;
    for n=1:N
        %A(n)=( V(side_mod(n+1,N),:)- V(n,:) )*( V(side_mod(n+2,N),:)- V(side_mod(n+1,N),:) ).';
         A(n)=( V(side_mod(n,N),:)- V(side_mod(n-1,N),:) )*( V(side_mod(n,N),:)- V(side_mod(n+1,N),:) ).'/(norm(V(side_mod(n,N),:)- V(side_mod(n-1,N),:))*norm(V(side_mod(n,N),:)- V(side_mod(n+1,N),:)));
    end
    A=acos(A);
    %A(N)=[ V(N,1) V(N,2)]*[ V(1,1); V(1,2) ];
    if min(abs(mean(A)-A))<abs(T)
        tf=true;
    else
        tf=false;
    end
end