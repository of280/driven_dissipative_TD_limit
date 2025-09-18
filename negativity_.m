%%% USING qubit4matlab found at https://www.gtoth.eu/qubit4matlab.html

k=(0:0.01:1);  %%%VALUE OF kappa/g
neg=zeros(1,length(k));
N=3; %%VALUE OF N
s=1; %%VALUE OF S
for i=1:length(k)
   neg(i)=negatgl(N,k(i),s);
end
plot(k,neg)

xlabel('\kappa/g') 
ylabel('Negativity') 

function neg = negatgl(k, g, s)
    ii = zeros(k^4+4*k^3-5*k^2+2*k,1);
    jj = zeros(k^4+4*k^3-5*k^2+2*k,1);
    ss = zeros(k^4+4*k^3-5*k^2+2*k,1);
    b = zeros(k^4+1,1);
    b(k^4+1) = 1;
    c=0;
    for i=1:k^4
        d=intoba(i-1,k);
        n=d(1);
        m=d(2);
        p=d(3);
        q=d(4);
        if n+m==p+q
            c=c+1;
            ii(c)=i;
            jj(c)=i;
            ss(c)=-(s*heaviside(k-n-1.1)*(n+1)+s*heaviside(k-p-1.1)*(p+1)+(m+q))/2;
            if m~=k-1
                if q~=k-1
                    c=c+1;
                    ii(c)=i;
                    jj(c)=i+k+k^3;
                    ss(c)=sqrt((m+1)*(q+1));
                end
            end
            if n~=0
                if p~=0
                    c=c+1;
                    ii(c)=i;
                    jj(c)=i-1-k^2;
                    ss(c)=s*sqrt(n*p);
                end
            end
            if p~=0
                if q~=k-1
                    c=c+1;
                    ii(c)=i;
                    jj(c)=i-k^2+k^3;
                    ss(c)=1j*g*sqrt(p*(q+1));
                end
            end
            if p~=k-1
                if q~=0
                    c=c+1;
                    ii(c)=i;
                    jj(c)=i+k^2-k^3;
                    ss(c)=1j*g*sqrt(q*(p+1));
                end
            end
            if n~=0
                if m~=k-1
                    c=c+1;
                    ii(c)=i;
                    jj(c)=i+k-1;                    
                    ss(c)=-1j*g*sqrt(n*(m+1));
                end
            end
            if n~=k-1
                if m~=0
                    c=c+1;
                    ii(c)=i;
                    jj(c)=i-k+1;                    
                    ss(c)=-1j*g*sqrt(m*(n+1));
                end
            end
            if n==p
                if m==q
                    c=c+1;
                    ii(c)=k^4+1;
                    jj(c)=i;                    
                    ss(c)=1;
                end
            end
        else
            c=c+1;
            ii(c)=i;
            jj(c)=i;            
            ss(c)=1;
        end
    end
    rr=sparse(ii,jj,ss,k^4+1,k^4);
    sol=rr\b;
    rho=zeros(k^2);
    for i=1:k^4
        if sol(i) ~= 0
            d2=intoba(i-1, k);
            rho(d2(1)+k*d2(2)+1,1+d2(3)+k*d2(4))=sol(i);
        end
    end
   neg = negativity(rho,1,k);
end 

function base = intoba(k, n)
    r1=floor(k/n^3);
    r2=floor((k-r1*n^3)/n^2);
    r3=floor((k-r1*n^3-r2*n^2)/n);
    r4=k-r1*n^3-r2*n^2-r3*n;
    base = [r4 r3 r2 r1];
end