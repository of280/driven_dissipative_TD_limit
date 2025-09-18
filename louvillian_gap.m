k=0:0.01:1; %%%VALUES OF kappa/g
s=1; %%VALUE OF s
N=3; %%VALUE OF N
lgap=zeros(1,length(k));
for i=1:length(k)
   lgap(i)=-real(dellgl(N,k(i),s));
end
plot(k,lgap)

xlabel('\kappa/g') 
ylabel('Liouvillian gap') 

function vp = dellgl(k, g, s)
    ii = zeros(7*k^4-12*k^3+6*k^2,1);
    jj = zeros(7*k^4-12*k^3+6*k^2,1);
    ss = zeros(7*k^4-12*k^3+6*k^2,1);
    c=0;
    for i=1:k^4
        d=intoba(i-1,k);
        n=d(1);
        m=d(2);
        p=d(3);
        q=d(4);
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
    end
    rr=sparse(ii,jj,ss,k^4,k^4);
   % sm6=eig(full(rr));
   sm6=eigs(rr,10,'largestreal');
  if isnan(sm6(2))
     sm7=eigs(rr,20,'largestreal');
     vp=sm7(2);
     if isnan(sm7(2))
          sm8=eigs(rr,30,'largestreal');
          vp=sm8(2);
     else
          vp=sm7(2);
     end
  else
      vp=sm6(2);
  end
end 
function base = intoba(k, n)
    r1=floor(k/n^3);
    r2=floor((k-r1*n^3)/n^2);
    r3=floor((k-r1*n^3-r2*n^2)/n);
    r4=k-r1*n^3-r2*n^2-r3*n;
    base = [r4 r3 r2 r1];
end