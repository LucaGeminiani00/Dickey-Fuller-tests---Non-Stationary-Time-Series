function[x,epsilon] = for_arma(T,sigma,phi,theta)

epsilon=sqrt(sigma) * randn(T, 1);
x=zeros(T,1);
x(1)=epsilon(1);
Q=length(theta);
P=length(phi);
A=zeros(P,1);
M=zeros(Q,1);

for t=2:T
    for p=1:P
        if t-1>=p
            A(p)=phi(p)*x(t-p);
        else
            A(p)=0;
        end
    end
    for q=1:Q
        if t-1>=q
            M(q)=theta(q)*epsilon(t-q);
        else
            M(q)=0;
        end
    end
    x(t)=sum(A)+sum(M)+epsilon(t);
end
end
