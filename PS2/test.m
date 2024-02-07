N=100;
f=zeros(N,N);
X=zeros(N,N);
Y=zeros(N,N);

ff=@(x,y) 4*x^2-2*x*y+3*y^2-log(y)+exp(x-3);
gf=@(x,y) [8*x-2*y+exp(x-3);-2*x+6*y-1/y]; 
Hf=@(x,y) [8+exp(x-3), -2; -2, 6+1/(y^2)];

for j=1:N
    for k=1:N
        x=[5*(j-1)/(N-1)-1.001;5*(k-1)/(N-1)+0.001];
        X(j,k)=x(1);
        Y(j,k)=x(2);
    end
end
F=4*X.^2-2*X.*Y+3*Y.^2-log(Y)+exp(X-3);
contourf(X,Y,F,15)

x0=[3;4];
hold on;
plot(x0(1),x0(2),'ok','MarkerFaceColor','k','MarkerSize',4)
w=0.8;
c=0.01;
for j=1:100
    %spk=-gf(x0(1),x0(2));
    %npk=-(inv(Hf(x0(1),x0(2)))*gf(x0(1),x0(2)));
    pk=-(inv(Hf(x0(1),x0(2)))*gf(x0(1),x0(2)));
    alpha=5;
    while x0(2)+alpha*pk(2)<=0
        alpha=alpha*w;
    end
    nf=ff(x0(1)+alpha*pk(1),x0(2)+alpha*pk(2));
    dd=gf(x0(1),x0(2))'*pk;
    while (nf>ff(x0(1),x0(2))+c*dd*alpha)
        alpha=alpha*w;
        nf=ff(x0(1)+alpha*pk(1),x0(2)+alpha*pk(2));
    end
    Hf
    nf
    x1=x0+alpha*pk
    plot(x1(1),x1(2),'ok','MarkerFaceColor','k','MarkerSize',4)
    plot([x0(1);x1(1)],[x0(2);x1(2)],'k','LineWidth',2)
    x0=x1;
    pause
end