%{
N=100;
f=zeros(N,N);
X=zeros(N,N);
Y=zeros(N,N);

ff=@(x,y) 100*(y-x^2)^2+(1-x)^2;
gf=@(x,y) [-400*x*(y-x^2)-2*(1-x);200*(y-x^2)]; 
Hf=@(x,y) [800*y^2-400*(y-x^2)+2 -400*x;-400*x 200];

for j=1:N
    for k=1:N
        x=[5*(j-1)/(N-1)-1.001;5*(k-1)/(N-1)+0.001];
        X(j,k)=x(1);
        Y(j,k)=x(2);
    end
end
F=100*(Y-X.^2)^2+(1-X)^2;
contourf(X,Y,F,15)
%}

%------------------------BFGS------------------------%
function []=BFGS(a, b)

ff=@(x,y) 100*(y-x^2)^2+(1-x)^2;
gf=@(x,y) [-400*x*(y-x^2)-2*(1-x);200*(y-x^2)]; 
Hf=@(x,y) [800*y^2-400*(y-x^2)+2 -400*x;-400*x 200];

x0=[a;b];
%hold on;
%plot(x0(1),x0(2),'ok','MarkerFaceColor','k','MarkerSize',4)
w=0.8;
c=0.1;
Bk=inv(Hf(x0(1),x0(2)));
determinant = det(Bk)
for j=1:15
    pk=-(Bk*gf(x0(1),x0(2)));
    alpha=1;
    while x0(2)+alpha*pk(2)<=0
        alpha=alpha*w;
    end
    nf=ff(x0(1)+alpha*pk(1),x0(2)+alpha*pk(2));
    dd=gf(x0(1),x0(2))'*pk;
    while (nf>ff(x0(1),x0(2))+c*dd*alpha)
        alpha=alpha*w;
        nf=ff(x0(1)+alpha*pk(1),x0(2)+alpha*pk(2));
    end
    alpha
    sk=alpha*pk;
    x1=x0+sk;
    %plot(x1(1),x1(2),'ok','MarkerFaceColor','k','MarkerSize',4)
    %plot([x0(1);x1(1)],[x0(2);x1(2)],'k','LineWidth',2)
    yk=gf(x1(1),x1(2))-gf(x0(1),x0(2));
    Bk= Bk-(Bk*(sk*sk')*Bk')/(sk'*Bk*sk)+(yk*yk')/(yk'*sk);
    x0=x1
    yk'*sk
    pause
end
%-----------------------------------------------------%