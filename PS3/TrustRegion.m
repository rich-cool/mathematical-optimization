%------------------Newton Method------------------%
function []=NewtonMethod(a, b)

N=100;
f=zeros(N,N);
X=zeros(N,N);
Y=zeros(N,N);

ff=@(x,y) 100*(y-x.^2).^2+(1-x).^2;
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

x0=[a;b];

hold on;
plot(x0(1),x0(2),'ok','MarkerFaceColor','k','MarkerSize',4)
p=inv(Hf(x0(1),x0(2)))*(-gf(x0(1),x0(2)));
m=@(p) ff(x0(1),x0(2))+gf(x0(1),x0(2))'*p+.5*p'*Hf(x0(1),x0(2))*p;
ff(x0(1),x0(2))-ff(x0(1)+p,x0(2)+p)
inv(m(0)-m(p))
row=inv((m(0)-m(p)))*(ff(x0(1),x0(2))-ff(x0(1)+p,x0(2)+p));
w=.25;
delta_max = 4;
n = .1;
for j=1:10
    %{
     while x0(2)+alpha*pk(2)<=0
        alpha=alpha*w;
     end
    %}
    x1=x0;
    delta = 1;
    while row<w
        delta=delta*w;
    end
    
    while (row>.75 & norm(p)==delta)
        delta = min(2*delta,delta_max);
    end
    
    if row>n
        x1=x0+p;
    end
    plot(x1(1),x1(2),'ok','MarkerFaceColor','k','MarkerSize',4)
    plot([x0(1);x1(1)],[x0(2);x1(2)],'k','LineWidth',2)
    x0=x1
    pause
end
%---------------------------------------------------%