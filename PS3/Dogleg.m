%------------------Newton Method------------------%
function []=Dogleg(a, b)

N=100;
X=zeros(N,N);
Y=zeros(N,N);

ff=@(x,y) 10*(y-x^2)^2+(1-x)^2;
gf=@(x,y) [-40*x*(y-x^2)-2*(1-x);20*(y-x^2)]; 
Hf=@(x,y) [80*y^2-40*(y-x^2)+2 -40*x;-40*x 20];
mm=@(x,y,z) ff(x,y) + gf(x,y)'*z + .5*z'*Hf(x,y)*z;
for j=1:N
    for k=1:N
        x=[5*(j-1)/(N-1)-1.001;5*(k-1)/(N-1)+0.001];
        X(j,k)=x(1);
        Y(j,k)=x(2);
    end
end

F=100*(Y-X.^2)^2+(1-X)^2;
contourf(X,Y,F,15)

x0 = [a,b];
hold on;
plot(x0(1),x0(2),'ok','MarkerFaceColor','k','MarkerSize',4)

delta_max = 4;
delta = 1;
eta = 0.20;

for k = 1:5 
    f = ff(x0(1), x0(2));
    g = gf(x0(1), x0(2));
    h = Hf(x0(1), x0(2));
    
    %pk = (-delta)*(g/norm(g));
    pu = g*((-g'*g)/(g'*h*g));
    x1 = x0;
    pt = [];
    
    if g'*h*g <= 0
        tau = 1;
    else
        tau = min((norm(g)^3)*inv((h*g)*(h*g)'));
    end
    %pc = tau * pk;
    while tau > 2
        tau = tau *.25;
    end
    if tau <=1
        pt = tau * pu;
    else
        if tau == 2
            pt = pu + (tau - 1)*(-h-pu);
        end
    end
    
    m = mm(x0(1), x0(2), pt);
    row = (f - ff(x0(1)+pt, x0(2)+pt))*inv(mm(x0(1),x0(2),0) - m);

    if row < 1/4
        delta = delta*.25;
    else
        if (row > .75) & (norm(pu) == delta)
            delta = min(2*delta,delta_max);
        end
    end
    
    if row > eta
        x1 = x0 + pt;
    end
    plot(x1(1),x1(2),'ok','MarkerFaceColor','k','MarkerSize',4)
    plot([x0(1);x1(1)],[x0(2);x1(2)],'k','LineWidth',2)
    x0 = x1
end

%{
ff=@(x,y) 100*(y-x^2)^2+(1-x)^2;
gf=@(x,y) [-400*x*(y-x^2)-2*(1-x);200*(y-x^2)]; 
Hf=@(x,y) [800*y^2-400*(y-x^2)+2 -400*x;-400*x 200];

x0=[a;b];
%hold on;
%plot(x0(1),x0(2),'ok','MarkerFaceColor','k','MarkerSize',4)
w=0.8;
c=0.1;
for j=1:17
    pk=-(inv(Hf(x0(1),x0(2)))*gf(x0(1),x0(2)));
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
    Hf(x0(1),x0(2))
    alpha
    x1=x0+alpha*pk
    %plot(x1(1),x1(2),'ok','MarkerFaceColor','k','MarkerSize',4)
    %plot([x0(1);x1(1)],[x0(2);x1(2)],'k','LineWidth',2)
    x0=x1;
    pause
end
%}
%---------------------------------------------------%