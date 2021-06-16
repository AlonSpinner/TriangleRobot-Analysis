function dL=DrawMechanism(ax,x,y,theta,d1,d2,d3)
load('Parameters.mat');
r=prm.r;
H=prm.H;
L=prm.L;

if isempty(ax)
    fig=figure('color',[1,1,1]);
    ax=axes(fig);
    grid(ax,'on'); hold(ax,'on'); axis(ax,'equal');
    ax.XLim=[0,5];
    ax.YLim=[-1,4];
end

D1=[d1,0];
D2=[d2,0];
D3=[d3,H];
A=[x,y];
B=A+r*[cos(theta),sin(theta)];
C=A+r*[cos(theta+pi/3),sin(theta+pi/3)];

xy_tri=[A;B;C];
htri=patch(ax,xy_tri(:,1),xy_tri(:,2),'red');

hbox=r/3; wbox=r/2;
hD1=drawbox(ax,D1,wbox,hbox);
hD2=drawbox(ax,D2,wbox,hbox);
hD3=drawbox(ax,D3,wbox,hbox);

hL1=plot([D1(1),A(1)],[D1(2),A(2)],'linewidth',2,'color',[0,0.5,0]);
hL2=plot([D2(1),B(1)],[D2(2),B(2)],'linewidth',2,'color',[0,0.5,0]);
hL3=plot([D3(1),C(1)],[D3(2),C(2)],'linewidth',2,'color',[0,0.5,0]);

dL1=abs(vecnorm(A-D1)-L);
dL2=abs(vecnorm(B-D2)-L);
dL3=abs(vecnorm(C-D3)-L);
dL=[dL1,dL2,dL3];

plot(ax,ax.XLim,[0,0],'linewidth',2,'color','k');
plot(ax,ax.XLim,[H,H],'linewidth',2,'color','k');
end

function hbox=drawbox(ax,P,w,h)
x=[-1,1,1,-1]*(w/2)+P(1);
y=[-1,-1,1,1]*(h/2)+P(2);
hbox=patch(ax,x,y,'blue');
end