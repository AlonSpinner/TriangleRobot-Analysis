%% Q1: Inverse Kinematics for [x,y,theta]=[2,1,20degree];
x=2;
y=1;
theta=deg2rad(20);

fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
elbows=[1,1,1;
    -1,1,1;
    1,-1,1;
    1,1,-1;
    -1,-1,1;
    1,-1,-1;
    -1,1,-1;
    -1,-1,-1];
for i=1:8
    ax=subplot(2,4,i,'parent',fig);
    hold('on'); grid('on'); axis('equal'); axis('manual');
    xlim([0,5]); ylim([-1,4]);
    [d1,d2,d3]=Inverse_kinematics(x,y,theta,elbows(i,:));
    DrawMechanism(ax,x,y,theta,d1,d2,d3);
    title(sprintf('[s_1,s_2,s_3]=[%g,%g,%g]\n[d_1,d_2,d_3]=[%.2g,%.2g,%.2g]',...
        [elbows(i,:),[d1,d2,d3]]));
end
sgtitle(sprintf('x=2    y=1     theta=20^o\n'),'fontsize',20);
%% Q2: Finding Direct Kinematics
syms d1 d2 d3 x y theta L r H real
Eq1=(x-d1)^2+y^2==L^2;
Eq2=(x+r*cos(theta)-d2)^2+(y+r*sin(theta))^2==L^2;
Eq3=(x+r*cos(theta+pi/3)-d3)^2+(y+r*sin(theta+pi/3)-H)^2==L^2;

Eq12=collect(Eq1-Eq2,[x,y]);
Eq13=collect(Eq1-Eq3,[x,y]);

[A,b] = equationsToMatrix([Eq12;Eq13],[x;y]);

X=linsolve(A,b);
xt=expand(X(1)); %Expand breaks down cos(theta+pi/3) -> sin(theta)*cos(pi/3)+cos(theta)*sin(pi/3)
yt=expand(X(2));

syms T; %T=tan(theta/2)
sinT=2*T/(1+T^2);
cosT=(1-T^2)/(1+T^2);

xT=subs(xt,[sin(theta),cos(theta)],[sinT,cosT]);
yT=subs(yt,[sin(theta),cos(theta)],[sinT,cosT]);

Eq1T=subs(Eq1,[x,y],[xT,yT]);
load('Parameters.mat');
q=[1.7,2.5,2]';
Eq1Tnum=subs(Eq1T-L^2,[d1,d2,d3,r,L,H],[q',prm.r,prm.L,prm.H]); %deduct L^2 from both sides so left side compared to zero

ExprT=collect(lhs(Eq1Tnum),T);
[polyT,~]=numden(ExprT);
coffs=double(coeffs(polyT,'all'));
TVec=roots(coffs);

sVec = 2*TVec./(1+TVec.^2);
cVec = (1-TVec.^2)./(1+TVec.^2);
thetaVec = double(atan2(sVec,cVec));

xTNum=subs(xT,[d1,d2,d3,r,L,H],[q',prm.r,prm.L,prm.H]);
yTNum=subs(yT,[d1,d2,d3,r,L,H],[q',prm.r,prm.L,prm.H]);
xVec=double(subs(xTNum,T,TVec));
yVec=double(subs(yTNum,T,TVec));

%Implemented in Forward_Kinematics function
%% Q2 testing
q=[1.7,2.5,2]';
[xVec,yVec,thetaVec]=Forward_kinematics(q(1),q(2),q(3));

fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
dL=zeros(6,3);
for k=1:6
    ax=subplot(2,3,k,'parent',fig);
    hold('on'); grid('on'); axis('equal'); axis('manual');
    xlim([0,5]); ylim([-1,4]);
    dL(k,:)=DrawMechanism(ax,xVec(k),yVec(k),thetaVec(k),q(1),q(2),q(3));
    title(sprintf('x=%.3g   y=%.3g   theta=%.3g^o\ndL_1=%.2g   dL_2=%.2g   dL_3=%.2g',...
        [xVec(k),yVec(k),rad2deg(thetaVec(k)),dL(k,:)]))
end
Title=sgtitle(sprintf('d_1=%g     d_2=%g      d_3=%g',q'));
Title.FontSize=20;
%% Q3 
syms d1 d2 d3 x y theta L r H real
% Eq1=(x-d1)^2+y^2==L^2;
% Eq2=(x+r*cos(theta)-d2)^2+(y+r*sin(theta))^2==L^2;
% Eq3=(x-+r*cos(theta+pi/3)-d3)^2+(y+r*sin(theta+pi/3)-H)^2==L^2;

f1=(x-d1)^2+y^2-L^2;
f2=(x+r*cos(theta)-d2)^2+(y+r*sin(theta))^2-L^2;
f3=(x+r*cos(theta+pi/3)-d3)^2+(y+r*sin(theta+pi/3)-H)^2-L^2;

F=[f1;f2;f3];
Jx=jacobian(F,[x,y,theta]);
Jq=-jacobian(F,[d1,d2,d3]);

syms s1 s2 s3 real
[d1s,d2s,d3s]=Inverse_kinematics(x,y,theta,[s1 s2 s3],'method','symbolic');
Jx_xs=subs(Jx,[d1,d2,d3],[d1s,d2s,d3s]);
Jq_xs=subs(Jq,[d1,d2,d3],[d1s,d2s,d3s]);
%% Q4 Jx
load('Parameters.mat');
x0=2; theta0=0;
Jx_xs_num=subs(Jx_xs,[r,L,H,x,theta],[prm.r,prm.L,prm.H,x0,theta0]);
Det=det(Jx_xs_num);

elbows=[1,1,1;
    -1,1,1;
    1,-1,1;
    1,1,-1;
    -1,-1,1;
    1,-1,-1;
    -1,1,-1;
    -1,-1,-1];

yEq=sym(zeros(8,1));
S=cell(8,1);

assume(y >= 0); 
assume(y <= 2); 
for i=1:8
   yEq(i)=subs(Det,[s1,s2,s3], elbows(i,:));
   S{i}=solve(yEq(i));
end
assume(y,'clear');

ValidSolVec={};
for k=1:length(S)
    for j=1:length(S{k})
        y0=double(S{k}(j));
        s=elbows(k,:);
        Jxnum=double(subs(Jx_xs_num,[y,s1,s2,s3],[y0,s]));
        [d1,d2,d3]=Inverse_kinematics(x0,y0,theta0,s);
        %check conditions
        IsImagFlag=any(imag(Jxnum)>1e-7,'all');
        LinkGoodFlag=checkSolutionValidity(x0,y0,theta0,d1,d2,d3);
        DetZeroFlag=abs(det(Jxnum))<1e-7;   
%         disp([IsImagFlag,~LinkGoodFlag,DetZeroFlag]);
        if ~IsImagFlag
            Jxnum=real(Jxnum);
        end
        if ~any([IsImagFlag,~LinkGoodFlag,~DetZeroFlag])
            ValidSolVec{end+1}={y0,d1,d2,d3,Jxnum,s};
        end
    end
end
%% Q4 plotting
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
for i=1:4
    y0=ValidSolVec{i}{1};
    d1=ValidSolVec{i}{2};
    d2=ValidSolVec{i}{3};
    d3=ValidSolVec{i}{4};
    Jxnum=ValidSolVec{i}{5};
    s=ValidSolVec{i}{6};
    
    FreeVec=null(Jxnum);
    FreeVec(3)=rad2deg(FreeVec(3));
    
    ax=subplot(2,2,i,'parent',fig);
    hold('on'); grid('on'); axis('equal'); axis('manual');
    xlim([0,5]); ylim([-1,4]);
    DrawMechanism(ax,x0,y0,theta0,d1,d2,d3);
    quiver(x0,y0,FreeVec(1),FreeVec(2),0,'color',[0.8500 0.3250 0.0980],'LineWidth',0.8)
    
    Line1='x=%.3g   y=%.3g   theta=%.3g^o   ';
    Line2='d_1=%.2g   d_2=%.2g   d_3=%.2g   ';
    Line3='s_1=%.2g   s_2=%.2g   s_3=%.2g\n';
    Line4='Free direction d[x,y,theta]/dt=[%.2gm/s,%.2gm/s,%.2gdeg/s]';
    Title=title(sprintf([Line1,Line2,Line3,Line4],...
        [x0,y0,rad2deg(theta0),d1,d2,d3,s,FreeVec']));
end
Title=sgtitle(sprintf('Singular Poses for det(Jx)==0, theta=0, x=2, 0<y<2'));
Title.FontSize=15;

fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
for i=5:8
    y0=ValidSolVec{i}{1};
    d1=ValidSolVec{i}{2};
    d2=ValidSolVec{i}{3};
    d3=ValidSolVec{i}{4};
    Jxnum=ValidSolVec{i}{5};
    s=ValidSolVec{i}{6};
    
    FreeVec=null(Jxnum);
    FreeVec(3)=rad2deg(FreeVec(3));
    
    ax=subplot(2,2,i-4,'parent',fig);
    hold('on'); grid('on'); axis('equal'); axis('manual');
    xlim([0,5]); ylim([-1,4]);
    DrawMechanism(ax,x0,y0,theta0,d1,d2,d3);
    quiver(x0,y0,FreeVec(1),FreeVec(2),0,'color',[0.8500 0.3250 0.0980],'LineWidth',0.8)
    
    Line1='x=%.3g   y=%.3g   theta=%.3g^o   ';
    Line2='d_1=%.2g   d_2=%.2g   d_3=%.2g   ';
    Line3='s_1=%.2g   s_2=%.2g   s_3=%.2g\n';
    Line4='Free direction d[x,y,theta]/dt=[%.2gm/s,%.2gm/s,%.2gdeg/s]';
    Title=title(sprintf([Line1,Line2,Line3,Line4],...
        [x0,y0,rad2deg(theta0),d1,d2,d3,s,FreeVec']));
end
Title=sgtitle(sprintf('Singular Poses for det(Jx)==0, theta=0, x=2, 0<y<2'));
Title.FontSize=15;
%% Q4 Jq
load('Parameters.mat');
x0=2; theta0=0;
Jq_xs_num=subs(Jq_xs,[r,L,H,x,theta],[prm.r,prm.L,prm.H,x0,theta0]);
Det=det(Jq_xs_num);

elbows=[1,1,1;
    -1,1,1;
    1,-1,1;
    1,1,-1;
    -1,-1,1;
    1,-1,-1;
    -1,1,-1;
    -1,-1,-1];

yEq=sym(zeros(8,1));
S=cell(8,1);

assume(y >= 0); 
assume(y <= 2); 
for i=1:8
   yEq(i)=subs(Det,[s1,s2,s3], elbows(i,:));
   S{i}=solve(yEq(i));
end
assume(y,'clear');

ValidSolVec={};
for k=1:length(S)
    for j=1:length(S{k})
        y0=double(S{k}(j));
        s=elbows(k,:);
        Jqnum=double(subs(Jq_xs_num,[y,s1,s2,s3],[y0,s]));
        [d1,d2,d3]=Inverse_kinematics(x0,y0,theta0,s);
        %check conditions
        IsImagFlag=any(imag(Jxnum)>1e-7,'all');
        LinkGoodFlag=checkSolutionValidity(x0,y0,theta0,d1,d2,d3);
        DetZeroFlag=abs(det(Jxnum))<1e-7;   
%         disp([IsImagFlag,~LinkGoodFlag,DetZeroFlag]);
        if ~IsImagFlag
            Jxnum=real(Jxnum);
        end
        if ~any([IsImagFlag,~LinkGoodFlag,~DetZeroFlag])
            ValidSolVec{end+1}={y0,d1,d2,d3,Jxnum,s};
        end
    end
end
%% Q4 Jq plotting
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
for i=1:16
    y0=ValidSolVec{i}{1};
    d1=ValidSolVec{i}{2};
    d2=ValidSolVec{i}{3};
    d3=ValidSolVec{i}{4};
    Jqnum=ValidSolVec{i}{5};
    s=ValidSolVec{i}{6};
    
    FreeVec=null(Jqnum);
    FreeVec(3)=rad2deg(FreeVec(3));
    
    ax=subplot(4,4,i,'parent',fig);
    hold('on'); grid('on'); axis('equal'); axis('manual');
    xlim([0,5]); ylim([-1,4]);
    DrawMechanism(ax,x0,y0,theta0,d1,d2,d3);
    
    Line1='x=%.3g   y=%.3g   theta=%.3g^o\n';
    Line2='d_1=%.2g   d_2=%.2g   d_3=%.2g\n';
    Line3='s_1=%.2g   s_2=%.2g   s_3=%.2g\n';
    Line4='Free direction d[x,y,theta]/dt =\n  [%.2gm/s,%.2gm/s,%.2gdeg/s]';
%     title(sprintf([Line1,Line2,Line3,Line4],...
%         [x0,y0,rad2deg(theta0),d1,d2,d3,s,FreeVec']));
end

Title=sgtitle(sprintf('Singular Poses for det(Jx)==0, theta=0, x=2, 0<y<2'));
Title.FontSize=15;