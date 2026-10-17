%**************************************************************************
%       User documentation of any code that uses this code must cite the
%       Author in following article:
%
%       S. R. Nekoo, "Output- and State-Dependent Riccati Equation: An
%       Output Feedback Controller Design,"
%       Aerospace Science and Technology, 126, 107649, 2022.
%
%       Permission to modify this file (ONLY) and to distribute modified
%       code is granted, provided the above notice is retained, and a
%       notice that the code was modified is included with the above
%       statement.
%**************************************************************************
clc
clear
close all
% Time---------------------------------------------------------------------
tf=10;              % Final time of the simulation
N=999*tf;
Nt=tf/N;
t=0:tf/N:tf;
% Characteristics of the robot---------------------------------------------
a1=2;               % Length of the first link
a2=2;               % Length of the second link
ac1=a1/2;
ac2=a2/2;
m1=2;               % Mass of the first link
m2=2;               % Mass of the second link
mb=5;               % Mass of the base
mp=1;               % Mass of the payload
Izz1=1/12*m1*a1^2;
Izz2=1/12*m2*a2^2;
W1=1;
W2=0.5;
Izzb=1/12*mb*(4*W1^2+4*W2^2);
ex=W1-0.1;
ey=W2-0.5;
% Initial and final points-------------------------------------------------
xi=4.5;             % End-effector
yi=0.5;             % End-effector
xib=1.5;            % Center of mass of base
yib=-0.5;           % Center of mass
psii=0.1;           % Initial angle of the base
% Inverse Kinematic--------------------------------------------------------
xin=xi-(xib+ex*cos(psii)+ey*sin(psii));  % 末端位置-机械臂安装位置
yin=yi-(yib-ey*cos(psii)+ex*sin(psii));
qin2=real(acos((xin^2+yin^2-a1^2-a2^2)/(2*a1*a2)));
qin1=real(atan2(yin,xin))-real(atan2(a2*sin(qin2),(a1+a2*cos(qin2))));
qin1=qin1-psii;
% Initial condition--------------------------------------------------------
qin=[qin1;qin2];
dqin=[-0.15;0.2];           % velocity of the first and second link
xbin=[xib;yib;psii];
dxbin=[-0.35;0.1;-0.1];     % velocity of the base, linear, dxc, dyc, and angular velocity d_psi
xin=[xbin;qin;dxbin;dqin];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution SDRE------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(:,1)=xin;
for i=1:1:length(t)
    Time=t(i)
    % States---------------------------------------------------------------
    xc=x(1,i);
    yc=x(2,i);
    psi=x(3,i);
    q1=x(4,i);
    q2=x(5,i);
    dxc=x(6,i);
    dyc=x(7,i);
    dpsi=x(8,i);
    dq1=x(9,i);
    dq2=x(10,i);

    q=[q1;q2];
    dq=[dq1;dq2];
    dqb=[dxc;dyc;dpsi];
    % Model----------------------------------------------------------------
    m11=m1+m2+mb+mp;
    m12=0;
    m13=ey*m1*cos(psi)+ey*m2*cos(psi)+ey*mp*cos(psi)-ex*m1*sin(psi)-ex*m2*sin(psi)-ex*mp*sin(psi)-ac2*m2*sin(psi+q1+q2)-a2*mp*sin(psi+q1+q2)-a1*m2*sin(psi+q1)-ac1*m1*sin(psi+q1)-a1*mp*sin(psi+q1);
    m14=-ac2*m2*sin(psi+q1+q2)-a2*mp*sin(psi+q1+q2)-a1*m2*sin(psi+q1)-ac1*m1*sin(psi+q1)-a1*mp*sin(psi+q1);
    m15=-sin(psi+q1+q2)*(ac2*m2+a2*mp);
    m22=m1+m2+mb+mp;
    m23=ex*m1*cos(psi)+ex*m2*cos(psi)+ex*mp*cos(psi)+ey*m1*sin(psi)+ey*m2*sin(psi)+ey*mp*sin(psi)+ac2*m2*cos(psi+q1+q2)+a2*mp*cos(psi+q1+q2)+a1*m2*cos(psi+q1)+ac1*m1*cos(psi+q1)+a1*mp*cos(psi+q1);
    m24=ac2*m2*cos(psi+q1+q2)+a2*mp*cos(psi+q1+q2)+a1*m2*cos(psi+q1)+ac1*m1*cos(psi+q1)+a1*mp*cos(psi+q1);
    m25=cos(psi+q1+q2)*(ac2*m2+a2*mp);
    m33=Izz1+Izz2+Izzb+a1^2*m2+ac1^2*m1+ac2^2*m2+a1^2*mp+a2^2*mp+ex^2*m1+ex^2*m2+ey^2*m1+ey^2*m2+ex^2*mp+ey^2*mp+2*ac2*ex*m2*cos(q1+q2)+2*a2*ex*mp*cos(q1+q2)-2*ac2*ey*m2*sin(q1+q2)-2*a2*ey*mp*sin(q1+q2)+2*a1*ac2*m2*cos(q2)+2*a1*a2*mp*cos(q2)+2*a1*ex*m2*cos(q1)+2*ac1*ex*m1*cos(q1)+2*a1*ex*mp*cos(q1)-2*a1*ey*m2*sin(q1)-2*ac1*ey*m1*sin(q1)-2*a1*ey*mp*sin(q1);
    m34=Izz1+Izz2+a1^2*m2+ac1^2*m1+ac2^2*m2+a1^2*mp+a2^2*mp+ac2*ex*m2*cos(q1+q2)+a2*ex*mp*cos(q1+q2)-ac2*ey*m2*sin(q1+q2)-a2*ey*mp*sin(q1+q2)+2*a1*ac2*m2*cos(q2)+2*a1*a2*mp*cos(q2)+a1*ex*m2*cos(q1)+ac1*ex*m1*cos(q1)+a1*ex*mp*cos(q1)-a1*ey*m2*sin(q1)-ac1*ey*m1*sin(q1)-a1*ey*mp*sin(q1);
    m35=Izz2+ac2^2*m2+a2^2*mp+ac2*ex*m2*cos(q1+q2)+a2*ex*mp*cos(q1+q2)-ac2*ey*m2*sin(q1+q2)-a2*ey*mp*sin(q1+q2)+a1*ac2*m2*cos(q2)+a1*a2*mp*cos(q2);
    m44=Izz1+Izz2+a1^2*m2+ac1^2*m1+ac2^2*m2+a1^2*mp+a2^2*mp+2*a1*ac2*m2*cos(q2)+2*a1*a2*mp*cos(q2);
    m45=mp*a2^2+a1*mp*cos(q2)*a2+m2*ac2^2+a1*m2*cos(q2)*ac2+Izz2;
    m55=mp*a2^2+m2*ac2^2+Izz2;
    % Matrices-------------------------------------------------------------
    Mb=[m11,m12,m13;...
        m12,m22,m23;...
        m13,m23,m33];

    Mba=[m14,m15;...
        m24,m25;...
        m34,m35];

    Ma=[m44,m45;...
        m45,m55];

    M=[Mb,Mba;...
        Mba',Ma];

    Me=Ma-Mba'*inv(Mb)*Mba;

    c1=-ac2*dpsi^2*m2*cos(psi+q1+q2)-ac2*dq1^2*m2*cos(psi+q1+q2)-ac2*dq2^2*m2*cos(psi+q1+q2)-a2*dpsi^2*mp*cos(psi+q1+q2)-a2*dq1^2*mp*cos(psi+q1+q2)-a2*dq2^2*mp*cos(psi+q1+q2)-a1*dpsi^2*m2*cos(psi+q1)-a1*dq1^2*m2*cos(psi+q1)-ac1*dpsi^2*m1*cos(psi+q1)-ac1*dq1^2*m1*cos(psi+q1)-a1*dpsi^2*mp*cos(psi+q1)-a1*dq1^2*mp*cos(psi+q1)-dpsi^2*ex*m1*cos(psi)-dpsi^2*ex*m2*cos(psi)-dpsi^2*ex*mp*cos(psi)-dpsi^2*ey*m1*sin(psi)-dpsi^2*ey*m2*sin(psi)-dpsi^2*ey*mp*sin(psi)-2*ac2*dpsi*dq1*m2*cos(psi+q1+q2)-2*ac2*dpsi*dq2*m2*cos(psi+q1+q2)-2*ac2*dq1*dq2*m2*cos(psi+q1+q2)-2*a2*dpsi*dq1*mp*cos(psi+q1+q2)-2*a2*dpsi*dq2*mp*cos(psi+q1+q2)-2*a2*dq1*dq2*mp*cos(psi+q1+q2)-2*a1*dpsi*dq1*m2*cos(psi+q1)-2*ac1*dpsi*dq1*m1*cos(psi+q1)-2*a1*dpsi*dq1*mp*cos(psi+q1);
    c2=dpsi^2*ey*m1*cos(psi)-ac2*dq1^2*m2*sin(psi+q1+q2)-ac2*dq2^2*m2*sin(psi+q1+q2)-a2*dpsi^2*mp*sin(psi+q1+q2)-a2*dq1^2*mp*sin(psi+q1+q2)-a2*dq2^2*mp*sin(psi+q1+q2)-a1*dpsi^2*m2*sin(psi+q1)-a1*dq1^2*m2*sin(psi+q1)-ac1*dpsi^2*m1*sin(psi+q1)-ac1*dq1^2*m1*sin(psi+q1)-a1*dpsi^2*mp*sin(psi+q1)-a1*dq1^2*mp*sin(psi+q1)-ac2*dpsi^2*m2*sin(psi+q1+q2)+dpsi^2*ey*m2*cos(psi)+dpsi^2*ey*mp*cos(psi)-dpsi^2*ex*m1*sin(psi)-dpsi^2*ex*m2*sin(psi)-dpsi^2*ex*mp*sin(psi)-2*ac2*dpsi*dq1*m2*sin(psi+q1+q2)-2*ac2*dpsi*dq2*m2*sin(psi+q1+q2)-2*ac2*dq1*dq2*m2*sin(psi+q1+q2)-2*a2*dpsi*dq1*mp*sin(psi+q1+q2)-2*a2*dpsi*dq2*mp*sin(psi+q1+q2)-2*a2*dq1*dq2*mp*sin(psi+q1+q2)-2*a1*dpsi*dq1*m2*sin(psi+q1)-2*ac1*dpsi*dq1*m1*sin(psi+q1)-2*a1*dpsi*dq1*mp*sin(psi+q1);
    c3=-ac2*dq1^2*ex*m2*sin(q1+q2)-ac2*dq2^2*ex*m2*sin(q1+q2)-a2*dq1^2*ex*mp*sin(q1+q2)-a2*dq2^2*ex*mp*sin(q1+q2)-a1*dq1^2*ey*m2*cos(q1)-ac1*dq1^2*ey*m1*cos(q1)-a1*dq1^2*ey*mp*cos(q1)-a1*ac2*dq2^2*m2*sin(q2)-a1*a2*dq2^2*mp*sin(q2)-a1*dq1^2*ex*m2*sin(q1)-ac1*dq1^2*ex*m1*sin(q1)-a1*dq1^2*ex*mp*sin(q1)-ac2*dq1^2*ey*m2*cos(q1+q2)-ac2*dq2^2*ey*m2*cos(q1+q2)-a2*dq1^2*ey*mp*cos(q1+q2)-a2*dq2^2*ey*mp*cos(q1+q2)-2*ac2*dpsi*dq1*ey*m2*cos(q1+q2)-2*ac2*dpsi*dq2*ey*m2*cos(q1+q2)-2*ac2*dq1*dq2*ey*m2*cos(q1+q2)-2*a2*dpsi*dq1*ey*mp*cos(q1+q2)-2*a2*dpsi*dq2*ey*mp*cos(q1+q2)-2*a2*dq1*dq2*ey*mp*cos(q1+q2)-2*ac2*dpsi*dq1*ex*m2*sin(q1+q2)-2*ac2*dpsi*dq2*ex*m2*sin(q1+q2)-2*ac2*dq1*dq2*ex*m2*sin(q1+q2)-2*a2*dpsi*dq1*ex*mp*sin(q1+q2)-2*a2*dpsi*dq2*ex*mp*sin(q1+q2)-2*a2*dq1*dq2*ex*mp*sin(q1+q2)-2*a1*dpsi*dq1*ey*m2*cos(q1)-2*ac1*dpsi*dq1*ey*m1*cos(q1)-2*a1*dpsi*dq1*ey*mp*cos(q1)-2*a1*ac2*dpsi*dq2*m2*sin(q2)-2*a1*ac2*dq1*dq2*m2*sin(q2)-2*a1*a2*dpsi*dq2*mp*sin(q2)-2*a1*a2*dq1*dq2*mp*sin(q2)-2*a1*dpsi*dq1*ex*m2*sin(q1)-2*ac1*dpsi*dq1*ex*m1*sin(q1)-2*a1*dpsi*dq1*ex*mp*sin(q1);
    c4=ac2*dpsi^2*ex*m2*sin(q1+q2)+a2*dpsi^2*ex*mp*sin(q1+q2)+a1*dpsi^2*ey*m2*cos(q1)+ac1*dpsi^2*ey*m1*cos(q1)+a1*dpsi^2*ey*mp*cos(q1)-a1*ac2*dq2^2*m2*sin(q2)-a1*a2*dq2^2*mp*sin(q2)+a1*dpsi^2*ex*m2*sin(q1)+ac1*dpsi^2*ex*m1*sin(q1)+a1*dpsi^2*ex*mp*sin(q1)+ac2*dpsi^2*ey*m2*cos(q1+q2)+a2*dpsi^2*ey*mp*cos(q1+q2)-2*a1*ac2*dpsi*dq2*m2*sin(q2)-2*a1*ac2*dq1*dq2*m2*sin(q2)-2*a1*a2*dpsi*dq2*mp*sin(q2)-2*a1*a2*dq1*dq2*mp*sin(q2);
    c5=(ac2*m2+a2*mp)*(dpsi^2*ey*cos(q1+q2)+dpsi^2*ex*sin(q1+q2)+a1*dpsi^2*sin(q2)+a1*dq1^2*sin(q2)+2*a1*dpsi*dq1*sin(q2));

    cb=[c1;c2;c3];
    ca=[c4;c5];
    ce=ca-Mba'*inv(Mb)*cb;
    c=[cb;ca];
    % Input----------------------------------------------------------------
    u(:,i)=[0;0];
    % System---------------------------------------------------------------
    x(:,i+1)=x(:,i)+[dqb;dq;inv(M)*([zeros(3,1);u(:,i)]-c)]*Nt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots--------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold all;grid on;box on
yyaxis left
plot(t,x(1,1:length(t)),'-k','Linewidth',1.5)
plot(t,x(2,1:length(t)),'--b','Linewidth',1.5)
xlabel('t (s)')
ylabel('base position (m)')
yyaxis right
plot(t,x(3,1:length(t)),'-.r','Linewidth',1.5)
ylabel('base rotation (rad)')
legend('x_c','y_c','\psi')

figure(2)
hold all;grid on;box on
yyaxis left
plot(t,x(4,1:length(t)),'-k','Linewidth',1.5)
plot(t,x(5,1:length(t)),'--b','Linewidth',1.5)
xlabel('t (s)')
ylabel('angle of link (rad)')
xlim([0,tf])
yyaxis right
plot(t,x(6,1:length(t)),'-.r','Linewidth',1.5)
plot(t,x(7,1:length(t)),':','Linewidth',1.5,'Color',[0.47,0.67,0.18])
ylabel('angular velocity of link (rad/s)')
legend('1','2','1','2')

xc=x(1,:);
yc=x(2,:);
psi=x(3,:);
q1=x(4,:);
q2=x(5,:);
dq1=x(6,:);
dq2=x(7,:);
xb=xc+ex*cos(psi)+ey*sin(psi);
yb=yc-ey*cos(psi)+ex*sin(psi);
xe1=xc+a1*cos(psi+q1)+ex*cos(psi)+ey*sin(psi);
ye1=yc+a1*sin(psi+q1)-ey*cos(psi)+ex*sin(psi);
xe2=xc+a1*cos(psi+q1)+ex*cos(psi)+ey*sin(psi)+a2*cos(psi+q1+q2);
ye2=yc+a1*sin(psi+q1)-ey*cos(psi)+ex*sin(psi)+a2*sin(psi+q1+q2);

figure(6)
hold all;grid on;box on
plot(xi,yi,'*r')
plot(xe2,ye2,'-k','Linewidth',1.5)
plot(xb,yb,'-r','Linewidth',1.5)
for i=1:round(length(t)/25):length(t)
    x11=xc(i)+W1.*cos(psi(i))-W2.*sin(psi(i));
    y11=yc(i)+W1.*sin(psi(i))+W2.*sin(psi(i)+pi/2);
    x21=xc(i)-W1.*cos(psi(i))+W2.*cos(psi(i)+pi/2);
    y21=yc(i)-W1.*sin(psi(i))+W2.*sin(psi(i)+pi/2);
    x31=xc(i)+W1.*cos(psi(i))+W2.*sin(psi(i));
    y31=yc(i)+W1.*sin(psi(i))-W2.*sin(psi(i)+pi/2);
    x41=xc(i)-W1.*cos(psi(i))-W2.*cos(psi(i)+pi/2);
    y41=yc(i)-W1.*sin(psi(i))-W2.*sin(psi(i)+pi/2);
    line([xb(i);xe1(i)],[yb(i);ye1(i)],'color',[0.3,0.75,1])
    line([xe1(i);xe2(i)],[ye1(i);ye2(i)],'color',[0.3,0.75,1])
    line([x11;x31],[y11;y31],'color',[0.9,0.7,0.15])
    line([x41;x31],[y41;y31],'color',[0.9,0.7,0.15])
    line([x21;x41],[y21;y41],'color',[0.9,0.7,0.15])
    line([x21;x11],[y21;y11],'color',[0.9,0.7,0.15])
end
plot(xe2,ye2,'-k','Linewidth',1.5)
plot(xb,yb,'-r','Linewidth',1.5)
xlabel('X(m)')
ylabel('Y(m)')
legend('start point','end-effector trajectory','base trajectory','links')
axis equal

figure(5000)
hp=round(length(t)/250);
newVid2=VideoWriter('Openloop_FFSM_2Links');
newVid2.FrameRate=10;
newVid2.Quality=100;
open(newVid2);
hold all; grid on; box on
plot(xi,yi,'*r')
plot(xe2,ye2,':k','Linewidth',1)
plot(xb,yb,':r','Linewidth',1)
xlabel('X(m)')
ylabel('Y(m)')
axis equal
for i=hp+1:hp:length(t)
    x11=xc(i)+W1.*cos(psi(i))-W2.*sin(psi(i));
    y11=yc(i)+W1.*sin(psi(i))+W2.*sin(psi(i)+pi/2);
    x21=xc(i)-W1.*cos(psi(i))+W2.*cos(psi(i)+pi/2);
    y21=yc(i)-W1.*sin(psi(i))+W2.*sin(psi(i)+pi/2);
    x31=xc(i)+W1.*cos(psi(i))+W2.*sin(psi(i));
    y31=yc(i)+W1.*sin(psi(i))-W2.*sin(psi(i)+pi/2);
    x41=xc(i)-W1.*cos(psi(i))-W2.*cos(psi(i)+pi/2);
    y41=yc(i)-W1.*sin(psi(i))-W2.*sin(psi(i)+pi/2);

    Line1(i)=line([xb(i);xe1(i)],[yb(i);ye1(i)],'color',[0.3,0.75,1],'Linewidth',1.5);
    Line2(i)=line([xe1(i);xe2(i)],[ye1(i);ye2(i)],'color',[0.3,0.75,1],'Linewidth',1.5);
    Line3(i)=line([x11;x31],[y11;y31],'color',[0.9,0.7,0.15],'Linewidth',1.5);
    Line4(i)=line([x41;x31],[y41;y31],'color',[0.9,0.7,0.15],'Linewidth',1.5);
    Line5(i)=line([x21;x41],[y21;y41],'color',[0.9,0.7,0.15],'Linewidth',1.5);
    Line6(i)=line([x21;x11],[y21;y11],'color',[0.9,0.7,0.15],'Linewidth',1.5);

    Plot1(i)=plot(xe2(i),ye2(i),'squarek');

    if i>hp+1
        delete(Line1(i-hp));
        delete(Line2(i-hp));
        delete(Line3(i-hp));
        delete(Line4(i-hp));
        delete(Line5(i-hp));
        delete(Line6(i-hp));
        delete(Plot1(i-hp));
    end

    ZZ=getframe(figure(5000));
    writeVideo(newVid2,ZZ)
end
close(newVid2)