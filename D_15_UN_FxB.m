% Analyze 15_Story W/O FIXED BASE

clear all
clc
%% EarthQuake
% E01:E22= Comp. 1 Far-Field, E61:E74= Comp. 1 Near-Field
NEQ=0;
for EQ=100%[1:22,61:74]
NEQ=NEQ+1;
if EQ==100 %ElCentro
dt=0.02;
tmax=20;
else
A = importdata('E:\\Codes\\EQ\\EQdata.txt');
dt=A(EQ,2);
tmax=A(EQ,3);
end
fopenT=sprintf('E:\\Codes\\EQ\\E%d.txt',EQ);
fileID=fopen(fopenT,'r');
formatSpec='%f';
G=fscanf(fileID,formatSpec);

d=15;

% First Story
ms(1)=450000;
cc(1)=26170;
st(1)=18050000;
% Other Story
 for i=2:15
     ms(i)=345600;
     cc(i)=293700;
     st(i)=340400000;
 end

for i=1:d
    m(i,i)=ms(i);
    if i<d
        c(i,i)=cc(i)+cc(i+1);
        c(i+1,i)=-cc(i+1);
        c(i,i+1)=-cc(i+1);
        k(i,i)=st(i)+st(i+1);
        k(i+1,i)=-st(i+1);
        k(i,i+1)=-st(i+1);
    else
        c(i,i)=cc(i);
        k(i,i)=st(i);
    end
end

t=0;
n=1;
e=(-1)*ones(d,1);
xn=zeros(d,1);
vn=zeros(d,1); 
an=e*G(1)*9.80665;     
Ga=0.5;
Be=1/4;
a1=((1/(Be*(dt^2)))*m)+((Ga/(Be*dt))*c);
a2=((1/(Be*dt))*m)+(((Ga/Be)-1)*c);
a3=(((1/(2*Be))-1)*m)+(dt*((Ga/(2*Be))-1)*c);
Tk=k+a1;
while t<tmax   
    pm=m*e*9.80665*G(n+1);
    TDpn=pm+a1*xn+a2*vn+a3*an;
    xm=Tk\TDpn;
    vm=((Ga/(Be*dt))*(xm-xn))+((1-(Ga/Be))*vn)+(dt*(1-(Ga/(2*Be)))*an);
    am=((1/(Be*(dt^2)))*(xm-xn))-((1/(Be*dt))*vn)-(((1/(2*Be))-1)*an);
    maxGdis(n)=max(abs(xn(1:d,1)));
    Mt(n,:)=t;
    Mx(n,:)=xn;
    Mv(n,:)=vn;
    Ma(n,:)=an;
    xn=xm;
    vn=vm;
    an=am;
    n=n+1;
    t=t+dt;
end

for i=1:15
    if i==1
        MaxDrift(1,i)=max(abs(Mx(:,1)));
    else
        MaxDrift(1,i)=max(abs(Mx(:,i)-Mx(:,i-1)));
    end
end
A_DRIFT(NEQ,1)=max(MaxDrift);

maxSdis=max(maxGdis);
RMS=rms(Mx);
maxGdis=0;
disp(['EQ',num2str(EQ),' : ',num2str(maxSdis)])
A_RESULT(NEQ,1)=maxSdis;
A_RMS(NEQ,1)=max(RMS);
disp('------')
clear Mt
clear Mx
clear Mv
clear Ma
end

%  [V,D]=eig(k,m);
%  [W,kk]=sort(diag(D));
%  V=V(:,kk);
%  Phi=V*inv(sqrt(diag(diag(V'*m*V))));
%  Omega=diag(sqrt(Phi'*k*Phi));
%  Omega(1:3)