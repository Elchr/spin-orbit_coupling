clc; clear;
k=1;
d=1;
n=400;
r=1;
T=2.4325;;
dt=T/400;
total=n*T;
nstep=fix(total/dt);

b=1;
heta=1/sqrt(0.001);
numm=numel(-20:1:20);
zold=zeros(nstep,numm);
pzold=zeros(nstep,numm);
znew=zeros(nstep,numm);
pznew=zeros(nstep,numm);
t=zeros(nstep,numm);

gold=zeros(nstep,numm);
pgold=zeros(nstep,numm);
gnew=zeros(nstep,numm);
pgnew=zeros(nstep,numm);
pth=zeros(nstep,numm);

for h=-20:1:20;
    zold(1,k)=pi/3+0.1879;
    pzold(1,k)=0;
    gold(1,k)=0;
    pgold(1,k)=2*b*h/20;
for i=1:1:nstep;
        
  t(i,k)=i*dt;
     
  pznew(i,k)=+pzold(i,k)-3*(1-(2-2*cos(zold(i,k)))^(-3/2))*sin(zold(i,k))*dt;
  znew(i,k)=zold(i,k)+pznew(i,k)*dt;
  
  pgnew(i,k)=+pgold(i,k)+(-b^2/2*sin(2*(gold(i,k)-d*zold(i,k)+d*zold(1,k))))*dt;
  gnew(i,k)=mod((gold(i,k)+pgnew(i,k)*dt+pi/2),pi)-pi/2;
  pth(i,k)=pgnew(i,k)/heta+1;  
    
  zold(i+1,k)=znew(i,k);
  pzold(i+1,k)=pznew(i,k);
  gold(i+1,k)=gnew(i,k);
  pgold(i+1,k)=pgnew(i,k);
  
 
end
k=k+1;
end 


alpha=abs(max(znew(:,1))-min(znew(:,1)))*(180/(2*pi));
  
z_f=zeros(n,numm);
pz_f=zeros(n,numm);
g_f=zeros(n,numm);
pth_f=zeros(n,numm);
pg_f=zeros(n,numm);
for jj=1:1:numm
for j=1:1:n
z_f(j,jj)=zold(j*(T/dt),jj);
pz_f(j,jj)=pzold(j*(T/dt),jj); 
g_f(j,jj)=gold(j*(T/dt),jj);
pth_f(j,jj)=pth(j*(T/dt),jj);
pg_f(j,jj)=pgold(j*(T/dt),jj);
end
end


figure(1)
plot (z_f,pz_f,'.k')
xlabel('\zeta')
ylabel('$\dot{\zeta}$','interpreter','latex')
grid on

%legend('a=0.7','b=0.04')

figure(2)

colorspec = 'kbgrymkbgrymkbgrymkbgrymkbgrymkbgrymkbgry';
cla;
 for i=1:1:numm 
scatter(g_f(:,i),pth_f(:,i),'.',colorspec(i))
hold on
 end
%plot (g_f,pth_f,'.k')
xlabel('\gamma')
ylabel('$\dot{\theta}/ \eta $','interpreter','latex')
grid on
title('Phase Diagram')

