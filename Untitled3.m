clc; clear;

n=400;
r=1;
T=0.0097;
dt=T/1000;
total=n*T;
nstep=fix(total/dt);

b=1;%10^(-1.5);
k=1;
numm=numel(-10:1:10);
zold=zeros(nstep,numm);
pzold=zeros(nstep,numm);
znew=zeros(nstep,numm);
pznew=zeros(nstep,numm);
t=zeros(nstep,numm);
per=zeros(nstep,1);
for h=-10:1:10;
    zold(1,k)=pi/3;
    pzold(1,k)=0;
    
for i=1:1:nstep;
        
  t(i,k)=i*dt;
     
  pznew(i,k)=+pzold(i,k)-3*(1-(2-2*cos(zold(i,k)))^(-3/2))*sin(zold(i,k))*dt;
  znew(i,k)=mod(zold(i,k)+pznew(i,k)*dt,2*pi);
   
  if pzold(i,k)>0 && pznew(i,k)<0
  per=t(i,k);
      
  end 
  
  zold(i+1,k)=znew(i,k);
  pzold(i+1,k)=pznew(i,k);
   
end
k=k+1;
end 

alpha=abs(max(znew(:,1))-min(znew(:,1)))*(180/(2*pi))
  
z_f=zeros(n,numm);
pz_f=zeros(n,numm);

for jj=1:1:numm
for j=1:1:n
z_f(j,jj)=zold(j*(T/dt),jj);
pz_f(j,jj)=pzold(j*(T/dt),jj); 

end
end


figure(1)
plot (z_f,pz_f,'.k')
xlabel('\zeta')
ylabel('$\dot{\zeta}$','interpreter','latex')
grid on


