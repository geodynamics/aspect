function matlab

%etastar = exp(a*depth);
%omega = exp(beta*depth)*Pi*depth; %rho*t;
global a;
global Ra;
global k;
global beta_;

Di = 0.5;
gamm = 1;

beta_ = Di/gamm;
k = 1*pi;
Ra =-1;

a=0;

%variables: u_z u_x eta_zz eta_xz
%initial values: 0 ? ? 0

tic
solinit = bvpinit(linspace(0,1,5000),[0 0 0 0],[]);
SOL = bvp4c(@myode,@bcfun,solinit, bvpset('RelTol', 1e-12,'AbsTol',1e-12));

toc

convx=[];
convy=[];
for j=[8 16 32 64 128]
    
    vv=load(sprintf('vel_%d.dat',j));
    vv=vv(2:size(vv,1),:);
    
    dat_calc = vv(:,3);
    dat_ref = interp1(SOL.x, SOL.y(2,:), vv(:,2)); % interpolate in 1d
    dat_ref = dat_ref .* sin(vv(:,1).*k); % and make it 2d
    
    l2error = sqrt( sum((dat_ref - dat_calc).^2 .* vv(:,5)) );
    
    convx=[convx; j];
    convy=[convy; l2error];
    
end

figure(1);clf;
loglog(convx,convy,'o-');
xlabel('1/h');
ylabel('L2 error of u_x');
xlim([4 256]);
set(gca,'XTick',[8 16 32 64 128]);

ylim([5e-10 5e-5]);

convorder=polyfit(log(convx),log(convy),1);
convorder=convorder(1);

[convx convy]

convorder

%  figure(2);clf;
%  plot(dat_calc);hold on;
%  plot(dat_ref,'r');

%  figure(1);clf;
%  plot(x,SOL.y(1,:),'r');
%  hold on;
%  plot(v(:,2),v(:,4));
%  %f=load('fort.100');
%  %plot(f(:,1),f(:,2),'k');
%
%  figure(3);clf;
%  plot(x,SOL.y(2,:),'r');
%  hold on;
%  %v(:,3)=v(:,3)./sin(k*v(:,1));
%  plot(v(:,2),v(:,3));
%  f=load('fort.101');
%  plot(f(:,1),f(:,2),'k');
%
%  U=linspace(0,1,33);
%  aa=round(size(SOL.x,2)/33);
%  sub = 1:aa:size(SOL.x,2);
%  V=SOL.x(sub);
%  DU=0.5*SOL.y(2,sub)'*sin(k*U);
%  DV=0.5*SOL.y(1,sub)'*cos(k*U);
%  figure(4);clf;
%  quiver(U,V,DU,DV);
%  hold on;
%  quiver(v(:,1),v(:,2),0.5*v(:,3),0.5*v(:,4),'r');
%  legend('matlab','deal');

end

function res = bcfun(ya,yb,p)
res = [ya(1);yb(1);ya(4);yb(4)];
end

function [r]=myode(x,Y,p)

global a;
global Ra;
global k;
global beta_;

etazero = 1;
etastar = exp(a*(1-x))/etazero;

omega=exp(beta_*(1-x))*sin(pi*(1-x));

M=[beta_ -k 0 0; k 0 0 2*k/etastar; 0 0 0 -k;-beta_*etastar 2*k*etastar k 0];
r = M * Y + [0;0;omega*Ra/(2*etazero*k);0];

end