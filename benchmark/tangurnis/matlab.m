function matlab

format LONGG

global a;
global Ra;
global k;
global beta_;

Ra =-1;

figure(1);clf;

for method={'tala','tala_c', 'ba'}

convx=[];
convy=[];
first=1;
for j=[8 16 32 64 128]
    
    vv=load(sprintf('%s/vel_%d.csv',method{1},j));
    if (first==1)
        first = 0;
        
        Di = vv(1,1); % 0.5
        gamm = vv(1,2); % 1
        beta_ = Di/gamm;

        k = vv(1,3)*pi; %1*pi
        a=vv(1,4); % 0 or 2
        fprintf('computing Di=%i, gamma = %i\n', Di, gamm)

        tic
        %variables: u_z u_x eta_zz eta_xz
        %initial values: 0 ? ? 0
        solinit = bvpinit(linspace(0,1,5000),[0 0 0 0],[]);
        SOL = bvp4c(@myode,@bcfun,solinit, bvpset('RelTol', 1e-12,'AbsTol',1e-12));
        toc
    end
    vv=vv(2:size(vv,1),:); %remove the first line
    
    %now vv is an array with x,y,u_x,u_y,jxw,p,T
    
    %we are looking at u_x:
    dat_calc = vv(:,3);
    % interpolate reference solution on points of the solution (in 1d)
    dat_ref = interp1(SOL.x, SOL.y(2,:), vv(:,2)); 
    dat_ref = dat_ref .* sin(vv(:,1).*k); % and make it 2d
    
    fprintf('%i,%i %i,%i\n',min(dat_calc),max(dat_calc),min(dat_ref),max(dat_ref))
    
    %scatter3(vv(:,1),vv(:,2),dat_ref,[],dat_ref)
    %scatter3(vv(:,1),vv(:,2),vv(:,3),[],vv(:,3))
    
    %compute \sqrt(\int (u-u_ref)^2)
    l2error = sqrt( sum((dat_ref - dat_calc).^2 .* vv(:,5)) );
    
    convx=[convx; j];
    convy=[convy; l2error];
    
end
figure(1);
loglog(convx,convy); hold on;
convorder=polyfit(log(convx),log(convy),1);
convorder=convorder(1);

display([convx convy]);
fprintf('convergence order: %f\n',-convorder);
end


xlabel('1/h');
ylabel('L2 error of u_x');
xlim([4 256]);
set(gca,'XTick',[8 16 32 64 128]);

ylim([5e-10 5e-5]);
legend('tala a=0','tala a=2', 'ba')



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