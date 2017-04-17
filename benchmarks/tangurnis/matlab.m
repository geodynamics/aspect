% tala   Di=0.5, gamma = 1, a= 0 heating terms: adab = 5.251998887401354e-03, shear = -5.306915879791999e-03
% tala_c Di=0.5, gamma = 1, a= 2 heating terms: adab = 1.718840927336709e-03, shear = -1.671297523892296e-03
% BA     Di=0.0, gamma = 1, a= 0 heating terms: adab = 0, shear = 0 ;-)

function matlab

format LONGG

global a;
global Ra;
global k;
global beta_;

Ra =-1;

figure(1);clf;

for method={'tala','tala_c', 'ba'}

fprintf('*** method %s ***\n', method{1})
    
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
        fprintf('computing Di=%i, gamma = %i, a= %i\n', Di, gamm, a)

        tic
        %variables: u_z u_x eta_zz eta_xz
        %initial values: 0 ? ? 0
        solinit = bvpinit(linspace(0,1,5000),[0 0 0 0],[]);
        SOL = bvp4c(@myode,@bcfun,solinit, bvpset('RelTol', 1e-12,'AbsTol',1e-12));
        toc
    end
    vv=vv(2:size(vv,1),:); %remove the first line
    
    %now vv is an array with x,z,u_x,u_z,jxw,p,T
    dat_calc_ux = vv(:,3);
    dat_calc_uz = vv(:,4);
    
    % compute reference u*:
    dat_ref_ux = interp1(SOL.x, SOL.y(2,:), vv(:,2)); 
    dat_ref_ux = dat_ref_ux .* sin(vv(:,1).*k); % and make it 2d
    dat_ref_uz = interp1(SOL.x, SOL.y(1,:), vv(:,2)); 
    dat_ref_uz = dat_ref_uz .* cos(vv(:,1).*k); % and make it 2d

    % compute |u-u*|_0 of the vector-valued function u=(u_x, u_z)
    % as sqrt( \int (u_x-u*_x)^2 + (u_z-u*_z)^2 )  
    l2error = sqrt( sum( ...
        ( (dat_ref_ux - dat_calc_ux).^2 + (dat_ref_ux - dat_calc_ux).^2 ) ...
          .* vv(:,5) ...
        ) );
    
    convx=[convx; j];
    convy=[convy; l2error];
    
end
figure(1);
loglog(convx,convy,'x-'); hold on;
convorder=polyfit(log(convx),log(convy),1);
convorder=convorder(1);

display([convx convy]);
fprintf('convergence order: %f\n',-convorder);

%%% compute heating terms

Uxs = interp1(SOL.x, SOL.y(2,:), vv(:,2));
Uzs = interp1(SOL.x, SOL.y(1,:), vv(:,2));
eta_xzs = interp1(SOL.x, SOL.y(4,:), vv(:,2));
adab_heating_int = 0;
shear_heating_int = 0;

aa = vv(:,1);

for i=1:size(vv,1)    
    x=vv(i,1);
    z=vv(i,2);
    Ux=Uxs(i);
    Uz=Uzs(i);
    JxW=vv(i,5);
    eta=exp(a*(1-z));
    Ts=0;
    eta_xz=eta_xzs(i)*2*k;
    adab_heating = Di*exp(beta_*(1-z))*Uz*cos(k*x)*(sin(pi*z)*cos(k*x)+Ts);
    % note: this formula is incorrect in the paper. 10/9 is the correct
    %factor, not 4/3. There is also a missing Di/Ra in the sin^2() term.
    shear_heating = Di/Ra * eta * (4*k^2*Ux^2+10/9*beta_^2*Uz^2-4*beta_*k*Ux*Uz) * cos(k*x)^2 + Di/Ra*1/eta*(eta_xz)^2*sin(k*x)^2;
    aa(i) = shear_heating;

    adab_heating_int = adab_heating_int + adab_heating*JxW;
    shear_heating_int = shear_heating_int + shear_heating*JxW;
end

fprintf('heating terms: adab = %.15d, shear = %.15d\n\n',adab_heating_int, shear_heating_int)

end


xlabel('1/h');
ylabel('L2 error of u_x');
xlim([4 256]);
set(gca,'XTick',[8 16 32 64 128]);

ylim([5e-10 5e-5]);
legend('tala a=0','tala a=2', 'ba')

end

function res = bcfun(ya,yb,p)
res = [ya(1);yb(1);ya(4);yb(4)];
end

function [r]=myode(x,Y,p)
% Uz Ux eta_zz eta_xz
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