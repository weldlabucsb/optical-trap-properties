% this code calculates the potential of an optical trap
function [p]=calculateCrossedODTpotential_Li(p)

% load constant parameters
p=loadparams(p);
p=calculatepotential(p);
p=displaypotential(p);

end

function[p]=loadparams(p)

%%%%%%%%%%%%%%
% parameters you might want to vary

% trap parameters
p.lambda=1064E-9; % wavelength of trapping light
p.ODTpower=5/5.5; % power of beam in Watts
p.ODTw=100E-6; % gaussian beam waist in m
p.XODTpower=7/5.5; % power of beam in Watts
p.XODTw=120E-6; % gaussian beam waist in m

% atomic parameters
p.m=7*1.67E-27; % mass of atom
p.lambda_0=670.9E-9;%wavelength of transition in m
p.gamma=2*pi*5.92E6;%linewidth in Radians*Hz
p.fudgefactor=110/70.9; % (should this be different?)
p.fudgefactor=1; % (should this be different?)


%%%%%
% shouldn't need to change parameters below here
%%%%%

%%%%%%%%%%
% physical constants
p.hbar=1.05457E-34;%hbar in SI
p.h=2*pi*p.hbar;
p.kB=1.38E-23;%boltzman constant in SI
p.c=2.99E8;%speed of light in m/s
p.g=9.8;

%%%%%%%%%%%%%%
% dependent parameters
p.ODTz0=pi*p.ODTw^2/p.lambda; % rayleigh range
p.XODTz0=pi*p.XODTw^2/p.lambda; % rayleigh range
p.k=2*pi/p.lambda;
p.freq_0=p.c/p.lambda_0;%frequency of transition in Hz
p.freq=p.c/p.lambda;%frequency of laser in Hz
p.detuning=2*pi*(p.freq_0-p.freq);%detuning in Rad*Hz (nb: positive if red)
p.Er=p.hbar^2*(p.k)^2/(2*p.m);
p.ODTI_avg=p.ODTpower/(pi*p.ODTw^2);
p.XODTI_avg=p.XODTpower/(pi*p.XODTw^2);
p.ODTI_peak=2*p.ODTI_avg;%For a gaussian beam, the peak intensity is twice the average intensity
p.XODTI_peak=2*p.XODTI_avg;%For a gaussian beam, the peak intensity is twice the average intensity
p.ODTpeakpot=p.fudgefactor*3*pi*p.c^2/(2*(2*pi*p.freq_0)^3)*p.gamma/p.detuning*p.ODTI_peak;
p.XODTpeakpot=p.fudgefactor*3*pi*p.c^2/(2*(2*pi*p.freq_0)^3)*p.gamma/p.detuning*p.XODTI_peak;
p.scatt=3*pi*p.c^2/(2*p.hbar*(2*pi*p.freq_0)^3)*(p.gamma/p.detuning)^2*(p.ODTI_peak+p.XODTI_peak);
end

function [p] = calculatepotential(p)
% currently this assumes the beam propagates in the z direction and gravity is in the y direction.

p.numpoints=301;
p.plotrangexy=3; % number of beam waists over which to calculate potential
p.plotrangez=3;  % number of rayleigh lengths over which to calculate potential
p.x=linspace(-p.plotrangexy*p.ODTw,p.plotrangexy*p.ODTw,p.numpoints);
p.y=linspace(-p.plotrangexy*p.ODTw,p.plotrangexy*p.ODTw,p.numpoints);
%p.z=linspace(-p.plotrangez*p.z0,p.plotrangez*p.z0,p.numpoints);
p.z=p.y;

[p.xx p.yy p.zz]=ndgrid(p.x,p.y,p.z);
p.pot=p.m.*p.g.*p.yy - p.ODTpeakpot .* (1./(1+(p.zz./p.ODTz0).^2)) .* exp(-2.*(p.xx.^2+p.yy.^2)./(p.ODTw.^2.*(1+(p.zz./p.ODTz0).^2)))...
                     - p.XODTpeakpot .* (1./(1+(p.xx./p.XODTz0).^2)) .* exp(-2.*(p.zz.^2+p.yy.^2)./(p.XODTw.^2.*(1+(p.xx./p.XODTz0).^2)));

% calculate the trap frequencies and trap center
p.fitvecsize=10;

% do y first since that's the one that might be offset by gravity
p.pot1D_y=p.pot(abs(p.xx)==min(min(min(abs(p.xx)))) & abs(p.zz)==min(min(min(abs(p.zz)))));
ymin=find(p.y==min(min(abs(p.y))));
p.pot1D_y_short=p.pot1D_y((ymin-p.fitvecsize):(ymin+p.fitvecsize));
p.y_short=p.y(ymin-p.fitvecsize:ymin+p.fitvecsize);

% the offset is the trap sag due to gravity
p.sag=p.y_short(p.pot1D_y_short==min(p.pot1D_y_short));
p.minpot=p.pot1D_y_short(p.pot1D_y_short==min(p.pot1D_y_short));

p.pot1D_x=p.pot(p.yy==p.sag & abs(p.zz)==min(min(min(abs(p.zz)))));
xmin=find(p.pot1D_x==min(p.pot1D_x));
p.pot1D_x_short=p.pot1D_x(xmin-p.fitvecsize:xmin+p.fitvecsize);
p.x_short=p.x(xmin-p.fitvecsize:xmin+p.fitvecsize);

p.pot1D_z=p.pot(abs(p.xx)==min(min(min(abs(p.xx)))) & p.yy==p.sag);
zmin=find(p.pot1D_z==min(p.pot1D_z));
p.pot1D_z_short=p.pot1D_z(zmin-p.fitvecsize:zmin+p.fitvecsize);
p.z_short=p.z(zmin-p.fitvecsize:zmin+p.fitvecsize);

p.xfit=fit(1E6*p.x_short',1E9/p.kB*p.pot1D_x_short,'poly2'); % these multipliers are because Princess Fit likes things just so
p.yfit=fit(1E6*p.y_short',1E9/p.kB*p.pot1D_y_short,'poly2');
p.zfit=fit(1E6*p.z_short',1E9/p.kB*p.pot1D_z_short,'poly2');

p.wx=1E6*sqrt(2*(p.kB/1E9)*p.xfit.p1/p.m);
p.wy=1E6*sqrt(2*(p.kB/1E9)*p.yfit.p1/p.m);
p.wz=1E6*sqrt(2*(p.kB/1E9)*p.zfit.p1/p.m);

% now calculate the trap depth with gravity
% start at the y position where we broke the last loop and keep moving down
% until we reach a maximum of the minimum potential at each y position

ii=find(p.y==p.sag); ii=ii-5; % make sure we are down in y from the trap bottom
prevy=p.y(ii);
prevpotvec=p.pot(abs(p.xx)==min(min(min(abs(p.xx)))) & p.yy==prevy);
prevpotmin=min(prevpotvec);
p.maxpot=sqrt(-1);
for jj=(ii-1):-1:1
    thisy=p.y(jj);
    thispotvec=p.pot(p.xx==min(min(abs(p.xx))) & p.yy==thisy);
    thispotmin=min(thispotvec);
    if(thispotmin<prevpotmin)
        p.maxpot=prevpotmin;
        p.maxpoty=prevy;
        break
    else
        prevy=thisy;
        prevpotmin=thispotmin;
    end
end

if(p.maxpot==sqrt(-1))
    beep;
    disp(['Potential maximum not reached -- try increasing plotrangexy for a more accurate trap depth.']);
    p.maxpot=thispotmin;
end

p.trapdepth=p.maxpot-p.minpot;
p.trapdepth_nK=p.trapdepth*1E9/p.kB;
p.trapdepth_rec=p.trapdepth/p.Er;

p.meanfreq=(p.wx*p.wy*p.wz)^(1/3);
p.harmonicosclength=sqrt(p.hbar/p.m/p.meanfreq);

end

function [p] = displaypotential(p)

inds=(abs(p.yy)==min(min(min(abs(p.yy)))));

figure(30);clf;
surf(1E6*squeeze(reshape(p.xx(inds),size(p.xx(1,:,:)))),1E6*squeeze(reshape(p.zz(inds),size(p.zz(1,:,:)))),(1E9/p.kB)*squeeze(reshape(p.pot(inds),size(p.pot(1,:,:)))));
%shading flat;
xlabel('Y Position (microns)');
ylabel('Z Position (microns)');
zlabel('Trap Potential (nK)');
set(gca,'FontSize',18);
view(-160,10);
colormap jet
shading flat

inds=(abs(p.xx)==min(min(min(abs(p.xx)))));

figure(31);clf;
surf(1E6*squeeze(reshape(p.yy(inds),size(p.yy(1,:,:)))),1E6*squeeze(reshape(p.zz(inds),size(p.zz(1,:,:)))),(1E9/p.kB)*squeeze(reshape(p.pot(inds),size(p.pot(1,:,:)))));
%shading flat;
xlabel('Y Position (microns)');
ylabel('Z Position (microns)');
zlabel('Trap Potential (nK)');
set(gca,'FontSize',18);
view(-160,10);
colormap jet
shading flat

figure(32);clf;
subplot(311);
plot(1E6*p.x,1E9/p.kB*p.pot1D_x,'k.');
yl=get(gca,'YLim');
hold on;
plot(1E6*p.x,(p.xfit.p1*(1E6*p.x).^2 + p.xfit.p2*(1E6*p.x) + p.xfit.p3), 'r');
ylim(yl);
legend('Potential',['Quadratic Fit (2\pi\times' num2str(.01*round(100*p.wx/2/pi)) ' Hz)'],'Location','SE');
xlabel('X Position (microns)');
ylabel('Potential (nK)');
set(gca,'FontSize',14);

subplot(312);
plot(1E6*p.y,1E9/p.kB*p.pot1D_y,'k.');
yl=get(gca,'YLim');
hold on;
plot(1E6*p.y,(p.yfit.p1*(1E6*p.y).^2 + p.yfit.p2*(1E6*p.y) + p.yfit.p3), 'r');
ylim(yl);
legend('Potential',['Quadratic Fit (2\pi\times' num2str(.01*round(100*p.wy/2/pi)) ' Hz)'],'Location','SE');
xlabel('Y Position (microns)');
ylabel('Potential (nK)');
set(gca,'FontSize',14);

subplot(313);
plot(1E6*p.z,1E9/p.kB*p.pot1D_z,'k.');
yl=get(gca,'YLim');
hold on;
plot(1E6*p.z,(p.zfit.p1*(1E6*p.z).^2 + p.zfit.p2*(1E6*p.z) + p.zfit.p3), 'r');
ylim(yl);
legend('Potential',['Quadratic Fit (2\pi\times' num2str(.01*round(100*p.wz/2/pi)) ' Hz)'],'Location','SE');
xlabel('Z Position (microns)');
ylabel('Potential (nK)');
set(gca,'FontSize',14);

disp(['ODT Power             : ' num2str(.001*round(1000*p.ODTpower)) ' W']);
disp(['XODT Power            : ' num2str(.001*round(1000*p.XODTpower)) ' W']);
disp(['ODT Beam waist        : ' num2str(.1*round(10*p.ODTw*1E6)) ' um']);
disp(['XODT Beam waist       : ' num2str(.1*round(10*p.XODTw*1E6)) ' um']);
disp(['Wavelength            : ' num2str(.1*round(10*p.lambda*1E9)) ' nm']);
disp(' ');
disp(['Gravitational sag     : ' num2str(.001*round(-1000*p.sag*1E6)) ' um']);
disp(['Trap depth            : ' num2str(.01*round(100*p.trapdepth_nK)) ' nK']);
disp(['Trap depth            : ' num2str(.01*round(100*p.trapdepth_rec)) ' recoils']);
disp(['Trap depth            : ' num2str(.1*round(10*p.trapdepth/p.hbar/p.wx)) ' times omega_x']);
disp(' ');
disp(['X Frequency           : 2pi * ' num2str(.01*round(100*p.wx/2/pi)) ' Hz']);
disp(['Y Frequency           : 2pi * ' num2str(.01*round(100*p.wy/2/pi)) ' Hz']);
disp(['Z Frequency           : 2pi * ' num2str(.01*round(100*p.wz/2/pi)) ' Hz']);
disp(' ');
disp(['Mean Frequency        : 2pi * ' num2str(.01*round(100*p.meanfreq/2/pi)) ' Hz']);
disp(['Harmonic Osc. Length  : ' num2str(.01*round(100*p.harmonicosclength*1E6)) ' um']);



end