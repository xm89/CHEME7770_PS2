%% Problem 1d 
% consider 1/KD=0.1 at k1=0.1
KD=1/0.1; k1=0.1;
thetaB1d=1/(1+KD);
xstar1d=((5.*thetaB1d.*k1+k1+1-5.*thetaB1d)-((5.*thetaB1d-5.*thetaB1d.*k1-k1-1).^2-4.*(1-5.*thetaB1d).*(5.*thetaB1d.*k1)).^(0.5))./(2.*(1-5.*thetaB1d));
ystar1d=((1+k1+10.*xstar1d.*k1-10.*xstar1d)-((10.*xstar1d-10.*xstar1d.*k1-k1-1).^2-4.*(1-10.*xstar1d).*(10.*xstar1d.*k1)).^(0.5))./(2.*(1-10.*xstar1d));
% consider 1/KD=0.15 at k1=0.1
KD2=1/0.15;
thetaB1d2=1/(1+KD2);
xstar1d2=((5.*thetaB1d2.*k1+k1+1-5.*thetaB1d2)-((5.*thetaB1d2-5.*thetaB1d2.*k1-k1-1).^2-4.*(1-5.*thetaB1d2).*(5.*thetaB1d2.*k1)).^(0.5))./(2.*(1-5.*thetaB1d2));
ystar1d2=((1+k1+10.*xstar1d2.*k1-10.*xstar1d2)-((10.*xstar1d2-10.*xstar1d2.*k1-k1-1).^2-4.*(1-10.*xstar1d2).*(10.*xstar1d2.*k1)).^(0.5))./(2.*(1-10.*xstar1d2));
% percentage change of thetaB, xstar, and ystar for k=0.1 respectivly:
thetaB_ch=abs((thetaB1d-thetaB1d2)./thetaB1d).*100
xstar_ch=abs((xstar1d-xstar1d2)./xstar1d).*100
ystar_ch=abs((ystar1d-ystar1d2)./ystar1d).*100

% consider 1/KD=0.1 at k1=10
KD=1/0.1; k2=10;
thetaB1d=1/(1+KD);
xstar1d=((5.*thetaB1d.*k2+k2+1-5.*thetaB1d)-((5.*thetaB1d-5.*thetaB1d.*k2-k2-1).^2-4.*(1-5.*thetaB1d).*(5.*thetaB1d.*k2)).^(0.5))./(2.*(1-5.*thetaB1d));
ystar1d=((1+k2+10.*xstar1d.*k2-10.*xstar1d)-((10.*xstar1d-10.*xstar1d.*k2-k2-1).^2-4.*(1-10.*xstar1d).*(10.*xstar1d.*k2)).^(0.5))./(2.*(1-10.*xstar1d));
% consider 1/KD=0.15 at k1=10
KD2=1/0.15;
thetaB1d2=1/(1+KD2);
xstar1d2=((5.*thetaB1d2.*k2+k2+1-5.*thetaB1d2)-((5.*thetaB1d2-5.*thetaB1d2.*k2-k2-1).^2-4.*(1-5.*thetaB1d2).*(5.*thetaB1d2.*k2)).^(0.5))./(2.*(1-5.*thetaB1d2));
ystar1d2=((1+k2+10.*xstar1d2.*k2-10.*xstar1d2)-((10.*xstar1d2-10.*xstar1d2.*k2-k2-1).^2-4.*(1-10.*xstar1d2).*(10.*xstar1d2.*k2)).^(0.5))./(2.*(1-10.*xstar1d2));
% percentage change of thetaB, xstar, and ystar for k=10 respectivly:
thetaB_ch2=abs((thetaB1d-thetaB1d2)./thetaB1d).*100
xstar_ch2=abs((xstar1d-xstar1d2)./xstar1d).*100
ystar_ch2=abs((ystar1d-ystar1d2)./ystar1d).*100