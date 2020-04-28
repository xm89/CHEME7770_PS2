%% Problem 1b
clear all; close all; clc;
% assume k1=k2=k3=k4=k,for first consideration, let k=k1=0.1
% gamma1*R_T/V2=5;
% gamma3*X_T/V4=10;
% V1/V2=(gamma1*R_T/V2)*thetaB;
% V3/V4=(gamma3*X_T/V4)*xstar;
% Assume k_D=0.05

k_D=0.05;k1=0.1;
for i=1:100 
thetaB(i)=1/(1+k_D(i));
xstar(i)=((5.*thetaB(i).*k1+k1+1-5.*thetaB(i))-sqrt((5.*thetaB(i)-5.*thetaB(i).*k1-k1-1).^2-4.*(1-5.*thetaB(i)).*(5.*thetaB(i).*k1)))./(2.*(1-5.*thetaB(i)));
% xstar=((1+k1+V1/V2*k2-V1/V2)+-sqrt((V1/V2-V1/V2*k2-k1-1)^2-4*(1-V1/V2)*(V1/V2*k2)))/(2*(1-V1/V2));
ystar(i)=((1+k1+10.*xstar(i).*k1-10.*xstar(i))-sqrt((10.*xstar(i)-10.*xstar(i).*k1-k1-1).^2-4.*(1-10.*xstar(i)).*(10.*xstar(i).*k1)))./(2.*(1-10.*xstar(i)));
% ystar=(1+k3+V3/V4*k4-V3/V4)+sqrt((V3/V4-V3/V4*k4-k3-1)^2-4*(1-V3/V4)*(V3/V4*k4))/(2*(1-V3/V4))
k_D(i+1)=k_D(i)+0.05;
end
kd=k_D(1:length(k_D)-1);
% plot figure at k=0.1
figure (1)
plot(1./kd,thetaB,1./kd,xstar,1./kd,ystar);
legend('thetaB','xstar','ystar');
xlabel('1/k_D');
title('PS 1b for k=0.1');

% assume k1=k2=k3=k4=k,for first consideration, let k=k2=10
k2=10;
k_D=0.05;
for i=1:100
thetaB2(i)=1/(1+k_D(i));
xstar2(i)=((5.*thetaB2(i).*k2+k2+1-5.*thetaB2(i))-sqrt((5.*thetaB2(i)-5.*thetaB2(i).*k2-k2-1).^2-4.*(1-5.*thetaB2(i)).*(5.*thetaB2(i).*k2)))./(2.*(1-5.*thetaB2(i)));
% xstar=((1+k1+V1/V2*k2-V1/V2)-sqrt((V1/V2-V1/V2*k2-k1-1)^2-4*(1-V1/V2)*(V1/V2*k2)))/(2*(1-V1/V2));
ystar2(i)=((1+k2+10.*xstar2(i).*k2-10.*xstar2(i))-sqrt((10.*xstar2(i)-10.*xstar2(i).*k2-k2-1).^2-4.*(1-10.*xstar2(i)).*(10.*xstar2(i).*k2)))./(2.*(1-10.*xstar2(i)));
% ystar=(1+k3+V3/V4*k4-V3/V4)+sqrt((V3/V4-V3/V4*k4-k3-1)^2-4*(1-V3/V4)*(V3/V4*k4))/(2*(1-V3/V4))
k_D(i+1)=k_D(i)+0.05;
end
kd=k_D(1:length(k_D)-1);
% plot figure at k=10
figure (2)
plot(1./kd,thetaB2,1./kd,xstar2,1./kd,ystar2);
legend('thetaB','xstar','ystar');
xlabel('1/k_D');
title('PS 1b for k=10');