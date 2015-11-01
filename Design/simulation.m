clear all;

load system;
load lqg_design;
load control_design;
%load simulation;

init_simulator;
sim('simulator');

figure; plot(t,N); title('Noise filter');
figure; plot(t,D); title('Disturbance on control input');

%check
mN = mean(N), mD = mean(D)
rmsN = std(N), rmsD = std(D)
varN = var(N), varD = var(D)

tfDY = feedback(P,F*K1); tfDU = feedback(P*F*K1,1);
tfNY = feedback(F,K1*P); tfNU = feedback(F*K1,P);
y1 = lsim(tfDY,D,t) + lsim(tfNY,N,t);
u1 = lsim(tfDU,D,t) + lsim(tfNU,N,t);

tfDY = feedback(P,F*K05); tfDU = feedback(P*F*K05,1);
tfNY = feedback(F,K05*P); tfNU = feedback(F*K05,P);
y2 = lsim(tfDY,D,t) + lsim(tfNY,N,t);
u2 = lsim(tfDU,D,t) + lsim(tfNU,N,t);

varU1=var(u1), varU2=var(u2)
varY1=var(y1), varY2=var(y2)

figure; 
subplot 211; plot(t,y1); title('y with K1');
subplot 212; plot(t,y2); title('y with K05');
figure; 
subplot 211; plot(t,u1); title('u with K1');
subplot 212; plot(t,u2); title('u with K05');