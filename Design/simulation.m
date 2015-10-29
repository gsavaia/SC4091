clear all;

load system;
load lqg_design;
load control_design;

RMSd = 4;
RMSn = 1;

t = (1:1:2000)';
N = sqrt(RMSn)*randn(2000,1);
D = sqrt(RMSd)*randn(2000,1);

%check
mN = mean(N), mD = mean(D)
varN = var(N), varD = var(D)


%% Dmin = 1;
x = 0.5;
Q1 = tf( NumQ_minima( Dmin_inv==x, : ), DenQ, -1 );
K1 = feedback(Q1, -P*F);
Q1
disp(['J = ', num2str(J( Dmin_inv==x ) )]);

%% Dmin = 0.5 (=> 1/Dmin = 2)
x = 2;
Q05 = tf( NumQ_minima( Dmin_inv==x, : ), DenQ, -1 );
K05 = feedback(Q05, -P*F);
Q05
disp(['J = ', num2str(J( Dmin_inv==x ))]);

