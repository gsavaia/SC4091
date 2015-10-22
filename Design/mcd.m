% MCD   Multicriteria controller design

RMS2_d=16;              % Squared RMS process noise.
RMS2_n=1;               % Squared RMS sensor noise.
rho=1e-4;
num_p=[0 1 .5];         % Numerator polynomial of the plant.
den_p=[1 -1.6 .7];      % Denomerator polynomial of the plant.

% State space representation of P.
[Ap,Bp,Cp,Dp]=tf2ss(num_p,den_p);

% Initial estimate LQG controller
F=dlqry(Ap,Bp,Cp,Dp,1,rho);           % Feedback gain.
L=dlqe(Ap,Bp,Cp,RMS2_d,RMS2_n);       % Kalman gain.
[Ac,Bc,Cc,Dc]=dreg(Ap,Bp,Cp,Dp,F,L);  % Compute controller.
[num_lqg,den_lqg]=ss2tf(Ac,Bc,Cc,Dc); % Convert to transfer function.
num_q=conv(num_lqg,den_p);
den_q=conv(den_lqg,den_p)+conv(num_lqg,num_p);

% State space representation of Q.
[Aq,Bq,Cq,Dq]=tf2ss(num_q,den_q);

% Solve Lyapunov equation.
At = [ Aq, zeros(4,6), -Bq*Cp
       zeros(4,4), Aq, Bq*Cp, -Bq*Dp*Cp
       zeros(2,8), Ap, -Bp*Cp
       zeros(2,10), Ap ];
Bt = [ -Bq*Dp, -Bq
       -Bq*Dp*Dp, -Bq*Dp
       -Bp*Dp, -Bp
       Bp, zeros(2,1) ];
Wt=diag([RMS2_d, RMS2_n]);
Xt=dlyap(At,Bt*Wt*Bt');

% Initialize optimization procedure and set options.
num_q0=num_q;                        % Initial guess.
opt=optimoptions('fmincon',...
                 'Algorithm','sqp');     % SQP algorithm.
opt=optimoptions(opt,'Display','iter');  % Display results.
opt=optimoptions(opt,'TolX',1e-6);       % Termination tolerance for parameters.
opt=optimoptions(opt,'TolFun',1e-6);     % Termination tolerance for objective
                                         % function.
opt=optimoptions(opt,'TolCon',1e-6);     % Termination tolerance for constraints.
opt=optimoptions(opt,'GradObj','on');    % Use gradient.
opt=optimoptions(opt,'GradConstr','on'); % Use constraint Jacobian.
%
% For a first run of this m-file, you may want to uncomment the following
% two lines, which cause extra diagnostic information to be printed and an
% extra derivative check to be run.
%
%%opt=optimoptions(opt,'Diagnostics','on');
%%opt=optimoptions(opt,'DerivativeCheck','on','FinDiffType','central');

% Compute the trade-off curve.
dv_list=[20:-1:3];                   % List of desired values for M11.
results=zeros(length(dv_list),8);    % Initialization of results.
f_counter=0;                         % Counter for function evaluations.
i_counter=0;                         % Counter for iterations.
for i=1:length(dv_list)
  num_q0=num_q;                      % Initial guess.
  dv_M11=dv_list(i);
  [num_q,J,exitflag(i),output]=...
      fmincon(@(num_q) mcd_J(num_q,den_q,Cp,Dp,Xt,Wt,rho),num_q0,...
              [],[],[],[],[],[],...
              @(num_q) mcd_con(num_q,den_q,dv_M11),opt);
  f_counter=f_counter+output.funcCount;
  i_counter=i_counter+output.iterations;
  infn_M11=mcd_con(num_q,den_q,dv_M11)+dv_M11;
  results(i,:)=[dv_M11,J,infn_M11,num_q];

  disp(['Step ',num2str(i)])
  disp(['Results: ',num2str(results(i,1:2))])
  disp(' ')
end

disp(['Total number of function evaluations: ',num2str(f_counter)]);
disp(['Total number of iterations: ',num2str(i_counter)]);

plot(results(:,1),results(:,2));
xlabel('1/D_{min}');
ylabel('J');
