function [J] = gradient_cost_function(x, u, ys)

J = zeros(size(u,1), size(x,1)); %J = 2000x6
    
J(:,1) = lsim( -tf([0 0 1], [1 x(2) x(3)], -1), u );
J(:,2) = lsim( tf([0 0 1 x(1) 0], [1 2*x(2) 2*x(3)+x(2)^2 2*x(2)*x(3) x(3)^2], -1), u);
J(:,3) = lsim( tf([0 0 0 1 x(1)], [1 2*x(2) 2*x(3)+x(2)^2 2*x(2)*x(3) x(3)^2], -1), u);

J(:,4) = lsim( -tf([1 x(6) 0], [x(4)^2 2*x(5)*x(4) x(5)^2], -1), ys);
J(:,5) = lsim( -tf([0 1 x(6)], [x(4)^2 2*x(5)*x(4) x(5)^2], -1), ys);
J(:,6) = lsim( tf([0 1],[x(4) x(5)], -1), ys);

end
