% Design Optimization
% wyyjtu@gmail.com (2012-11-09)
function optimization()
    clear
    clc
    % initialize
    p1 = [2900, 200, 1700];
    p2 = [6200, 0, 1500];
    p3 = [5900, 2600, 1200];
    p4 = [6100, -1700, 1100];
    alpha = 135;
    beta = 135;
    % find a constrained minimum of the objective function
    x0 = [0.5, 0.5];
    f = @(x)obj_fun(p1, p2, p3, p4, alpha, beta, x);
    options = optimset('Algorithm', 'sqp', 'TolX', 1e-9);
    [x, fval] = fmincon(f, x0, [], [], [], [], [0, 0], [1, 1], [], options);
    % output
    obj_fun(p1, p2, p3, p4, alpha, beta, x, 1);
    % graph
    graph_objective(p1, p2, p3, p4, alpha, beta);
end
% objective function
function val = obj_fun(p1, p2, p3, p4, alpha, beta, x, varargin)
    % m = x(1) = L_a1 / L_12
    % n = x(2) = L_b4 / L_34
    m = x(1);
    n = x(2);
    % nodes pa, pb
    pa = m * p2 + (1 - m) * p1;
    pb = n * p3 + (1 - n) * p4;
    % vectors
    if m <= 0
        p1 = (m - 1) * p2 + (2 - m) * p1;
    end
    if n <= 0
        p4 = (n - 1) * p3 + (2 - n) * p4;
    end
    pa_pb = pb - pa;
    pa_p1 = p1 - pa;
    pb_p4 = p4 - pb;
    % alpha, beta
    alpha2 = acos(dot(pa_p1, pa_pb) / (norm(pa_p1) * norm(pa_pb)));
    beta2 = acos(dot(pb_p4, -pa_pb) / (norm(pb_p4) * norm(pa_pb)));
    alpha2 = alpha2 / pi * 180;
    beta2 = beta2 / pi * 180;
    % the value of the objective function
    val = ((alpha - alpha2)^2 + (beta - beta2)^2)^(1/2);
    % output
    if length(varargin) == 1
        fprintf('When (m, n) = (%15.12f, %15.12f):\nf(m, n) = %f\n\n', x, val)
        fprintf('pa = (%f, %f, %f)\npb = (%f, %f, %f)\n\n', pa, pb)
        fprintf('alpha = %f\nbeta = %f\n\n', alpha2, beta2)
    end
end
% graph
function graph_objective(p1, p2, p3, p4, alpha, beta)
    [X, Y] = meshgrid(0:0.05:1);
    Z = zeros(size(X));
    for i = 1:size(Z,1)
        for j = 1:size(Z,2)
            mn(1) = X(i, j);
            mn(2) = Y(i, j);
            Z(i, j) = obj_fun(p1, p2, p3, p4, alpha, beta, mn);
        end
    end
    surf(X, Y, Z)
    xlabel('m = L_{a1} / L_{12}'), ylabel('n = L_{b4} / L_{34}')
    zlabel('f(m, n)')
end
