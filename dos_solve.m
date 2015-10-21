function sol = dos_solve(obj)
% Solve the 2D direct object scattering problem according to the parameters
% passed in obj using the Nystrom method described in Colton & Kress. 
% The structure "obj" contains the following members:
%
% obj.x:    parameterization of 2D object boundary. This should be a
%           function which takes a column vector of parameters and
%           returns 2 columns for of the same length as the parameter
%           column vector. In other words, 
%           x([1;2;3]) = [x_1(1), x_2(1); x_1(2), x_2(2); x_1(3), x_2(3)];
% obj.dx:   parameterization of the derivative of 2D object boundary. 
%           Should be in same format as obj.x
% obj.ddx:  parameterization of 2nd derivative of 2D object boundary.
%           Should be in same format as obj.x
% obj.f     Boundary condition for the problem. I.e., on the boundary
%           of the scattering object, the field should equal this function.
%           Should be a function in the same format as obj.x but with only
%           one column.
% obj.n     Number of points to divide the problem into using (2*n-1).
% obj.k     Wave number of the problem.
% obj.eta   Coupling parameter for the potential. Suggestion is eta=k
%           which will be used if this member doesn't exist.
%
% The solution is returned in the structure "sol". It has the following
% members:
%
% sol.t             Parameter grid that the problem was solved on.
% sol.phi           Solution density, phi
% sol.phiint        Function for interpolating the solution
% sol.find_u        Function for finding the potential at pairs (x,y)
% sol.find_u_list   Function for finding the potential at pairs x(j,:)
% sol.A             Discretized integral operator.
% sol.R             Quadrature matrix, mostly for debug
% sol.g             Discretized right-hand side from the boundary condition
% sol.eta           Eta repeated in case not in obj

    if ~isfield(obj, 'eta')
        obj.eta = obj.k;
    end

    [phi, solver] = solve_for_phi(obj);
    sol.phi = phi;
    sol.solver = solver;
    sol.eta = obj.eta;
    sol.phiint = @(t) phiint(t, obj, sol);
    sol.find_u = @(x, y, t) find_u(x, y, obj, sol, t);
    sol.find_u_list = @(x, t) find_u_list(x, obj, sol, t);
    sol.find_uinf = @(angles) find_far_field(obj, sol, angles);
    sol.update_bdy = @(objnew) solve_for_phi_new_bdy(sol, obj, objnew);
end


function y = H1_0(x)
% Helper function to compute the 0th order Hankel function of the first 
% kind with notation closer to the given formulas.
    y = besselh(0, 1, x); 
end

function y = H1_1(x)
% Helper function to compute the 1st order Hankel function of the first 
% kind with notation closer to the given formulas.
    y = besselh(1, 1, x); 
end

function y = J_0(x)
% Helper function to compute the 0th order Bessel function with notation 
% closer to the given formulas.
    y = besselj(0, x); 
end

function y = J_1(x)
% Helper function to compute the 1st order Bessel function with notation 
% closer to the given formulas.
    y = besselj(1, x); 
end



function [k, n, x, dx, ddx, f] = get_param(obj)
% Helper function to retrieve variables from the obj structure.
% [k, n, x, dx, ddx, f] is the order they're retrieved in.
    k = obj.k;
    n = obj.n;
    x = obj.x;
    dx = obj.dx;
    ddx = obj.ddx;
    f = obj.f;
end

function [t, tau, matsize] = make_param_vecs(t, tau)
% Helper function to create vectors which form all subtraction pairs 
% for computation to be reshaped into the discretized kernel matrix.
% SHOULD BE CALLED BEFORE CALLING KERNEL FUNCTIONS!!!!
% 
% t - discretized independent variable (column)
% tau - discretized integration variable (row)
% matsize - size vector for reshaping later

    s = size(t);
    if (s(2) > s(1))
        t = t';
    end
    s = size(tau);
    if (s(1) > s(2))
        tau = tau';
    end
    l1 = length(t);
    l2 = length(tau);

    % Repeat so vector
    s = size(t);
    t = repmat(t, size(tau));
    tau = repmat(tau, s);

    % Put parameter indices in vector form for easier calculation
    matsize = size(t);
    vecsize = [matsize(1)*matsize(2) 1];

    tau = reshape(tau, vecsize);
    t = reshape(t, vecsize);
end


function [v, si] = Ldiag(t, tau, obj)
% Compute the diagonal terms of L2 and L

    [k, n, x, dx, ddx] = get_param(obj);

    % Find singular indices, replace with deduced diagonal values
    si = (t == tau);
    tsi = t(si);

    dxt = dx(tsi);
    % Rotated second derivative
    rotdxt = ddx(tsi) * [0 -1; 1 0];
    % Squared norm of dx(t)
    nrmdxt2 = sum(dxt.^2, 2);
    v = 1/(2*pi) * sum(dxt .* rotdxt, 2) ./ nrmdxt2;
end


function v = L(t, tau, obj)
% Compute the integral kernel L from Colton and Kress. SHOULD NOT BE CALLED
% DIRECTLY!
%
% t - OUTPUT FROM make_param_vecs!!!
% tau - OUTPUT FROM make_param_vecs!!!
% obj - scattering object structure 

    [k, n, x, dx, ddx] = get_param(obj);

    xtau_t = x(tau) - x(t);
    % Rotate derivative by 90 degrees
    rotdxtau = dx(tau) * [0 -1; 1 0];
    % Compute norm of x(t) - x(tau) at each param point
    nrmxt_tau =  sqrt(sum(xtau_t.^2, 2)); 
    v = i*k/2 * sum(rotdxtau.*xtau_t, 2) .* H1_1(k*nrmxt_tau) ./ nrmxt_tau;

    % Replace diagonal values 
    [diagv, si] = Ldiag(t, tau, obj);
    v(si) = diagv;
    
end


function v = L1(t, tau, obj)
% Compute the integral kernel L1 from Colton and Kress. SHOULD NOT BE CALLED
% DIRECTLY!
%
% t - OUTPUT FROM make_param_vecs!!!
% tau - OUTPUT FROM make_param_vecs!!!
% obj - scattering object structure 

    [k, n, x, dx, ddx] = get_param(obj);

    xt_tau = x(t) - x(tau);
    % Rotate derivative by 90 degrees
    rotdxtau = dx(tau) * [0 -1; 1 0];
    % Compute norm of x(t) - x(tau) at each param point
    nrmxt_tau =  sqrt(sum(xt_tau.^2, 2)); 
    v = k/(2*pi) * sum(rotdxtau.*xt_tau, 2) .* J_1(k*nrmxt_tau) ./ nrmxt_tau;

    % Replace diagonal values 
    % WARNING - I am assuming L1 should be 0 on the diagonal without much
    % thought. It's not explicitly stated in Colton & Kress, but it looks
    % like that may be correct.
    [diagv, si] = Ldiag(t, tau, obj);
    v(si) = zeros(size(diagv));
end


function v = L2(t, tau, obj)
% Compute the integral kernel L2 from Colton and Kress. SHOULD NOT BE CALLED
% DIRECTLY!
%
% t - OUTPUT FROM make_param_vecs!!!
% tau - OUTPUT FROM make_param_vecs!!!
% obj - scattering object structure 

    [k, n, x, dx, ddx] = get_param(obj);
    v = L(t, tau, obj) - L1(t, tau, obj).*log(4*sin((t-tau)/2).^2);
    
    % Replace diagonal values 
    [diagv, si] = Ldiag(t, tau, obj);
    v(si) = diagv;
end


function [v, si] = Mdiag(t, tau, obj)
% Compute the diagonal terms of M2.

    [k, n, x, dx, ddx] = get_param(obj);

    % Find singular indices, replace with deduced diagonal values
    si = (t == tau);
    tsi = t(si);

    dxt = dx(tsi);
    nrmdxt  = sqrt(sum(dxt.^2, 2));

    % "Euler's constant"
    C = 0.57721;
    v = (i/2 - C/pi - 1/pi * log(k/2 * nrmdxt)).*nrmdxt;
end


function v = M(t, tau, obj)
% Compute the integral kernel M from Colton and Kress. SHOULD NOT BE CALLED
% DIRECTLY!
%
% t - OUTPUT FROM make_param_vecs!!!
% tau - OUTPUT FROM make_param_vecs!!!
% obj - scattering object structure 

    [k, n, x, dx, ddx] = get_param(obj);

    xt_tau = x(t) - x(tau);
    nrmxt_tau = sqrt(sum(xt_tau.^2, 2));
    dxtau = dx(tau);
    nrmdxtau = sqrt(sum(dxtau.^2, 2));

    v = i/2 * H1_0(k*nrmxt_tau) .* nrmdxtau;
end


function v = M1(t, tau, obj)
% Compute the integral kernel M1 from Colton and Kress. SHOULD NOT BE CALLED
% DIRECTLY!
%
% t - OUTPUT FROM make_param_vecs!!!
% tau - OUTPUT FROM make_param_vecs!!!
% obj - scattering object structure 

    [k, n, x, dx, ddx] = get_param(obj);

    xt_tau = x(t) - x(tau);
    nrmxt_tau = sqrt(sum(xt_tau.^2, 2));
    dxtau = dx(tau);
    nrmdxtau = sqrt(sum(dxtau.^2, 2));

    v = -1/(2*pi) * J_0(k*nrmxt_tau) .* nrmdxtau;
end


function v = M2(t, tau, obj)
% Compute the integral kernel M2 from Colton and Kress. SHOULD NOT BE CALLED
% DIRECTLY!
%
% t - OUTPUT FROM make_param_vecs!!!
% tau - OUTPUT FROM make_param_vecs!!!
% obj - scattering object structure 

    [k, n, x, dx, ddx] = get_param(obj);

    v = M(t, tau, obj) - M1(t, tau, obj) .* log(4*sin((t-tau) / 2).^2);

    % Replace diagonal values 
    [diagv, si] = Mdiag(t, tau, obj);
    v(si) = diagv;
end

function [K1, K2] = getK(t, obj)
% Compute the discretized kernels for computing the density.
% t - vector of values to evaluate independent variable at
% obj - structure for scattering problem

    [k, n, x, dx, ddx] = get_param(obj);
    eta = obj.eta;

    % This is the variable which is integrated over.
    tau = (0:(2*n-1)) * pi/n;

    % Create vector of all pairs of t and tau for easy computation
    [t, tau, s] = make_param_vecs(t, tau);

    K1 = L1(t, tau, obj) + i*eta*M1(t, tau, obj);
    K2 = L2(t, tau, obj) + i*eta*M2(t, tau, obj);

    K1 = reshape(K1, s);
    K2 = reshape(K2, s);
end

function g = getg(t, obj)
% Discretize the right-hand-side of the integral equation for the
% phi. 
% t - column vector of points to evaluate the density at
%
    s = size(t);
    if s(2) > s(1)
        t = t';
    end
    [k, n, x, dx, ddx, f] = get_param(obj);

    g = 2*f(t);
end


function solver = setup_solver(obj)
    [k, n, x, dx, ddx] = get_param(obj);
    t = (0:(2*n-1))' * pi/n;

    [K1, K2] = getK(t, obj);
    R = dos_quad(n, t);
    A = R.*K1 + pi/n*K2;
    g = getg(t, obj);

    % Finished matrix
    OP = eye(length(A)) - A;
    solver.OP = OP;
    solver.t = t;
    solver.g = g;
    solver.R = R;
    solver.A = A;
end

function [phi, solver] = solve_for_phi(obj)
    solver = setup_solver(obj);

    phi = solver.OP \ solver.g; 

end

function [sol, obj] = solve_for_phi_new_bdy(sol, obj, objnew)
% Use to update the solution for a new boundary without rebuilding the
% discretization and what not. For use with inverse solvers. 
%
% obj - old object structure
% objnew - new object structure, should be same as obj except with
%          updated x, dx, and ddx
% sol - previous solution
% 
    solver = sol.solver;
    t = solver.t;
    % Update the boundary info
    g = getg(t, objnew);
    solver.g = g;
    % Solve for new density
    phi = solver.OP \ solver.g;
    % Update solution structure
    sol.solver = solver;
    sol.phi = phi;
    sol.phiint = @(t) phiint(t, obj, sol);
    sol.find_u = @(x, y, t) find_u(x, y, obj, sol, t);
    sol.find_u_list = @(x, t) find_u_list(x, obj, sol, t);
    sol.find_uinf = @(angles) find_far_field(obj, sol, angles);
    sol.update_bdy = @(objnew) solve_for_phi_new_bdy(sol, obj, objnew);
    obj = objnew;
end

function phii = phiint(t, obj, sol)
    phi = sol.phi;
    n = obj.n; 
    [K1, K2] = getK(t, obj);
    R = dos_quad(n, t); 
    A = R.*K1 + pi/n*K2;
    g = getg(t, obj);

    phii = A*phi + g;
end

function v = Loffbdy(xv, tau, obj)
% Compute the integral kernel L from Colton and Kress for values
% x off the boundary. SHOULD NOT BE CALLED  DIRECTLY!
%
% xv - OUTPUT FROM make_param_vecs!!! These are the points to compute
%      at
% tau - OUTPUT FROM make_param_vecs!!! These are the integration
%      variable discretized points
% obj - scattering object structure 

    [k, n, x, dx, ddx] = get_param(obj);

    xtau_t = x(tau) - xv;
    % Rotate derivative by 90 degrees
    rotdxtau = dx(tau) * [0 -1; 1 0];
    % Compute norm of x(t) - x(tau) at each param point
    nrmxt_tau =  sqrt(sum(xtau_t.^2, 2)); 

    % Note - factor of 2 put back in (i.e., i*k/4 not i*k/2). 
    % The numerical treatment in 
    % Colton & Kress writes these operators so the 2 needed for
    % the density equation is included.
    v = -i*k/4 * sum(rotdxtau.*xtau_t, 2) .* H1_1(k*nrmxt_tau) ./ nrmxt_tau;

    % There should be no diagonal values since xv must be off the
    % boundary.
    %
    % Replace diagonal values 
    %[diagv, si] = Ldiag(t, tau, obj);
    %v(si) = diagv;
    
end

function v = Moffbdy(xv, tau, obj)
% Compute the integral kernel L from Colton and Kress for values
% x off the boundary. SHOULD NOT BE CALLED  DIRECTLY!
%
% xv - OUTPUT FROM make_param_vecs!!! These are the points to compute
%      at
% tau - OUTPUT FROM make_param_vecs!!! These are the integration
%      variable discretized points
% obj - scattering object structure 

    [k, n, x, dx, ddx] = get_param(obj);

    xt_tau = xv - x(tau);
    nrmxt_tau = sqrt(sum(xt_tau.^2, 2));
    dxtau = dx(tau);
    nrmdxtau = sqrt(sum(dxtau.^2, 2));

    % Note - factor of 2 put back in. The numerical treatment in 
    % Colton & Kress writes these operators so the 2 needed for
    % the density equation is included. Also, put back in the negative
    % sign which is added there when constructing the discretized
    % integral operator (i.e., I - A, whereas here, we're only
    % interested in A, in fact 2B = A).
    v =  -i/4 * H1_0(k*nrmxt_tau) .* nrmdxtau;
end


function [v, xv, yv] = find_u(x, y, obj, sol, t)
% Compute the potential u from the solution phi at the points
% (x,y) with phi being interpolated on t.
% 
% x - vector of x-values to compute 
% y - vector of y-values to compute
% Note: all pairs of (x(i), y(j)) are computed
% obj - scattering object structure
% sol - solution structure
% t - optional parameter vector to interpolate phi on. If not
%     used then the default is used.
%
% v - solution, u, evaluated at pairs (x,y)
% xv - x-coordinate of pairs, (x,y)
% yv - y-coordinate of pairs, (x,y)

    if length(t) == length(sol.solver.t) 
        if t == sol.solver.t
            phi = sol.phi;
        else
            phi = phiint(t, obj, sol);
        end
    else
        phi = phiint(t, obj, sol);
    end
    n = length(t)/2;
    %[k, n, x, dx, ddx] = get_param(obj);
    eta = obj.eta;
    % Make all pairs of x and y into a vector
    [xv, yv, matsize1] = make_param_vecs(x, y);

    % Now create pairs of the (x,v) pairs and the parameter
    % vector t (store in tauv - tv not used).
    [xvv, tauv, s1] = make_param_vecs(xv, t);
    [yvv, tauv, s2] = make_param_vecs(yv, t);
    vn = [xvv, yvv];

    % Compute the discretized kernels
    Koffbdy = (Loffbdy(vn, tauv, obj) + i*eta*Moffbdy(vn, tauv, obj))*pi/n;

    % Put them back in a matrix form
    Koffbdy = reshape(Koffbdy, s1);

    % Multiply the matrix by the density, which is the integration
    % needed to find u at the pairs, (x,y)
    v = Koffbdy * phi;
end

function [v, xv, yv] = find_u_list(x, obj, sol, t)
% Compute the potential u from the solution phi at the points
% x (which is a Nx2 matrix of (x,y) values to evaluate
% at) with phi being interpolated on t.
% 
% x - Nx2 matrix of N (x,y) values to compute u on
% obj - scattering object structure
% sol - solution structure
% t - optional parameter vector to interpolate phi on. If not
%     used then the default is used.
%
% v - solution, u, evaluated at pairs (x(j,1),x(j,2)), all j
% xv - x-coordinate, i.e., x(:,1)
% yv - y-coordinate, i.e., x(:,2)

    if length(t) == length(sol.solver.t) 
        if t == sol.solver.t
            phi = sol.phi;
        else
            phi = phiint(t, obj, sol);
        end
    else
        phi = phiint(t, obj, sol);
    end
    n = length(t)/2;
    %[k, n, x, dx, ddx] = get_param(obj);
    eta = obj.eta;

    s = size(x);
    if (s(2) ~= 2)
        disp('Error using find_u_list: x should be an Nx2 matrix.');
    end
    xv = x(:,1);
    yv = x(:,2);
    % Now create pairs of the (x,v) pairs and the parameter
    % vector t (store in tauv - tv not used).
    [xvv, tauv, s1] = make_param_vecs(xv, t);
    [yvv, tauv, s2] = make_param_vecs(yv, t);
    vn = [xvv, yvv];

    % Compute the discretized kernels
    Koffbdy = (Loffbdy(vn, tauv, obj) + i*eta*Moffbdy(vn, tauv, obj))*pi/n;

    % Put them back in a matrix form
    Koffbdy = reshape(Koffbdy, s1);

    % Multiply the matrix by the density, which is the integration
    % needed to find u at the pairs, (x,y)
    v = Koffbdy * phi;
end
