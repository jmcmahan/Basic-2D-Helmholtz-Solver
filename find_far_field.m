function [uinf, nu, y, xhat, phi] = find_far_field(obj, sol, angles, N)
% Compute the 2D far field solution using the trapezoid rule.
% obj - object structure for the scattering problem
% sol - solution structure for the scattering problem
% angles - list of angles indicating directions to compute the far-field
%          solution at (row vector)
% N - optional parameter indicating number of points to interpolate the
%     solution on in order to compute the trapezoid rule (i.e., 2*N-1).
%     If not set, it is automatically set to the problem's discretization
%     level.


    if nargin < 4
        N = obj.n;
    end
    x   = obj.x;
    dx  = obj.dx;
    k   = obj.k;

    t = (0:2*N-1)'*pi/N;
    % Boundary coordinates used in integration
    y = x(t);
    % Used in computing infinitesimal arc length (which I think is needed?)
    dy = dx(t);
    ds = sqrt(sum(dy.^2, 2));

    % Correct if angles aren't in column format
    s = size(angles);
    if s(1) > s(2)
        angles = angles';
    end

    % Compute angles to unit directions
    xhat = [cos(angles); sin(angles)];

    % Get unit outward normals on the boundary
    nu = unit_normal(obj, t);
    eta = sol.eta;

    % Interpolate to find phi if needed
    if N ~= obj.n
        phi = sol.phiint(t);
    else
        phi = sol.phi;
    end

    uinf = zeros(size(angles'));
    % This is the constant in front of the integral
    a = exp(-i*pi/4) / sqrt(8*pi*k);
    for j = 1:length(xhat)
        xh = xhat(:,j);
        K = (k*nu*xh + eta) .* exp(-i*k*y*xh) .* phi .* ds;
        %uinf(j) = a * trapz(t, K); 
        uinf(j) = a * sum(K * pi/N);
    end
end


function v = unit_normal(obj, t)
% obj - object structure for scattering problem
% t - parameter for the parameterized boundary. Should be a column

    dx = obj.dx;
    s = size(t);
    if s(2) > s(1)
        t = t';
    end
    % Take tangent vector, rotate it 90 degrees to get outward
    % pointing normal.
    v = dx(t) * [0 -1; 1 0];
    r = sqrt(sum(v.^2, 2));
    v = v ./ repmat(r, [1 2]);
end
