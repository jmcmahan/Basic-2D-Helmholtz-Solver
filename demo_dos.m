clear all;
fignum = 1;
n = 64;
k = 2;
x   = @(t) [-0.65 + cos(t) + 0.65*cos(2*t),  1.5*sin(t)];
dx  = @(t) [      - sin(t) - 1.30*sin(2*t),  1.5*cos(t)];
ddx = @(t) [      - cos(t) - 2.60*cos(2*t), -1.5*sin(t)];

r = @(t) sqrt( sum(x(t).^2,2) );
dr = @(t) cos(t) + sin(t);

% Basis functions. Alpha should be a row vector of positive
% values, t should be a column.
cosat = @(alpha,t) cos(t*alpha)/sqrt(pi);
sinat = @(alpha,t) sin(t*alpha)/sqrt(pi);
const = @(t) ones(size(t))/sqrt(2*pi);

% derivatives of basis vectors
dcosat = @(alpha,t) -repmat(alpha, size(t)) .* sinat(alpha,t);
dsinat = @(alpha,t) repmat(alpha, size(t)) .* cosat(alpha,t);
dconst = @(t) zeros(size(t));

% second derivatives of basis vectors
ddcosat = @(alpha,t) -repmat(alpha.^2, size(t)) .* cosat(alpha,t);
ddsinat = @(alpha,t) -repmat(alpha.^2, size(t)) .* sinat(alpha,t);
ddconst = @(t) zeros(size(t));



%x   = @(t) [ sqrt(1)*cos(t),  sqrt(2)*sin(t)];
%dx  = @(t) [-sqrt(1)*sin(t),  sqrt(2)*cos(t)];
%ddx = @(t) [-sqrt(1)*cos(t), -sqrt(2)*sin(t)];

% Angle of incident wave
inc_a = 0;
% Direction of incident wave
inc_d = [cos(pi * inc_a/180); sin(pi * inc_a/180)];
% Boundary condition due to incident wave
f     = @(t) -exp(i*k * x(t)*inc_d); 


obj.x   = x;
obj.dx  = dx; 
obj.ddx = ddx;
obj.f   = f;
obj.n   = n;
obj.k   = k;

disp('Time for solving equation:')
tic;
sol = dos_solve(obj);
toc;

% This is for verifying the same solution listed in Colton & Kress
%angles = linspace(0, pi, 2)';

angles = linspace(0, 2*pi, 129)';
angles = angles(1:end-1);

uinf = sol.find_uinf(angles);
N = 256;
disp('Time for computing far-field:')
tic;
[uinf, nu, v, xhat, phiN] = find_far_field(obj, sol, angles);
toc;

vn = x(sol.solver.t);
xn = vn(:,1);
yn = vn(:,2);
vn = dx(sol.solver.t);
dxn = vn(:,1);
dyn = vn(:,2);
phin = sol.phi;

xN = v(:,1);
yN = v(:,2);
dxN = nu(:,1);
dyN = nu(:,2);

objnew = obj;
%objnew.x   = @(t) [ sqrt(1)*cos(t),  sqrt(2)*sin(t)];
%objnew.dx  = @(t) [-sqrt(1)*sin(t),  sqrt(2)*cos(t)];
%objnew.ddx = @(t) [-sqrt(1)*cos(t), -sqrt(2)*sin(t)];
objnew.x = @(t) 1.01*x(t);
objnew.dx = @(t) 1.01*dx(t);
objnew.ddx = @(t) 1.01*ddx(t);

disp('Time for new boundary update');
tic;
[sol, obj] = sol.update_bdy(objnew);
toc;
disp('Time for new boundary far-field');
tic;
[uinf2, nu2, v2, xhat2, phi2] = find_far_field(obj, sol, angles);
toc;

% Cheesy test of inverse problem - not especially realistic conditions
% since this basis is more-or-less perfect for the obstacle, but we'll
% just see how it does with a basic Newton solve. 

t = sol.solver.t;
t = [t; t(1)];
%t = linspace(0,2*pi,513)';
% Basis for 1st and 2nd coordinates
xfreqs = [1 2];
yfreqs = [1 2];
b1 = [const(t), cosat(xfreqs, t)];
b2 = [sinat(yfreqs, t)];

% Find the coordinate representation of the surface in the given basis
% (this is the exact coordinates)
v = x(t);
x1 = v(:,1); 
x2 = v(:,2);
%c1 = trapz(repmat(t, size([1, xfreqs])), ...
%           b1 .* repmat(x1, size([1, xfreqs]) ))';
%c2 = trapz(repmat(t, size([yfreqs])), ...
%           b2 .* repmat(x2, size([yfreqs]) ))';
c1 = trapz((t), ...
           b1 .* repmat(x1, size([1, xfreqs]) ))';
c2 = trapz((t), ...
           b2 .* repmat(x2, size([yfreqs]) ))';


bdyx = @(c1, c2, t) [[const(t), cosat(xfreqs, t)]*c1, [sinat(yfreqs, t)]*c2];


if 1 == 0
scale = linspace(1, 2, 512)';
for j=1:length(scale)
    objnew.x = @(t) scale(j)*x(t);
    objnew.dx = @(t) scale(j)*dx(t);
    objnew.ddx = @(t) scale(j)*ddx(t);
    [sol, obj] = sol.update_bdy(objnew);
    [uinf2, nu2, v2, xhat2, phi2] = find_far_field(obj, sol, angles);
    plot(sol.solver.t, real(uinf), sol.solver.t, real(uinf2))
    axis([0 2*pi -3 3]);
    %plot(real(uinf), imag(uinf), real(uinf2), imag(uinf2));
    %axis([ -3 3 -3 3]);
    pause(0.01);
end
end

% Activate this to visualize the far-field pattern as a function of
% the scatterer shape

% NOTE: Maybe a better way to do this is to find the average radius
% of the object, make that the radius for the far-field shape, and
% put them on separate subplots. That way they're more comparable. 
% What might be an even better idea would be to find the area of the
% object and determine a radius that makes the other object have
% the same area. Then the transformation from object to scattered
% shape is isometric and the comparison should be easier to visualize.
% This requires approximating the area, though, and is possibly easier
% for some things than others, but since the function is extremely
% smooth, it should work pretty well.
uinforig = uinf;



if 1 == 1
    th = angle(xn + i*yn);
    th(th<0) = th(th<0) + 2*pi;
    A = trapz(th, 0.5*(xn.^2 + yn.^2));
    objeffr = sqrt(A/pi);

    maxur = max(real(uinf));
    minur = min(real(uinf));
    maxui = max(imag(uinf));
    minui = min(imag(uinf));

    uinf = objeffr*(real(uinf) - minur) / (maxur-minur) + ...
            objeffr*i*(imag(uinf) - minui) / (maxui-minui);
    vangs = x(angles);
    xangs = vangs(:,1);
    yangs = vangs(:,2);
    % Maximum radius of object
    objR = sqrt(max( (xangs.^2 + yangs.^2) ));
    % Set the far-field visualization shape so it doesn't intersect
    % the object shape.
    farmin = min(0, min(min(imag(uinf)), min(real(uinf))));
    % The objR*c gives a bit of cushion if needed
    R = objR - farmin + objR*0.1;

    %figure(fignum);
    %fignum = fignum + 1;
    %tryme(R, uinf, sol.solver.t)
    %uinfr = @(R) A - trapz(th, 0.5*(R + real(uinf).^2));
    R = fsolve(@(R) Rsolve(R, uinf, sol.solver.t, A), sqrt(A/pi));
    %R = abs(min(real(uinf)))
    Ri = fsolve(@(R) Rsolve(R, uinf, sol.solver.t, A, 'im'), sqrt(A/pi));

    Rmax = R + abs(max(max(real(uinf)), max(imag(uinf))));
    figure(fignum);
    fignum = fignum + 1;
    uinfx = zeros(length(angles)+1,1);
    uinfy = uinfx;
    % Reference circle
    ucirx = uinfx;
    uciry = uinfx;
    % Imaginary part, just for kicks
    uimgx = uinfx;
    uimgy = uinfx;

    % Visualize far-field as deviation of radius of sphere.
    for j = 1:length(angles)
        uinfx(j) = (R + real(uinf(j)))*cos(angles(j));
        uinfy(j) = (R + real(uinf(j)))*sin(angles(j));
        uimgx(j) = (Ri + imag(uinf(j)))*cos(angles(j));
        uimgy(j) = (Ri + imag(uinf(j)))*sin(angles(j));
        %ucirx(j) = (R+objeffr*0.5)*cos(angles(j));
        %uciry(j) = (R+objeffr*0.5)*sin(angles(j));
        ucirx(j) = (objeffr)*cos(angles(j));
        uciry(j) = (objeffr)*sin(angles(j));
    end
    uinfx(end) = uinfx(1);
    uinfy(end) = uinfy(1);
    uimgx(end) = uimgx(1);
    uimgy(end) = uimgy(1);
    ucirx(end) = ucirx(1);
    uciry(end) = uciry(1);
    plot([xn; xn(1)], [yn; yn(1)], ...
         uinfx, uinfy,...
         ucirx, uciry, '+',...
         uimgx, uimgy,...
         -[inc_d(1)*(Rmax + 0.4*R), inc_d(1)*(Rmax + 0.1*R)], ...
         -[inc_d(2)*(Rmax + 0.4*R), inc_d(2)*(Rmax + 0.1*R)], 'k',...
         -inc_d(1)* (Rmax + 0.1*R), -inc_d(2)*(Rmax + 0.1*R), ...
         'kx', 'LineWidth', 2, 'MarkerSize', 16);
    legend('Object', 'Re(u_\infty)', 'Reference', 'Im(u_\infty)', ...
            'Incident', 'Incident');
end


% Activate this to test for convergence to the far-field.
if 1 == 0

    % Same angles that uinf is computed at
    xdirs = [cos(angles), sin(angles)];
    %xs = R*[cos(angles), sin(angles)];
    tc = sol.solver.t;
    %tc = (0:(2*256-1))'*pi/256;
    j=1;
    Rn = linspace(3, 500, 500);
    errs = zeros(length(Rn),1);
    figure(fignum);
    fignum = fignum + 1;
    disp('Time for verification of far-field convergence');
    tic;
    for R = Rn;
        xs = R*xdirs;
        %disp('Time for computing potential at point list.')
        [uinf_est, xnf, ynf] = sol.find_u_list(xs,tc); 
        % This should retrieve the far-field portion of the potential

        uinf_est = uinf_est / exp(i*k*R) * sqrt(R);
        errs(j) = norm(uinf_est - uinf);
        plot(angles, uinf, angles, uinf_est);
        axis([0 2*pi -2.5 2.5]);
        pause(0.01);
        j = j + 1;
    end
    toc;
end
    
