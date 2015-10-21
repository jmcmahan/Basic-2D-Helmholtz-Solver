fignum = 1;

% **********************************************************************
%                       Setting up the parameters
% **********************************************************************

% The scattering object is represented by a radius function, which is
% a smooth, periodic function giving the radius of the object on some
% parametrization t between 0 and 2*pi. 
% In the code, this is given by a parametrized function:
% x(t) = [ x_coord(t),  y_coord(t) ]
% The function should be written so that if t is a column of M samples
% of the polar angle (theta), then x(t) will be a M by 2 matrix with
% the k^th row [x(t(k)), y(t(k))]. Some examples:
%
% A unit circle:
% x(t) = @(t) [cos(t), sin(t)];
%
% An ellipse stretched in the x-direction:
% x(t) = @(t) [5*cos(t), sin(t)]; 
% 
% You'll also need the first and second derivatives of this function
% with respect to t. For the circle,
%
% dx(t) = @(t) [-sin(t), cos(t)];
% ddx(t) = @(t) [-cos(t), -sin(t)];
%
% The example is a "kite" shape used in Colton and Kress, as well as 
% several papers involving scattering:

x   = @(t) [-0.65 + cos(t) + 0.65*cos(2*t),  1.5*sin(t)];
dx  = @(t) [      - sin(t) - 1.30*sin(2*t),  1.5*cos(t)];
ddx = @(t) [      - cos(t) - 2.60*cos(2*t), -1.5*sin(t)];

% The solver converges exponentially (with respect to n, defined below) 
% for analytic boundaries. One way to handle more complicated boundaries
% while retaining this convergence is to approximate them with a 
% trigonometric series. It should also be possible to base x, dx, and ddx on 
% some sort of (at least twice continuously differentiable) splines and 
% still get reasonable results, it'll just converge more slowly. 

% The solver is for the Helmholtz equation at a single wave number:
% u'' + k^2*u = 0
%
% so you need to input the wave number:

k = 2;

% You also need to provide the value of the incident wave on the boundary.
% This is a frequency-domain based solver, so the easiest things to do
% are plane waves in imaginary form. Here's an example:

% Angle of the incident wave in degrees
inc_ang = 0;    % Points horizontally to the right

% Compute the unit direction from this angle. I am writing it as a column
% since it will be multiplied on the left by a row in the next step.
inc_dir = [cos(pi * inc_ang/180); sin(pi * inc_ang/180)];

% The boundary function for the plane wave with this direction is:
f = @(t) -exp(i*k * x(t)* inc_dir);
% Note this involves the boundary shape function x(t) above. The negative
% sign is needed due to how the problem is set up. The problem
% is linear with respect to the function, f, so it should be possible to
% make f the sum of several plane waves.

% You also need a discretization level for the solver. I don't remember
% if it's necessary, but I always make this a power of two. Larger values
% correspond to more accurate solutions.

n = 32;

% Finally, you can optionally provide a coupling parameter. This has to
% do with the solution representation. If you don't provide it, it will
% be set to k automatically, which was what was suggested in Colton and
% Kress.

%eta = k;

% Stuff our definitions into a structure. Note inc_ang and inc_dir are
% not included, as these were just helper functions for defining f.

% Object shape function
obj.x = x;
obj.dx = dx;
obj.ddx = ddx;

% Boundary function
obj.f = f;

% Wave number
obj.k = k;

% Number of discretization points
obj.n = n;

% Optional coupling parameter
%obj.eta = eta;


% **********************************************************************
%                       Solving for the density
% **********************************************************************

% With the above parameters defined, you compute the boundary density like
% so:

disp('Time for solving for the potential')
tic;
sol = dos_solve(obj);
toc;

% This provides you with the solution structure "sol", which has these
% fields (these are in the comments at the top of  "dos_solve")

% sol.t             Parameter grid that the problem was solved on.
% sol.phi           Solution density, phi
% sol.phiint        Function for interpolating the solution
% sol.find_u        Function for finding the potential at pairs (x,y)
% sol.find_u_list   Function for finding the potential at pairs x(j,:)
% sol.A             Discretized integral operator.
% sol.R             Quadrature matrix, mostly for debug
% sol.g             Discretized right-hand side from the boundary condition
% sol.eta           Eta repeated in case not in obj


% **********************************************************************
%                       Computing the far field
% **********************************************************************

% The far-field depends on angle only. So provide the angles you want
% it solved at in a column vector. 

angles = linspace(0, 2*pi, 129)';
angles = angles(1:end-1);

% Then use this function to get the far field:
disp('Time for computing the far field')
tic;
uinf = find_far_field(obj, sol, angles);
toc;

% You can optionally provide an N after angles, which relates to how
% accurately the integration used to compute the far field is approximated.
% By default it determines this from the discretization in "sol". 

% The real part of the far field potential is what you are interested
% in:

figure(fignum);
fignum = fignum + 1;
plot(180*angles/pi, real(uinf));

% This is the potential, as a function of angle, of the wave field when
% you are far from the scatterer and without the inverse square decay
% exhibited by the field.


% **********************************************************************
%                       Computing the actual field
% **********************************************************************

% You can also compute the potential at arbitrary locations. There are
% two interfaces for this. One allows you to specify the coordinates of
% particular points you want the field at. For more information on that,
% see the comments of the function "find_u_list" in the "dos_solve.m"
% file (last function in the file, currently). You can call it from
% sol.find_u_list(xs, ys, ts). For instance, if you used the ts defined
% below and the values
% xs = [1; 2]; ys = [3; 4]
% the solver would compute the potential at coordinates (1, 3) and (2, 4).


% For the other interface, you specify an x-coordinate
% grid and a y-coordinate grid (column vectors):

xn = 200;
yn = 200;
xs = linspace(-8, 8, xn)';
ys = linspace(-8, 8, yn)';

% You also specify the parameter samples of the boundary used in
% approximating the integrals involved in this. This is similar to how the
% optional N part in the far field computation can compute the solution
% more precisely. For simplicity, you can just use the discretization the
% solver used:

ts = sol.solver.t;

% Now compute the values at each possible combination of xs and ys with

disp('Time for solving potential on grid')
tic;
[v, xv, yv] = sol.find_u(xs, ys, ts);
toc;


% Here v is a column vector of the potential at the grid points and 
% xv, yv, are column vectors indicating the x & y position of the
% sample in each row. For instance, v(1) is the potential at the
% x-coordinate xv(1) and y-coordinate yv(1). You could, for example, do

figure(fignum)
fignum = fignum + 1;
plot3(xv, yv, v, '.');

% although this isn't a great visualization

% This is the quasi-static solution for the scattered field, only. To see 
% it in the time-domain, define a sound-speed

c0 = 1;

% and determine your frequency from this and the wave number

w = k*c0;

% Pick some time values to sample this at. For example, 16 samples over
% two periods:
numperiods = 2;
times = linspace(0, numperiods*(2*pi/w), 16)';

% We'll reshape the solution vector in order to display it as an image:

vm = reshape(v, xn, yn)';

% The time-domain solution is found by modulating the quasi-static solution
% by exp(i*w*times(1)) and taking the real part. 

figure(fignum)
fignum = fignum + 1;
it = 1;
for j = 1:4
    for k = 1:4
        subplot(4,4,it);
        imagesc(xs, ys, real(exp(1i*w*times(it))*vm));
        colormap(gray);
        title(sprintf('Time=%f', times(it)));
        it = it + 1;
    end
end

% Note that the calculations in the interior of the obstacle are not valid.
% Also note that this is only the scattered field. To view the total field
% you have to add in the incident field, as well.


% Calculate the field at the same (x,y) coordinates as the scattered field:
ui = exp(i*k*[xv, yv]*inc_dir);

% Reshape into matrix for display
ui = reshape(ui, xn, yn)';


figure(fignum)
fignum = fignum + 1;
it = 1;
for j = 1:4
    for k = 1:4
        subplot(4,4,it);
        % Note the addition of the incident field 
        imagesc(xs, ys, real(exp(1i*w*times(it))*(vm + ui)));
        colormap(gray);
        title(sprintf('Time=%f', times(it)));
        it = it + 1;
    end
end        


% It's easier to see in motion

doit = input('Show field in motion? (y/n)', 's')
if doit == 'y' | doit == 'Y'
    

    figure(fignum);
    fignum = fignum + 1;
    timesnew = linspace(times(1), times(end), 128);
    colormap(gray);

    for j = 1:length(timesnew);
        subplot(1,2,1);
        tds = exp(1i*w*timesnew(j))*vm;
        tdi = exp(1i*w*timesnew(j))*ui;
        imagesc(xs, ys, real(tds + tdi), [-2 2]);
        title('Total');
        subplot(1,2,2);
        imagesc(xs, ys, real(tds), [-1.5 1.5]);
        
        title('Scattered');
        pause(0.1);
    end;
end
