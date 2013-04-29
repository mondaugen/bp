% Solve Basis Pursuit Denoising as described by CDS (1996)
clear;
% Length of signals in dictionary
N = 8;
% How overdetermined the dictionary is
L = 2;
% Size of dictionary
P = N*L;

% The dictionary
Phi = fourier_dict(N,L);

% The signal to decompose
% Secret coefficients so that our signal is build from dictionary (for
% testing)
scoeffs = zeros(P,1);
scoeffs(2) = 0.5;
sceoffs(6) = 0.75;
s = Phi*scoeffs;

% The signal concatenated with a negated copy
b = s;

% The dictionary concatenated with a negated copy
A = [Phi -Phi];

% The regularization parameter lambda
lambda = 0.01;

% The coefficient vector associated therewith
c = lambda*ones(2*P,1);

% Log-barrier parameter
mu = 0.5;

% Regularization parameters
gamma = 1e-4;
delta = 1e-4;

% Compute initial guess for x, y, z as described by Sardy, Bruce and Tseng
% (1998)

% Use ridge estimate (or Tikhonov regularization) to get a closed form
% initial guess
alpha_ridge = (Phi'*Phi + lambda*eye(P))\(Phi'*b);

alpha_plus = max(alpha_ridge,0);
alpha_minus = max(-alpha_ridge,0);
x = [alpha_plus; alpha_minus] + (0.1)*ones(2*P,1);
q  = A*sign(x);
w  = 1.1*norm(Phi'*q,inf);
y = (lambda/w)*q;
z = lambda*ones(2*P,1) - ((A')*y);

% Max number of iterations
K = 50;

for k=(1:K),
    t = c - (gamma^2)*x - z - A'*y;
    r = b - A*x - (delta^2)*y;
    v = mu*ones(2*P,1) - z.*x;
    D = inv(diag(z./x) + (gamma^2)*eye(2*P));

    % This system is better solved approximately using conjugate gradient
    % methods
    dy = (A*D*A' + (gamma^2)*eye(N))\(r - A*D*(v./x - t));

    dx = D*A'*dy + D*(v./x - t);
    dz = v./x - (z./x).*dx;

    % Calculate primal and dual step sizes
    pp = 0.99*p_max_find(x,dx);
    pd = 0.99*p_max_find(z,dz);

    % Update variables
    x = x + pp*dx;
    y = y + pd*dy;
    z = z + pd*dz;

    mu = (1 - min([pp pd 0.99]))*mu;
    
    % Terminate if following conditions are satisfied
    if ((norm(r,2)/(1 + norm(x,2)) < gamma) ...
            && (norm(t,2)/(1 + norm(y,2)) < gamma) ...
            && (z'*t/(1 + norm(z,2)*norm(x,2)) < delta)),
        break
    end
end
display(x);


