clear, clc
% Methods
% In this file, I used policy function iteration to solve the model. The
% grid search codes are commented out below the PFI part.

% PARAMETERS
alpha = 1/3; % Cobb-Douglas production function parameter
beta = .99; % discount factor
sigma = 2; % coefficient of risk aversion
rho = 0.5; % coefficient of log(z_t)
sigma_z = 0.2; % standard deviation of epsilon in AR-1 of log(z_t)
delta = 0.025; % depreciation rate

m = 5;
[lnz, PI] = rouwenhorst(rho, sigma_z, m);
z = exp(lnz');

tol = 1;
% i = 0; % iteration counter
while (tol > 1e-6)
    PI_in = PI * PI;
    tol = max(max(abs(PI - PI_in)));
    PI = PI_in;
    %     i = i + 1
end

N_s = z' * PI_in(1, :)'; % compute the aggregate effective labor supply
N_d = N_s;

a_lo = 0; % lower bound of grid points
a_hi = 80; % upper bound of grid points
% How to decide the size of asset vector?
num_a = 500;
a = linspace(a_lo, a_hi, num_a); % asset (row) vector

K_max = 50;
K_min = 20;

% ITERATE OVER ASSET PRICES
PI_num = 30;
aggsav = 1 ;
K_tol = 1;
while abs(K_tol) >= .01
    
    K_guess = (K_min + K_max) / 2;
    % Calculate the factor prices
    r = alpha * (K_guess / N_d) ^ (alpha - 1) + (1 - delta);
    w = (1 - alpha) * (K_guess / N_d) ^ alpha;
    
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, r * a', a);
    cons = bsxfun(@plus, cons, permute(z' * w, [1 3 2])); % permute(z, [1 3 2]) * w?
    ret = (cons .^ (1-sigma)) ./ (1-sigma); % current period utility
    ret(cons<0) = -Inf;
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(m, num_a);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%policy function iteration %%%%
    v_tol = 1;
    while v_tol >.0001;
   % CONSTRUCT TOTAL RETURN FUNCTION
    v_mat = ret + beta * ...
           repmat(permute(PI * v_guess, [3 2 1]), [num_a 1 1]);
   
   % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
    [vfn, pol_indx] = max(v_mat, [], 2);
    vfn = permute(vfn, [3 1 2]);
    pol_indx = permute(pol_indx, [3 1 2]);
    v_tol = abs(max(v_guess(:) - vfn(:)));
    v_guess = vfn;
    Q = makeQmatrix(pol_indx, PI);
    % update value functions
    
    % I can't efficiently transform pol_indx into the dimension I want.
    % However repeat the following command for 3 times gives me what I
    % want.
    pol_indx = permute(pol_indx, [3 1 2]);
    pol_indx = permute(pol_indx, [3 1 2]);
    pol_indx = permute(pol_indx, [3 1 2]);
    
    % KEEP DECSISION RULE
    pol_fn = a(pol_indx);
    cons_pi = bsxfun(@minus, r * a, pol_fn);
    cons_pi = bsxfun(@plus, cons_pi, w*z); 
    ret_pi = (cons_pi .^ (1-sigma)) ./ (1-sigma); % current period utility
    ret_pi_v = ret_pi(:);
    ret_pi_con = v_guess(:);
        for pii = 1:PI_num
            ret_pi_new = ret_pi_v+beta*Q*ret_pi_con;
            ret_pi_con = ret_pi_new;
        end;
    v_guess = reshape(ret_pi_con, m, num_a);
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%Grid search%%%%%%%%%%
    %v_tol = 1;
    %while v_tol >.0001;
   % CONSTRUCT TOTAL RETURN FUNCTION
     %   v_mat = ret + beta * ...
    %     repmat(permute(PI * v_guess, [3 2 1]), [num_a 1 1]);
   
   % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
     %[vfn, pol_indx] = max(v_mat, [], 2);
    % vfn = permute(vfn, [3 1 2]);
   
    % v_tol = abs(max(v_guess(:) - vfn(:)));
    
    % v_guess = vfn; % update value functions
    %end;
    %pol_indx = permute(pol_indx, [3 1 2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % KEEP DECSISION RULE
    pol_fn = a(pol_indx);
    % SET UP INITITAL DISTRIBUTION
    Mu = ones(m, num_a); % alternative initial guess: same mass in all states
    Mu = Mu / sum(Mu(:)); % normalize total mass to 1
    
    % ITERATE OVER DISTRIBUTIONS
    % loop over all non-zeros states
    mu_tol = 1;
    while mu_tol > 1e-08
        [z_ind, a_ind] = find(Mu > 0); % find non-zero indices
        MuNew = zeros(size(Mu));
        for ii = 1:length(z_ind)
            apr_ind = pol_indx(z_ind(ii), a_ind(ii));
            MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
                (PI(z_ind(ii), :) * Mu(z_ind(ii), a_ind(ii)) )';
        end
        mu_tol = max(abs(MuNew(:) - Mu(:)));
        Mu = MuNew ;
    end
    
    % CHECK AGGREGATE DEMAND
    aggsav = sum( pol_fn(:) .* Mu(:) ); % Aggregate future assets
    K_tol = aggsav - K_guess;
    if K_tol > 0 ;
        K_min = K_guess ;
    end ;
    if K_tol < 0;
        K_max = K_guess ;
    end ;
    
    % I can't really approach the tolerence, but K should be around 30.336
    % based on my limited times of iteration.
    display (['K = ', num2str(K_guess)])
    display (['Aggregate desired wealth = ', num2str(aggsav)]);
    display (['New Kmin is ', num2str(K_min), ', new Kmax is ', num2str(K_max)]);
    display (['New K is ', num2str((K_max + K_min)/2)]);
    display (['Tolerance is ', num2str(K_tol)]);
    display (' ') ;
    
end

% interest rate in complete market
r_cm = 1/beta;

% plot the value function
plot(a, vfn)


% plot the policy function
figure
plot(a, pol_fn)


