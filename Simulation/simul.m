% Simulation parameters
num_est = 7;
num_freq = 8;
freq = [100,5,2,2,3,5,2,2];

% Time parameters
sec_day = 86400;                                                            % Number of seconds in each day
ticks_sec = 100;                                                            % Number of observations per second
tra_d = 252;                                                                % Number of trading days
N = sec_day*ticks_sec;                                                      % Number of observations per day
r = 0.2/(252*N);                                                            % Annual return
r = 0.2/tra_d;                                                              % Annual return
dt = 1/N;                                                                   % Step size
num_days = 1;

% Initial parameters
%V0 = ((100:200:900).^2)/(252*N);                                            % Initial variance
V0 = (0.04:0.02:0.1)/tra_d;
s0 = 100;                                                                   % Initial stock price

% Variance Estimates
var_gbm = NaN(num_days,num_est,num_freq);
var_garch = var_gbm;
var_h = var_gbm;

% Actual variances and stock prices
true_gbm = NaN(num_days,2);
true_garch = true_gbm;
true_h = true_gbm;

%for k = 1
    
    v0 = V0(k);
    % GARCH(1,1) parameters
    alpha1 = 0.4;                                                           % Coefficient of innovations
    beta = 0.5;                                                             % Coefficient of past variances
    alpha0 = v0*(1-alpha1-beta);                                            % Constant: Set so that the long run variance is equal to v0
    
    % Daily starting values
    v0_garch = v0;
    v0_h = v0;
    v_gbm = v0;                                                             % GBM has constant variance
    
    s0_gbm = s0;
    s0_garch = s0;
    s0_h = s0;
    
    % Daily final values
    sf_gbm = s0;
    sf_garch = s0;
    sf_h = s0;
    
    true_gbm(:,2) = ones(num_days,1)*v0;
    
    for d = 1:num_days
        
        fprintf('k = %d, d = %d\n',k,d)
        % Heston parameters
        rho = -0.7;                                                         % Correlation between the two wiener processes
        kappa = 60.21;%6.21;                                                    % Speed of mean-variance reversion
        xi = 0.1;%0.61;                                                  % Variance of variance
        theta = v0;                                                         % Long-run variance, set so it equals the initial
        
        % Pre-allocation of variances
        v_garch = NaN(N,1);
        v_h = v_garch;
        
        % Pre-allocation of stock prices
        s_gbm = v_garch;
        s_garch = v_garch;
        s_h = v_garch;
        
        % Initialisation
        v_garch(1) = v0_garch;
        v_h(1) = v0_h;
        
        s_gbm(1) = s0_gbm;
        s_garch(1) = s0_garch;
        s_h(1) = s0_h;
        
        % GBM innovations
        Z_gbm = sqrt(v_gbm)*sqrt(dt)*randn(N,1);
        
        % GARCH(1,1) innovations
        Z_g = randn(N,1).^2;                                                % Generate GARCH innovations
        Z_w = sqrt(dt)*randn(N,1);                                          % Generate wiener steps
        
        % Heston innoavtions
        Z_s = sqrt(dt)*randn(N,1);                                          % Generate wiener steps for Heston process as
        Z_v = xi*sqrt(dt)*(rho*Z_s + sqrt(1-rho*rho)*randn(N,1));           % correlated normally distributed variables
        
        % Common operations
        mu = (1 + r*dt);
        
        % Generate price paths
        for i = 2:N
            s_gbm(i) = mu*s_gbm(i-1) + s_gbm(i-1)*Z_gbm(i-1);
            
            v_garch(i) = alpha0 + alpha1*v_garch(i-1)*Z_g(i-1) + beta*v_garch(i-1);
            s_garch(i) = mu*s_garch(i-1) + sqrt(v_garch(i-1))*s_garch(i-1)*Z_w(i-1);
            
            v_h(i) = v_h(i-1) + kappa*(theta - max(v_h(i-1),0))*dt + sqrt(max(v_h(i-1),0))*Z_v(i-1);
            s_h(i) = mu*s_h(i-1) + sqrt(max(v_h(i-1),0))*s_h(i-1)*Z_s(i-1);
        end
        
        % Store the closing stock price and the true variance for that
        % day
        true_gbm(d,1) = s_gbm(end);
        true_garch(d,:) = [s_garch(end) mean(v_garch)];
        true_h(d,:) = [s_h(end) mean(v_h)];
        
        % Calculate variance at each frequency
        var_gbm(d,:,:) = freq_vol(s_gbm,s0_gbm,freq);
        var_garch(d,:,:) = freq_vol(s_garch,s0_garch,freq);
        var_h(d,:,:) = freq_vol(s_h,s0_h,freq);
        
        %[k d V0(k) real_vol([1 N],s_gbm) mean(v_garch) real_vol([1 N],s_garch) mean(v_h) real_vol([1 N],s_h)]
        
        % Store final values 
        sf_gbm = s_gbm(end);
        sf_garch = s_garch(end);
        sf_h = s_h(end);
        
        % Calculate initial values for next day
        s0_gbm = mu*sf_gbm + sf_gbm*Z_gbm(end);
        
        v0_garch = alpha0 + alpha1*v_garch(end)*Z_g(end) + beta*v_garch(end);
        s0_garch = mu*sf_garch + sqrt(v_garch(end))*sf_garch*Z_w(end);
        
        v0_h = v_h(end) + kappa*(theta - max(v_h(end),0))*dt + sqrt(max(v_h(end),0))*Z_v(end);
        s0_h = mu*sf_h + sqrt(max(v_h(end),0))*sf_h*Z_s(end);
        
    end
%end

%save(sprintf('sim%d.mat',k))