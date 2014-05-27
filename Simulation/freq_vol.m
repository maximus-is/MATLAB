function A = freq_vol(s,s0,obs)

%Inputs
%   "obs"
%   s
%   s0

% Output
% Calculate variance at each frequency
N = size(s,1);
num_freq = size(obs,2);

high = s;
low = s;
op = s;
cl = s;
cl_1 = s;

num_est = 8;
alpha = 0.5;

e0 = 1;
e1 = 1.30775226685107;
e2 = 0.0559573802612633;
e3 = -0.00181469471720417;

A = NaN(1,num_est,num_freq);

for i = 1:num_freq
    J = N/prod(obs(1:i));
    
    high = max(reshape(high,obs(i),N/prod(obs(1:i))));
    low = min(reshape(low,obs(i),N/prod(obs(1:i))));
    
    op = reshape(op,obs(i),N/prod(obs(1:i)));
    op = op(1,:);
    
    cl = reshape(cl,obs(i),N/prod(obs(1:i)));
    cl = cl(end,:);
    
    cl_1 = [s0 cl(1:end-1)];
    
    n_obs = size(high,2);
    
    p = log(cl);
    r = p - log(cl_1);
    
    rv = r*r';
    
    ret = p(2:end) - p(1);
    b_ret = p(end) - p(1:end-1);
    
    cub_e = e0 + e1*log(J) + e2*log(J).^2 + e3*log(J).^3;

    cub = exp(cub_e);

    max_log = alpha*var(ret./(1:n_obs-1)) ...
        + (1-alpha)*var(b_ret./(n_obs-1:-1:1));
    
    max_log = max_log*cub;
    
    acc = sum((log(cl./cl_1).^2) - ...
        (log(cl(end)./cl_1(1)).^2)/...
        (n_obs*(n_obs-1)));                                                 % Close-close adjusted for drift
    
    park = (1/(4*log(2)))*log(high./low)*log(high./low)';                   % Parkinson's high-low
    
    gk = 0.5*sum(0.511*(log(high./low).^2) - ...                            % Garman and Klass high-low-open-close
        0.019*log(cl./op).*log((high.*low)./(op.^2)) - ...
        2*log(high./op).*log(low./op));
    
    gk_2 = 2*log(2)*park - (2*log(2)-1)*rv;                                 % Garman and Klass "practical" high-low-open-close
    
    rs = sum(log(high./cl).*log(high./op) + ...                             % Rogers and Satchell high-low-open-close
        log(low./cl).*log(low./op));
    
    yz = yang_zhang(high',low',op',cl',cl_1');                              % Yang and Zhang minimum variance
    
    A(1,:,i) = [acc park gk gk_2 rs yz rv max_log];
end