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
V0 = (0.01:0.01:0.1)/tra_d;
num_vol = length(V0);

s0 = 100;
s0_gbm = s0;
sf_gbm = s0;

mu = (1 + r*dt);

alpha = 0.5;

enum = [20:5:195,200:50:1950,2000:500:19500, 20000:5000:195000];%, 200000:50000:1950000, 2000000:500000:20000000];
rat = NaN(num_vol,length(enum));

var_gbm = NaN(length(enum),num_vol,num_days);
tic;
for d = 1:num_days
    for k = 1:num_vol
        s_gbm = NaN(N,1);
        s_gbm(1) = s0_gbm;
        Z_gbm = sqrt(V0(k))*sqrt(dt)*randn(N,1);
        
        for i = 2:N
            s_gbm(i,:) = mu*s_gbm(i-1) + s_gbm(i-1)*Z_gbm(i-1);
        end
        
        for j = 1:length(enum)
            
            idx = round([1, linspace(N/enum(j), N, enum(j))]);
            
            r = log(s_gbm(idx(2:end))) - log(s_gbm(1,:));                      % Calculate the expanding return
            back_r = log(s_gbm(end,:)) - log(s_gbm(idx(1:end-1),:));             % Calculate the backwards expanding return
            
            var_gbm(j,k,d) = alpha*var(r./[enum(j):-1:1]') ...
                + (1-alpha)*var(back_r./[1:enum(j)]');
            fprintf('d = %.4d - j = %.3d  ',d,j);
            
        end
        
        sf_gbm = s_gbm(end);
        s0_gbm = mu*sf_gbm + sf_gbm.*Z_gbm(end);
    end
end
toc;
% %%
% iter = 0;
% for J = enum
%     fprintf('\nJ = %d\n',J)
%     iter = iter + 1;
%     
%     times = [1:J]'/J;
%     
%     %min_p = NaN(M,1);
%     vala = NaN(num_est,num_bm,M);
%     
%     parfor k = 1:M
%         
%         paths = zeros(J,num_bm);
%         
%         for i = 1:num_bm
%             obj = gbm(daily_drift/J, daily_vol(i)/sqrt(J),'StartState',p0);
%             paths(:,i) = obj.simByEuler(J-1);
%         end
%         
%         Y = paths;
%         %min_p(k) = min(min(Y));
%         
%         %vala(:,:,k) = max_vol([1 J],Y,0.5,'log');
%         
%         vala(:,:,k) = vol_est_intra(Y,times,sub_samp);
%     end
%     
%     mean_var = mean(vala,3);
%     ratio(:,:,iter) = repmat(daily_vol.^2,8,1)./mean_var;
% end
% 
% %fixedInterval = seconds2wall(wall2seconds(93000):300:wall2seconds(160000));
% %RV = realized_variance(Y(:,1),times,'unit','CalendarTime',1/J)
% 
% msgbox('Done')
% 
% %% Plots
% est_names = {'Close-Close','Adjusted Close-Close','Parkinsons High-Low',...
%     'Garman and Klass High-Low-Open-Close','Rogers and Satchell',...
%     'Yang and Zhang','Realized Volatility','Max'};
% est_fig_names = {'cc','acc','p','gk','rs','yz','rv','max'};
% path = pwd;
% 
% for i = 1:num_est
%     figure; hold on;
%     for j = 1:num_bm
%         scatter(enum,squeeze(ratio(i,j,:)),'b.')
%     end
%     title(sprintf('%s',est_names{i}));
%     outputname = [path '\fig\fig_ratio_' sprintf('%s',est_fig_names{i}) '.eps'];
%     %print('-dpsc2',outputname)
% end
% 
% close all;
% 
% diary off;
