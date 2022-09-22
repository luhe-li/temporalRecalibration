% t05_MCD_tv_changes_lag

% this scripts trys to reproduce supplementary figure 8(c)

clear all; close all; clc; rng(1);

k = 2/3:0.02:3/2;
para_default = [0.0873, 0.0684, 0.7859]; % default tv, ta, tav
tvs = para_default(1).*k;
tas = para_default(2)./k;
tavs1 = para_default(3).*k;
tavs2 = para_default(3)./k;

%% reproduce supp fig 8a

clm = gray;
clm_idx = floor(linspace(1,200,numel(k)));

% plot the default filter
x = 0:0.1:5;

for i = 1:numel(k)
    tv = tvs(i);
    ta = tas(i);
    tav1 = tavs1(i);
    tav2 = tavs2(i);
    fv(i,:) = x.*exp(-x./tv);
    fa(i,:) = x.*exp(-x./ta);
    fav1(i,:) = x.*exp(-x./tav1);
    fav2(i,:) = x.*exp(-x./tav2);
end

clm_i = floor(linspace(1,250,numel(k)));
clm = gray;

% when tau gets larger, color gets lighter
% when k gets larger, ta gets larger (lighter), tv gets smaller(darker)
figure;
sgtitle('original exponential function');
for i = 1:numel(k)

    subplot 221; hold on
    plot(x, fa(i,:), 'Color', clm(clm_i(i),:));
    title('visual filter')

    subplot 223; hold on
    plot(x, fv(i,:), 'Color', clm(clm_i(i),:))
    title('auditory filter')

    subplot 222; hold on
    plot(x, fav1(i,:), 'Color', clm(clm_i(i),:));
    title('visual - multisensory filter')

    subplot 224; hold on
    plot(x, fav2(i,:), 'Color', clm(clm_i(i),:));
    title('auditory - multisensory filter')

end

% plot the default parameters
for i = 1:numel(para_default)
    filter(i,:) = x.*exp(-x./(para_default(i)));
end

f1 = subplot(2,2,1);
plot(x, filter(1,:),'m','LineWidth',2);
f2 = subplot(2,2,2);
plot(x, filter(3,:),'m','LineWidth',2);
f3 = subplot(2,2,3);
plot(x, filter(2,:),'c','LineWidth',2);
f4 = subplot(2,2,4);
plot(x, filter(3,:),'c','LineWidth',2);
linkaxes([f1, f2, f3, f4], 'xy')

%% normalize all distribution by dividing its peak?

fv = fv./sum(fv, 2);
fa = fa./sum(fa, 2);
fav1 = fav1./sum(fav1, 2);
fav2 = fav2./sum(fav2, 2);
filter = filter./sum(filter, 2);

% when tau gets larger, color gets lighter
% when k gets larger, ta gets larger (lighter), tv gets smaller(darker)
figure;
sgtitle('normalized exponential function');
for i = 1:numel(k)

    subplot 221; hold on
    plot(x, fa(i,:), 'Color', clm(clm_i(i),:));
    title('visual filter')

    subplot 223; hold on
    plot(x, fv(i,:), 'Color', clm(clm_i(i),:))
    title('auditory filter')

    subplot 222; hold on
    plot(x, fav1(i,:), 'Color', clm(clm_i(i),:));
    title('visual - multisensory filter')

    subplot 224; hold on
    plot(x, fav2(i,:), 'Color', clm(clm_i(i),:));
    title('auditory - multisensory filter')

end

% plot the default parameters
f1 = subplot(2,2,1);
plot(x, filter(1,:),'m','LineWidth',2);
f2 = subplot(2,2,2);
plot(x, filter(3,:),'m','LineWidth',2);
f3 = subplot(2,2,3);
plot(x, filter(2,:),'c','LineWidth',2);
f4 = subplot(2,2,4);
plot(x, filter(3,:),'c','LineWidth',2);
linkaxes([f1, f2, f3, f4], 'xy')

%% scale all distribution by dividing its peak?

fv = fv./max(fv, [], 2);
fa = fa./max(fa, [], 2);
fav1 = fav1./max(fav1, [], 2);
fav2 = fav2./max(fav2, [], 2);
filter = filter./max(filter, [], 2);

% when tau gets larger, color gets lighter
% when k gets larger, ta gets larger (lighter), tv gets smaller(darker)
figure;
sgtitle('scaled exponential function');
for i = 1:numel(k)

    subplot 221; hold on
    plot(x, fa(i,:), 'Color', clm(clm_i(i),:));
    title('visual filter')

    subplot 223; hold on
    plot(x, fv(i,:), 'Color', clm(clm_i(i),:))
    title('auditory filter')

    subplot 222; hold on
    plot(x, fav1(i,:), 'Color', clm(clm_i(i),:));
    title('visual - multisensory filter')

    subplot 224; hold on
    plot(x, fav2(i,:), 'Color', clm(clm_i(i),:));
    title('auditory - multisensory filter')

end

% plot the default parameters
f1 = subplot(2,2,1);
plot(x, filter(1,:),'m','LineWidth',2);
f2 = subplot(2,2,2);
plot(x, filter(3,:),'m','LineWidth',2);
f3 = subplot(2,2,3);
plot(x, filter(2,:),'c','LineWidth',2);
f4 = subplot(2,2,4);
plot(x, filter(3,:),'c','LineWidth',2);
linkaxes([f1, f2, f3, f4], 'xy')

%% reproduce figure bc

params = [tvs; tas; tavs1; tavs2]';
% add default parameters to the end
params = [params; 0.0873, 0.0684, 0.7859, 0.7859];

% set parameters
fs                          = 1e3; % hz
stim_dura                   = 0.033; % in sec
ts_duration                 = 16; % in sec, 7 x 2 padding + 2 stimulus duration

% for a range of soa
soas                        = [-2:0.1:2]; % in sec

% initiate
figure; hold on
[all_MCD_corr, all_MCD_lag] = deal(NaN(size(params, 1), numel(soas)));

for j = 1:size(params, 1)
    % initiate
    nsample                     = ts_duration * fs + 1;
    t                           = linspace(0, ts_duration, nsample);

    % low-pass filters (Equation 1)
    param                       = params(j,:);
    fv                          = fft(t.*exp(-t./param(1)));
    fa                          = fft(t.*exp(-t./param(2)));
    fav1                         = fft(t.*exp(-t./param(3)));
    fav2                         = fft(t.*exp(-t./param(4)));

    for i                       = 1:numel(soas)

        % flip the sign
        soa                         = -soas(i);

        %% reconstruct signals

        % initiate empty time seris
        stimV                       = zeros([1, nsample]); % reset v signals
        stimA                       = zeros([1, nsample]); % reset v signals
        midpoint                    = ts_duration/2 * fs + 1; % soa are centered around this midpoint

        % note that soa_m are rounded to be split into half and be centered around the midpoint
        half_soa                    = round(soa/2 * fs);

        % find out index of a/v stimulus onset and offset
        t_v_onset                   = midpoint + half_soa;
        t_v_offset                  = t_v_onset + stim_dura * fs;

        t_a_onset                   = midpoint - half_soa;
        t_a_offset                  = t_a_onset + stim_dura * fs;

        % assign 1 to stimulus inpulse
        stimV(t_v_onset:t_v_offset) = 1;
        stimA(t_a_onset:t_a_offset) = 1;

        %     % plot to check signals
        %     figure; hold on
        %     plot(stimA); plot(stimV);

        % early filtering
        st_v                        = fft(stimV).*fv;
        st_a                        = fft(stimA).*fa;

        % late filtering
        st_v_av                     = ifft(st_v.*fav1);
        st_a_av                     = ifft(st_a.*fav2);

        % xcorrelate
        u1                          = st_a_av.*ifft(st_v);         % Equation 2
        u2                          = st_v_av.*ifft(st_a);         % Equation 3

        % MCD correlation detector output
        MCD_corr_signal             = u2.*u1;         % Equation 4
        MCD_corr(i)                 = mean(MCD_corr_signal); % Equation 6

        % MCD lag detector output
        MCD_lag_signal              = u2-u1;           % Equation 5
        MCD_lag(i)                  = mean(MCD_lag_signal);   % Equation 7

    end

    % save mcd output by different k
    all_MCD_corr (j,:) = MCD_corr;
    all_MCD_lag (j,:) = MCD_lag;

    if j ~= size(params,1)
    %% plotting
    subplot(1, 2, 1); hold on
    plot(soas, MCD_corr,'LineWidth',1, 'Color', clm(clm_i(j),:));

    subplot(1, 2, 2); hold on
    plot(soas, MCD_lag,'LineWidth',1,  'Color', clm(clm_i(j),:));

    else % the default parameter
    subplot(1, 2, 1); hold on
    set(gca, 'LineWidth', 1, 'FontSize', 10)
    set(gcf, 'Position',[10 10 500 400])
    xlabel('SOA (s)')
    title('corr')
    plot(soas, MCD_corr,'LineWidth',1.5, 'Color', 'r');

    subplot(1, 2, 2); hold on
    set(gca, 'LineWidth', 1, 'FontSize', 10)
    set(gcf, 'Position',[10 10 500 400])
    xlabel('SOA (s)')
    title('lag')
    plot(soas, MCD_lag,'LineWidth',1.5, 'Color', 'r');

    end
end

%% plot delta_lag against soa
delta_MCD_lag = all_MCD_lag(1:(end-1),:) - all_MCD_lag(end,:);
figure; hold on
set(gca, 'LineWidth', 1, 'FontSize', 10)
set(gcf, 'Position',[10 10 500 400])
for j = 1:size(delta_MCD_lag, 1)
    plot(soas, delta_MCD_lag(j,:),'LineWidth',1.5,'Color', clm(clm_i(j),:));
end
yline(0,'r')
ylabel('\Delta_{lag}')
xlabel('SOA(s)')