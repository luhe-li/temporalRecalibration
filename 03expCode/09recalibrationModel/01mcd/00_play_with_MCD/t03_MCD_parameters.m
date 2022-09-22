% This script tests the combination of three taus in MCD, over a small
% range of AV lag (centered around 0), extracts the max mcd_coor and plots

clear all; close all; clc; rng('shuffle');

%% set parameters
fs          = 1e3; % hz
stim_dura   = 0.033; % in sec
ts_duration = 33; % in sec, 7 x 2 padding + 2 stimulus duration

%% loop through temporal constants
% % param=[0.0873, 0.0684, 0.7859]; % tv, ta, tav
% % a broader range based on min max +- se
% tv = linspace(0.0991, 0.1559, 10); %0.0873;
% ta = linspace( 0.0448, 0.1152, 10); %0.0684;
% tav = linspace(0.3104, 1.2086, 10);%0.7859;

% use a broad range of tv, ta, tav
ta = linspace(0.01, 0.2, 20); %0.0684;
tv = linspace(0.01, 0.2, 20); %0.0873;
tav = linspace(0.1, 2, 20);%0.7859;

% use a small range of soa centered around 0
soas = [-1:0.1:1]; % in sec

% initiate corr matrix
corr = NaN(20, 20, 20);

for i = 1:numel(ta)
    ta_i = ta(i);

    for j = 1:numel(tv)
        tv_i = tv(j);

        for v = 1:numel(tav)
            tav_i = tav(v);

            % initiate
            nsample      = ts_duration * fs + 1;
            t            = linspace(0, ts_duration, nsample);

            % low-pass filters (Equation 1)
            fa=fft(t.*exp(-t./ta_i));
            fv=fft(t.*exp(-t./tv_i));
            fav=fft(t.*exp(-t./tav_i));

            % loop 
            parfor ii = 1:numel(soas)
                soa = soas(ii);

                % reconstruct signals

                % initiate empty time seris
                stimV = zeros([1, nsample]); % reset v signals
                stimA = zeros([1, nsample]); % reset a signals
                midpoint = ts_duration/2 * fs + 1; % soa are centered around this midpoint

                % note that soa_m are rounded to be split into half and be centered around the midpoint
                half_soa = round(soa/2 * fs);

                % find out index of a/v stimulus onset and offset
                t_v_onset = midpoint + half_soa;
                t_v_offset = t_v_onset + stim_dura * fs;

                t_a_onset = midpoint - half_soa;
                t_a_offset = t_a_onset + stim_dura * fs;

                % assign 1 to stimulus inpulse
                stimV(t_v_onset:t_v_offset) = 1;
                stimA(t_a_onset:t_a_offset) = 1;

                %     % plot to check signals
                %     figure; hold on
                %     plot(stimA); plot(stimV);

                % early filtering
                st_v=fft(stimV).*fv;
                st_a=fft(stimA).*fa;

                % late filtering
                st_v_av=ifft(st_v.*fav);
                st_a_av=ifft(st_a.*fav);

                % xcorrelate
                u1=st_a_av.*ifft(st_v);         % Equation 2
                u2=st_v_av.*ifft(st_a);         % Equation 3

                % MCD correlation detector output
                MCD_corr_signal=u2.*u1;         % Equation 4
                MCD_corr(ii)=mean(MCD_corr_signal); % Equation 6

                % MCD lag detector output
                MCD_lag_signal=u2-u1;           % Equation 5
                MCD_lag(ii)=mean(MCD_lag_signal);   % Equation 7

            end

            % pick the max corr among the soas for each parameter
            % combination tested
            corr_max = max(MCD_corr);
            corr(i, j, v) = corr_max;
        end
    end

end

%% post-processing
threshold = 1;
bool_corr = corr <= threshold;

%% plot

figure; hold on
filename = "corr.gif";
set(gca, 'LineWidth', 2, 'FontSize', 15)
set(gcf, 'Position',[10 10 900 600])
box off

for v = 1:20
    imagesc(squeeze(bool_corr(:,:,v)))
    colorbar
%     pause(0.02)
    title(['max corr within SOA between [-1:0.1:1] s, tau_{AV} = ' num2str(tav(v))])
    xlabel('tau_{A} = 0.01:0.01:0.2')
    ylabel('tau_{V} = 0.01:0.01:0.2')
    drawnow
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if v == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
end

% %% export GIF
% % 1. Create the initial image file
% gifFile = 'corr.gif';
% exportgraphics(obj, gifFile);
% 
% % 2. Within a loop, append the gif image
% for i = 1:20
%       %   %   %   %   %   %    % 
%       % Update the figure/axes %
%       %   %   %   %   %   %    % 
%       exportgraphics(obj, gifFile, Append=true);
% end
