function Video_fMRI_any(fMRI_signal, TR, file_label, opts)

%%%%%
%
%  VIDEO OF fMRI IN SELECTED BAND
%  Bandpass filter (optional) the fMRI signal and generate a video of its
%  temporal evolution. Works for a single slice or a full volume.
%
%  Inputs:
%   fMRI_signal : X x Y x Z x T array. fMRI data is always 4D; a single
%                 slice is simply an array where one of X, Y, Z equals 1.
%   TR          : repetition time, in seconds
%   file_label  : base name used for the saved video file
%   opts        : struct with fields
%                   band_pass, high_pass, low_pass
%                   save_video, Figures_and_Videos_folder, video_acceleration
%                   select_colormap
%
%  scripts by Joana Cabral, July 2026
%  joanabcabral@tecnico.ulisboa.pt
%
%%%%%%%%%%%

[X_size, Y_size, Z_size, Tmax]=size(fMRI_signal);
fMRI_signal=double(reshape(fMRI_signal,[X_size*Y_size*Z_size, Tmax]));
fMRI_signal=fMRI_signal-mean(fMRI_signal,2); % Remove the mean for filter

if opts.band_pass
    disp(['    Now bandpass filtering ' num2str(opts.high_pass) '-' num2str(opts.low_pass) ' Hz'])
    fMRI_signal = bandpass_fft(fMRI_signal', [opts.high_pass opts.low_pass], 1/TR)';
    Boundary=round(10/TR); % To cut the first and last 10 seconds of the signals after band pass
    fMRI_signal=fMRI_signal(:,Boundary:end-Boundary);
    Tmax=size(fMRI_signal,2);
end

fMRI_signal=reshape(fMRI_signal,[X_size, Y_size, Z_size, Tmax]);

% fMRI data is always 4D (X,Y,Z,T). A single-slice ("2D") scan is simply
% one where one of the first three dimensions equals 1 - could be X, Y or
% Z depending on the slice orientation (SAG/COR/AX), so we detect it
% directly from the data shape instead of a manual flag.
is_volume = ~any([X_size Y_size Z_size] == 1);

%% Generate video of signals
if ~is_volume
    figure('Position',[158 103 1133 884])
else
    figure('units','normalized','outerposition',[0 0 1 1]) % full screen
end
colormap(opts.select_colormap)
colorlimitbar=5*std(fMRI_signal(:));

if opts.save_video
    if ~opts.band_pass
        videoModes = VideoWriter([opts.Figures_and_Videos_folder '/' file_label],'MPEG-4');
    else
        % This is to include the frequency band in the video name without dots
        high_pass_label=num2str(opts.high_pass);
        high_pass_label(high_pass_label=='.')='p';
        low_pass_label=num2str(opts.low_pass);
        low_pass_label(low_pass_label=='.')='p';
        videoModes = VideoWriter([opts.Figures_and_Videos_folder '/' file_label '_Filt' high_pass_label '-' low_pass_label],'MPEG-4');
    end
    videoModes.FrameRate = round(opts.video_acceleration*1/TR);
    videoModes.Quality = 100;
    open(videoModes);
end

if ~is_volume
% Run this if just one slice
% squeeze(fMRI_signal(:,:,:,t)) collapses whichever spatial dimension is
% singleton, regardless of orientation (SAG/COR/AX)

for t=1:Tmax

    imagesc(imresize(squeeze(fMRI_signal(:,:,:,t)),2)',[-colorlimitbar colorlimitbar])
    axis image
    axis off
    axis xy
    title(['T= ' num2str(t*TR,'%2f') ' secs'])
    colormap(opts.select_colormap)  % reapply after imagesc resets it
    drawnow                    % force complete rendering before capture
    if opts.save_video
        frame = print(gcf, '-RGBImage', '-r0');
        [h, w, ~] = size(frame);
        h2 = ceil(h/16)*16;
        w2 = ceil(w/16)*16;
        padded = 255*ones(h2, w2, 3, 'uint8');   % white pad, matches figure background
        padded(1:h, 1:w, :) = frame;
        writeVideo(videoModes, padded);
    else
        pause(0.1)
    end
end

else
% Run this if you have a 3D volume over time

n_slices=3;

for t=1:Tmax
    for s=1:n_slices
        subplot(n_slices,3,(s-1)*3+1)
        imagesc(squeeze(fMRI_signal(:,:,round(Z_size/(n_slices+3)*(s+1)),t)),[-colorlimitbar/2 colorlimitbar/2])
        axis image
        axis off
        if s==1
        title(['T= ' num2str(round(t*TR*10)/10,'%2f') ' secs'])
        end

        subplot(n_slices,3,(s-1)*3+2)
        imagesc(squeeze(fMRI_signal(:,round(Y_size/(n_slices+3)*(s+1)),:,t))',[-colorlimitbar/2 colorlimitbar/2])
        axis xy
        axis image
        axis off

        subplot(n_slices,3,(s-1)*3+3)
        imagesc(squeeze(fMRI_signal(round(X_size/(n_slices+3)*(s+1)),:,:,t))',[-colorlimitbar/2 colorlimitbar/2])
        axis xy
        axis image
        axis off
    end
    colormap(opts.select_colormap)  % reapply after imagesc resets it
    drawnow                    % force complete rendering before capture
    if opts.save_video
        frame = print(gcf, '-RGBImage', '-r0');
        [h, w, ~] = size(frame);
        h2 = ceil(h/16)*16;
        w2 = ceil(w/16)*16;
        padded = 255*ones(h2, w2, 3, 'uint8');   % white pad, matches figure background
        padded(1:h, 1:w, :) = frame;
        writeVideo(videoModes, padded);
    else
        pause(0.1)
    end
end

end

if opts.save_video
    close(videoModes);
    disp('    Video saved.')
end

close(gcf)

end
