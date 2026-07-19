function Video_fMRI_any(fMRI_signal, TR, file_label, opts, is_volume, planes)

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
%   is_volume   : true if fMRI_signal is a full volume, false if a single
%                 slice (one of X, Y, Z equals 1) - computed once by the
%                 caller from the data shape.
%   planes      : struct array (name, dim) of the planes to show for a
%                 volume - one plane if voxels are anisotropic, three
%                 (Sagittal/Axial/Coronal) if isotropic. Unused if
%                 is_volume is false. See get_plane_image.m.
%
%  scripts by Joana Cabral, July 2026
%  joanabcabral@tecnico.ulisboa.pt
%
%%%%%%%%%%%

[X_size, Y_size, Z_size, Tmax]=size(fMRI_signal);
fMRI_signal=double(reshape(fMRI_signal,[X_size*Y_size*Z_size, Tmax]));
fMRI_signal=detrend(fMRI_signal')'; % linear detrend along time (removes trend and mean) before filtering

if opts.band_pass
    disp(['    Now bandpass filtering ' num2str(opts.high_pass) '-' num2str(opts.low_pass) ' Hz'])
    fMRI_signal = bandpass_fft(fMRI_signal', [opts.high_pass opts.low_pass], 1/TR)';
    Boundary=round(10/TR); % To cut the first and last 10 seconds of the signals after band pass
    fMRI_signal=fMRI_signal(:,Boundary:end-Boundary);
    Tmax=size(fMRI_signal,2);
end

fMRI_signal=reshape(fMRI_signal,[X_size, Y_size, Z_size, Tmax]);

%% Generate video of signals
if ~is_volume
    figure('Position',[158 103 1133 884])
else
    figure('units','normalized','outerposition',[0 0 1 0.5]) % keep width, half height
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
    title(['T= ' num2str(t*TR,'%.2f') ' secs'],'FontSize',22)
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

elseif numel(planes) == 1
% Anisotropic volume: only the acquired plane is meaningful. Show
% several equidistant slices of it in a single row.

n_show = opts.n_slices_single_plane;
dim = planes(1).dim;
sizes3=[X_size Y_size Z_size];

for t=1:Tmax
    vol_t = squeeze(fMRI_signal(:,:,:,t));
    for s=1:n_show
        idx = round(sizes3(dim)/(n_show+3)*(s+1));

        subplot(1,n_show,s)
        imagesc(get_plane_image(vol_t, dim, idx),[-colorlimitbar/2 colorlimitbar/2])
        if dim ~= 3, axis xy; end
        axis image
        axis off
        if s==1
            title(['T= ' num2str(t*TR,'%.2f') ' secs'],'FontSize',22)
        end
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

else
% Isotropic volume: n_slices depth samples (rows) for each of the 3
% planes (columns) - unchanged from before.

n_slices=3;
n_planes=numel(planes);
sizes3=[X_size Y_size Z_size];

for t=1:Tmax
    vol_t = squeeze(fMRI_signal(:,:,:,t));
    for s=1:n_slices
        for p=1:n_planes
            dim = planes(p).dim;
            idx = round(sizes3(dim)/(n_slices+3)*(s+1));

            subplot(n_slices,n_planes,(s-1)*n_planes+p)
            imagesc(get_plane_image(vol_t, dim, idx),[-colorlimitbar/2 colorlimitbar/2])
            if dim ~= 3, axis xy; end
            axis image
            axis off
            if s==1 && p==1
                title(['T= ' num2str(t*TR,'%.2f') ' secs'],'FontSize',22)
            end
        end
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
