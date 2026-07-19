function Video_fMRI_any

%%%%%
%
%  VIDEO OF fMRI IN SELECTED BAND
%  Choose the scan and the band
%  Load the fMRI data, filter and generate video.
%
%  Works for Single Slice or Volume
%
%  scripts by Joana Cabral, November 2025
%  joana.barbosa.cabral@tecnico.ulisboa.pt
%
%
%%%%%%%%%%%

% USER: Determine folder and scan name:

general_path='/Users/user/Documents/Research/CSF-MIND/';

tag_name='ep2d_bold_FLIP_18_SAG';
N_dimensions=2; % Or 2D or 3D if volume

% Find the fMRI scan in format .nii with these properties in the NIFTI folder
file_name = dir(fullfile(general_path, ['/**/*' tag_name '*.nii.gz']));

json_name = strrep(file_name.name, '.nii.gz', '.json');
json_path = fullfile(file_name.folder, json_name);


fid = fopen(json_path);
if fid == -1
    warning('Could not find JSON file: %s', json_path);
    TR = input('Please enter manually the TR in seconds: ');
else
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    Scan_info = jsondecode(str);
    TR = Scan_info.RepetitionTime;
end

% USER: Chose to save video (only after debugging)
save_video=1; % Choose 1 to save, otherwise 0
if save_video
    video_acceleration=1; % If 1, video is saved at TR resolution. faster >1, slower <1
    % The folder where the videos will be saved
    Video_folder='/Users/user/Documents/Research/CSF-MIND/Videos/';
end

% USER: Choose to bandpass filter

band_pass=1;
if band_pass
    % Band pass filter the signals into FREQUENCY BANDS
    % [0.005-0.1]; [0.2-0.3]
    high_pass=0.01; % Lowest frequency boundary
    low_pass=0.1; % HIGHEST FREQUENCY BOUNDARY
end
% 
addpath(genpath(general_path))

select_colormap='jet'; % 'jet'; %'bipolar'; % 'redblue'

if band_pass
    % This is to include the frequency band in the video name without dots
    high_pass_label=num2str(high_pass);
    high_pass_label(high_pass_label=='.')='p';
    low_pass_label=num2str(low_pass);
    low_pass_label(low_pass_label=='.')='p';
end


%% Read fMRI file, remove mean of each signal and compute the power spectrum
disp('%%%%% Video fMRI %%%%% ')
disp(['- Now reading file ' file_name.name])
disp(['    TR = ' num2str(TR) ' seconds.'])
if band_pass
    disp(['   Bandpass fitering ' num2str(high_pass) '-' num2str(low_pass) ' Hz'])
else
    disp('    No bandpass fitering applied.')
end

if save_video
    disp(['    Video saved to ' Video_folder])
else
    disp('    Video not saved.')
end

fMRI_signal=niftiread([file_name.folder '/' file_name.name]);

[X_size, Y_size, Z_size, Tmax]=size(fMRI_signal);
fMRI_signal=double(reshape(fMRI_signal,[X_size*Y_size*Z_size, Tmax]));
fMRI_signal=fMRI_signal-mean(fMRI_signal,2); % Remove the mean for filter


if band_pass
    disp(['    Now bandpass filtering ' num2str(high_pass) '-' num2str(low_pass) ' Hz'])
    fMRI_signal = bandpass_fft(fMRI_signal', [high_pass low_pass], 1/TR)';
    % for n=1:size(fMRI_signal,1)
    %     fMRI_signal(n,:)=bandpass(fMRI_signal(n,:),[high_pass low_pass],1/TR);
    % end
    Boundary=round(10/TR); % To cut the first and last 10 seconds of the signals after band pass
    fMRI_signal=fMRI_signal(:,Boundary:end-Boundary);
    Tmax=size(fMRI_signal,2);
end

fMRI_signal=reshape(fMRI_signal,[X_size, Y_size, Z_size, Tmax]);

% Create figure showing A) Mean Signal in each Voxel, STD in 
figure('Position',[1 1 586 884])
colormap(select_colormap)
colorlimitbar = 5 * std(fMRI_signal(:));


%% Generate video of signals
figure('Position',[  1     1   586   884])
colormap(select_colormap)
colorlimitbar=5*std(fMRI_signal(:));

if save_video
    if ~band_pass
    videoModes = VideoWriter([Video_folder '/' file_name.name],'MPEG-4');
    else
    videoModes = VideoWriter([Video_folder '/' file_name.name '_Filt' high_pass_label '-' low_pass_label],'MPEG-4');
    end
    videoModes.FrameRate = round(video_acceleration*1/TR);
    videoModes.Quality = 100;
    open(videoModes);
end

% Run this if just one slice
if N_dimensions==2

for t=1:Tmax

    imagesc(imresize(squeeze(fMRI_signal(:,:,t)),2)',[-colorlimitbar colorlimitbar])
    axis image
    axis off
    axis xy
    title(['T= ' num2str(t*TR,'%2f') ' secs'])
    colormap(select_colormap)  % reapply after imagesc resets it
    drawnow                    % force complete rendering before capture
    if save_video
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

if save_video
    close(videoModes);
    disp('    Video saved.')
end

elseif N_dimensions==3
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
    colormap(select_colormap)  % reapply after imagesc resets it
    drawnow                    % force complete rendering before capture
    if save_video
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

if save_video
    close(videoModes); %close the file
end

end

