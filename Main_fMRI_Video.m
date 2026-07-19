function Main_fMRI_Video

%%%%%
%
%  MAIN PIPELINE: fMRI VIDEO GENERATION
%
%  1) Define folder, files, TR and options
%  2) Plot and save mean, std and power spectrum of the raw signal
%  3) Bandpass filter (optional) and generate the video
%
%  Loops over all files matching tag_name in the NIFTI folder.
%
%  scripts by Joana Cabral, July 2026
%  joanabcabral@tecnico.ulisboa.pt
%
%%%%%%%%%%%

%% 1) USER: Define folder, files, TR and options

general_path='/Users/user/Documents/Research/CSF-MIND/';
tag_name='ep2d_bold_FLIP_18_SAG';

addpath(genpath(general_path))

% USER: Choose to save video (only after debugging)
opts.save_video=1; % Choose 1 to save, otherwise 0
if opts.save_video
    opts.video_acceleration=1; % If 1, video is saved at TR resolution. faster >1, slower <1
end

% The folder where the figures and videos will be saved
opts.Figures_and_Videos_folder='/Users/user/Documents/Research/CSF-MIND/Figures_and_Videos/';
if ~exist(opts.Figures_and_Videos_folder,'dir')
    mkdir(opts.Figures_and_Videos_folder)
end

% USER: Choose to bandpass filter
opts.band_pass=1;
if opts.band_pass
    % Band pass filter the signals into FREQUENCY BANDS
    % [0.005-0.1]; [0.2-0.3]
    opts.high_pass=0.02; % Lowest frequency boundary
    opts.low_pass=0.12;  % HIGHEST FREQUENCY BOUNDARY
end

opts.select_colormap='jet'; % 'jet'; %'bipolar'; % 'redblue'

% Find all fMRI scans in format .nii with these properties in the NIFTI folder
file_list = dir(fullfile(general_path, ['/**/*' tag_name '*.nii.gz']));

disp('%%%%% Video fMRI %%%%% ')
disp(['Found ' num2str(numel(file_list)) ' file(s) matching "' tag_name '"'])

for f = 1:numel(file_list)

file_name = file_list(f);

json_name = strrep(file_name.name, '.nii.gz', '.json');
json_path = fullfile(file_name.folder, json_name);

fid = fopen(json_path);
if fid == -1
    warning('Could not find JSON file: %s', json_path);
    TR = input(['Please enter manually the TR in seconds for ' file_name.name ': ']);
else
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    Scan_info = jsondecode(str);
    TR = Scan_info.RepetitionTime;
end

disp(['- Now reading file ' file_name.name ' (' num2str(f) '/' num2str(numel(file_list)) ')'])
disp(['    TR = ' num2str(TR) ' seconds.'])
if opts.band_pass
    disp(['    Bandpass filter ' num2str(opts.high_pass) '-' num2str(opts.low_pass) ' Hz'])
else
    disp('    No bandpass filtering applied.')
end
if opts.save_video
    disp(['    Video saved to ' opts.Figures_and_Videos_folder])
else
    disp('    Video not saved.')
end

fMRI_signal=single(niftiread([file_name.folder '/' file_name.name]));

%% 2) Static figure: Mean signal, signal variance and power spectrum

[X_size, Y_size, Z_size, Tmax]=size(fMRI_signal);

fig_static = figure('Position',[ 428   159   978   831]);
colormap(jet)

subplot(2,2,1)
imagesc(imresize(mean(squeeze(fMRI_signal),3)',2))
    axis image
    axis off
    axis xy
    title('Mean signal in each voxel')

subplot(2,2,2)
image_to_plot=imresize(std(squeeze(fMRI_signal),[],3),2)';
imagesc(image_to_plot,[0 5*std(image_to_plot(:))])
    axis image
    axis off
    axis xy
    title('Signal variance in each voxel')

fMRI_signal_detrended=double(reshape(fMRI_signal,[X_size*Y_size*Z_size, Tmax]));
fMRI_signal_detrended=fMRI_signal_detrended-mean(fMRI_signal_detrended,2);

N = size(fMRI_signal_detrended, 2);      % number of time points
fs = 1/TR;                               % sampling frequency in Hz
freq = (0:N-1) * (fs/N);                 % frequency for each FFT bin, 0 to fs
Nhalf = floor(N/2) + 1;
freq_half = freq(1:Nhalf);

subplot(2,2,3:4)
powerSpectrum = abs(fft(fMRI_signal_detrended, [], 2)).^2;   % FFT along time (dim 2)
plot(freq_half, powerSpectrum(:,1:Nhalf)');
Nyquist = fs/2;                 % = 1/(2*TR)
xlim([0 Nyquist/2]);
set(gca,'YScale', 'linear' );
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum from each voxel');

static_fig_name = strrep(file_name.name, '.nii.gz', '_Mean_STD_PSD.jpeg');
exportgraphics(fig_static, fullfile(opts.Figures_and_Videos_folder, static_fig_name), 'Resolution', 150);
close(fig_static)

%% 3) Bandpass filter (optional) and generate the video

Video_fMRI_any(fMRI_signal, TR, file_name.name, opts);

end

end
