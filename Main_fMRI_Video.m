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

% USER: threshold for treating a volume's voxels as isotropic. If the
% ratio between the largest and smallest voxel dimension exceeds this,
% the scan is treated as anisotropic (thick slices, as in a typical
% single-band 2D EPI acquisition) and only the plane that was actually
% acquired is shown - the other two views would mix data recorded at
% different times (no slice-timing correction) and need interpolation
% across large gaps between slices.
opts.max_voxel_anisotropy = 1.5;

% USER: how many equidistant slices to show for an anisotropic volume
% (single acquired plane), both in the static figure and the video.
opts.n_slices_single_plane = 6;

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

% Read the header first (fast - no need to load the full data yet) to
% report the scan's shape and voxel size, and to decide the layout of
% the figures below.
nii_info = niftiinfo([file_name.folder '/' file_name.name]);
img_size = nii_info.ImageSize;
X_size = img_size(1); Y_size = img_size(2); Z_size = img_size(3); Tmax = img_size(4);
voxel_size = nii_info.PixelDimensions(1:3);

% fMRI data is always 4D (X,Y,Z,T). A single-slice ("2D") scan is simply
% one where one of the first three dimensions equals 1 - could be X, Y or
% Z depending on the slice orientation (SAG/COR/AX), so we detect it
% directly from the header instead of a manual flag.
is_volume = ~any([X_size Y_size Z_size] == 1);
if is_volume
    type_label = 'Volume (3D)';
else
    type_label = 'Single slice (2D)';
end

% Sagittal/Axial/Coronal plane definitions: which array dimension is
% held fixed for each view (see get_plane_image.m for orientation).
plane_defs = struct('name',{'Sagittal','Axial','Coronal'},'dim',{1,3,2});

if is_volume
    [~, thick_dim] = max(voxel_size);
    is_isotropic = (max(voxel_size)/min(voxel_size)) < opts.max_voxel_anisotropy;
    if is_isotropic
        planes = plane_defs;
    else
        planes = plane_defs([plane_defs.dim]==thick_dim);
    end
else
    planes = plane_defs([]); % unused for a single slice
end

disp(['- Now reading file ' file_name.name ' (' num2str(f) '/' num2str(numel(file_list)) ')'])
disp(['    TR = ' num2str(TR,'%.3f') ' seconds.'])
disp(['    Type: ' type_label])
disp(['    Voxel size: ' num2str(voxel_size,'%.2f  ') ' mm'])
if is_volume && ~is_isotropic
    disp(['    Anisotropic voxels - showing only the ' planes(1).name ' plane.'])
end
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

if ~is_volume
    fig_static = figure('Position',[ 428   159   978   831]);
else
    fig_static = figure('units','normalized','outerposition',[0 0 1 1]); % full screen
end
colormap(jet)

mean_img = squeeze(mean(fMRI_signal,4));   % average over time
std_img  = squeeze(std(fMRI_signal,[],4)); % std over time

if ~is_volume
    subplot(2,2,1)
    imagesc(imresize(mean_img',2))
        axis image
        axis off
        axis xy
        title('Mean signal in each voxel')

    subplot(2,2,2)
    image_to_plot=imresize(std_img',2);
    imagesc(image_to_plot,[0 5*std(image_to_plot(:))])
        axis image
        axis off
        axis xy
        title('Signal variance in each voxel')

    psd_subplot = {2,2,3:4};
elseif numel(planes) == 1
    % Anisotropic: show several equidistant slices of the one acquired plane
    dim = planes(1).dim;
    n_show = opts.n_slices_single_plane;
    sizes3 = [X_size Y_size Z_size];

    for s = 1:n_show
        idx = round(sizes3(dim)/(n_show+3)*(s+1));

        subplot(3,n_show,s)
        imagesc(get_plane_image(mean_img, dim, idx))
        if dim ~= 3, axis xy; end
        axis image; axis off
        if s==1
            title({'Mean signal in each voxel', [planes(1).name ' ' num2str(s)]})
        else
            title([planes(1).name ' ' num2str(s)])
        end

        subplot(3,n_show,n_show+s)
        img = get_plane_image(std_img, dim, idx);
        imagesc(img,[0 5*std(img(:))])
        if dim ~= 3, axis xy; end
        axis image; axis off
        if s==1
            title({'STD in each voxel', ['STD ' num2str(s)]})
        else
            title(['STD ' num2str(s)])
        end
    end

    psd_subplot = {3,n_show,(2*n_show+1):(3*n_show)};
else
    % Isotropic: middle slice of each of the 3 planes
    n_planes = numel(planes);
    mid_idx = round([X_size Y_size Z_size]/2);

    % Row 1: mean, Row 2: std, one column per plane shown.
    % Orientation (transpose / axis xy) matches the volume video below.
    for p = 1:n_planes
        dim = planes(p).dim;
        idx = mid_idx(dim);

        subplot(3,n_planes,p)
        imagesc(get_plane_image(mean_img, dim, idx))
        if dim ~= 3, axis xy; end
        axis image; axis off
        title(['Mean - ' planes(p).name])

        subplot(3,n_planes,n_planes+p)
        img = get_plane_image(std_img, dim, idx);
        imagesc(img,[0 5*std(img(:))])
        if dim ~= 3, axis xy; end
        axis image; axis off
        title(['STD - ' planes(p).name])
    end

    psd_subplot = {3,n_planes,(2*n_planes+1):(3*n_planes)};
end

fMRI_signal_detrended=double(reshape(fMRI_signal,[X_size*Y_size*Z_size, Tmax]));
fMRI_signal_detrended=fMRI_signal_detrended-mean(fMRI_signal_detrended,2);

N = size(fMRI_signal_detrended, 2);      % number of time points
fs = 1/TR;                               % sampling frequency in Hz
freq = (0:N-1) * (fs/N);                 % frequency for each FFT bin, 0 to fs
Nhalf = floor(N/2) + 1;
freq_half = freq(1:Nhalf);

subplot(psd_subplot{:})
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

Video_fMRI_any(fMRI_signal, TR, file_name.name, opts, is_volume, planes);

end

end
