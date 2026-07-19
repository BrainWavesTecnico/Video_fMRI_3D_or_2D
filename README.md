# Video_fMRI_3D_or_2D

MATLAB scripts to visualize the temporal evolution of fMRI signals from NIfTI files: a static summary figure (mean, std, power spectrum) and a video of the signal over time, with an optional bandpass filter. Works for both single-slice (2D) and full-volume (3D) scans, and automatically adapts to anisotropic voxels (e.g. thick-slice 2D EPI without slice-timing correction).

## Requirements

- MATLAB with the Image Processing Toolbox (`niftiread`, `niftiinfo`, `imresize`) and the Signal Processing Toolbox is **not** required - bandpass filtering is done with a custom FFT-based implementation (`bandpass_fft.m`), no toolbox dependency.
- Each `.nii.gz` file should have a matching BIDS-style `.json` sidecar (as produced by `dcm2niix`) in the same folder, with a `RepetitionTime` field. If the JSON is missing, the script prompts for the TR manually instead of failing.

## Files

| File | Role |
|---|---|
| `Main_fMRI_Video.m` | Entry point. Configure paths/options here and run this. Finds all matching scans, plots+saves a static summary figure for each, then calls `Video_fMRI_any` to generate its video. |
| `Video_fMRI_any.m` | Bandpass-filters (optional) a single scan and writes its video. Called once per file by `Main_fMRI_Video`. |
| `bandpass_fft.m` | Vectorized FFT-based bandpass filter, applied to all voxels at once. |
| `get_plane_image.m` | Shared helper that extracts an oriented 2D slice (sagittal/coronal/axial) from a 3D volume; used by both the static figure and the video so their slicing/orientation stay consistent. |

## Usage

Open `Main_fMRI_Video.m` and edit the settings under **section 1** (marked `USER:`), then run the function.

### Key options

- `general_path` / `tag_name` - folder to search (recursively) and a substring to match in the scan's filename. Every `*.nii.gz` file matching `tag_name` is processed in one run.
- `opts.save_video` - if `0`, figures are shown live instead of being written to a video file (useful for debugging without waiting for a full render).
- `opts.Figures_and_Videos_folder` - single output folder for both the static `.jpeg` figures and the `.mp4` videos. Created automatically if it doesn't exist.
- `opts.band_pass`, `opts.high_pass`, `opts.low_pass` - whether to bandpass filter, and the frequency band (Hz). The first/last 10 seconds are trimmed after filtering to remove edge artifacts.
- `opts.select_colormap` - colormap used for the video (default `'jet'`).
- `opts.max_voxel_anisotropy` - ratio threshold (default `1.5`) between the largest and smallest voxel dimension, above which a volume is treated as anisotropic (see below).
- `opts.n_slices_single_plane` - how many equidistant slices to show for an anisotropic volume (default `6`).

## Pipeline

For each matching file:

1. **Read TR and header.** TR comes from the JSON sidecar; voxel size and image dimensions come from the NIfTI header (`niftiinfo`), read before the (slower) full data load so this info can be logged immediately: file name, TR, scan type, voxel size, and filter settings.
2. **Static figure** (`<file>_Mean_STD_PSD.jpeg`): mean signal, temporal std, and power spectrum (after linear detrending) across voxels. Layout depends on the scan:
   - **Single slice:** 2×2 grid (mean, std, power spectrum spanning the bottom row).
   - **Volume, isotropic voxels:** 3×3 grid - mean and std for the middle Sagittal/Axial/Coronal slice, power spectrum below.
   - **Volume, anisotropic voxels:** 3×N grid (`N = opts.n_slices_single_plane`) - mean and std for N equidistant slices of *only* the plane that was actually acquired (see below), power spectrum below.
3. **Video** (`<file>.mp4`, or `<file>_Filt<low>-<high>.mp4` if bandpassed): the signal is linearly detrended, optionally bandpass-filtered (`bandpass_fft`), then rendered frame by frame. Each frame's title shows the filter band (or "No filter") and the elapsed time. Layout mirrors the static figure's logic (single slice / isotropic 3-plane / anisotropic N-slice single-plane).

### Single slice vs. volume

fMRI data is always stored as 4D (X, Y, Z, T). A "single slice" scan is simply one where one of X, Y, Z equals 1 - which dimension depends on acquisition orientation (sagittal/coronal/axial), so this is detected from the data shape rather than a manual flag.

### Isotropic vs. anisotropic volumes

A typical multi-slice 2D EPI acquisition (no slice-timing correction) has much thicker voxels through the slice direction than in-plane, and slices within a volume are acquired at different times. Reformatting such data into the other two anatomical planes would mix data from different acquisition times and require interpolating across large gaps - so when `max(voxel size)/min(voxel size)` exceeds `opts.max_voxel_anisotropy`, only the plane that was actually acquired is shown (whichever of Sagittal/Axial/Coronal corresponds to the thick voxel dimension - not hardcoded to axial). Isotropic volumes show all three planes as usual.
