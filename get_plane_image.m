function img = get_plane_image(vol, dim, idx)
% Extract a 2D image from a 3D volume (X x Y x Z) by fixing array
% dimension `dim` at index `idx`, oriented to match the sagittal (dim 1)
% / coronal (dim 2) / axial (dim 3) display convention used throughout
% this repo: sagittal and coronal are transposed, axial is not. Callers
% should apply "axis xy" after imagesc when dim ~= 3, to match.

switch dim
    case 1
        img = squeeze(vol(idx,:,:))';
    case 2
        img = squeeze(vol(:,idx,:))';
    case 3
        img = squeeze(vol(:,:,idx));
end
end
