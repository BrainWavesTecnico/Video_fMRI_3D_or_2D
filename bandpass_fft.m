function X = bandpass_fft(Y, freqrange, fres)
% Combines bandpass.m + convert_back_to_time.m into one call,
% vectorized across all columns of Y at once.
%   Y         : Nsamples x Nsignals matrix (each column a separate signal)
%   freqrange : [low high] band edges in Hz
%   fres      : sampling frequency (Hz)

[Nsamples, Nsignals] = size(Y);
Nunique_points = ceil((Nsamples+1)/2);
fHz = (0:Nunique_points-1)*fres/Nsamples;
freq_ind = find(fHz >= freqrange(1) & fHz <= freqrange(2));

FY = fft(Y);                        % one FFT call, all voxels at once
full_FY = zeros(Nsamples, Nsignals);
full_FY(freq_ind,:) = FY(freq_ind,:);

tmp = full_FY;
tmp(2:end,:) = tmp(2:end,:) + full_FY(end:-1:2,:);

X = real(ifft(tmp,'symmetric'));    % one IFFT call, all voxels at once
end