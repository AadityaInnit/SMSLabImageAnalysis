%% ============================================================
%  spectral_kymograph.m
%  Chowdhury Lab — Single-Molecule Spectroscopy, IIT Bombay
%  ============================================================
%
% ── MICROSCOPE LAYOUT ───────────────────────────────────────────
%
%  A transmission grating in the emission path splits each CCD
%  frame into two side-by-side images:
%
%  ┌──────────────────────────────────────────────────────────┐
%  │  ZEROTH ORDER (left) │  FIRST ORDER / spectra (right)   │
%  │  bright vertical     │  one horizontal streak per QD    │
%  │  line — each QD is   │  x-axis = wavelength             │
%  │  y-axis = spatial    │  y-axis = same spatial row       │
%  └──────────────────────────────────────────────────────────┘
%
% ── PIPELINE (5 steps) ──────────────────────────────────────────
%
%  1. LOAD       Read all TIFF frames into a 3-D array.
%
%  2. LOCATE     Define the zeroth-order search range as the
%                user-specified column interval
%                [ZORDER_COL_MIN, ZORDER_COL_MAX].  All particles
%                whose zeroth-order spots fall within this band are
%                eligible for detection.  The first-order region
%                begins at ZORDER_COL_MAX + ZORDER_GAP so that
%                every particle's spectrum is cleanly separated
%                from its zeroth-order spot regardless of where
%                within the band it sits.
%
%  3. DETECT     Build a 1-D spatial profile by taking the MAX
%                intensity across ALL columns in the zeroth-order
%                range, per row, using the MAXIMUM-INTENSITY
%                PROJECTION (MIP) across frames.  This captures
%                every QD at its peak brightness regardless of
%                blinking duty cycle and regardless of which
%                column within the range it occupies.
%
%                Use findpeaks() for NMS.  For each surviving peak,
%                find the particle's individual x coordinate (the
%                brightest column in the range at that row in the
%                MIP), then confirm with a 1-D Gaussian shape test
%                on the 5x5 block centred on that (row, x) position.
%
%                When two confirmed peaks are within
%                MIN_SEPARATION rows, keep the brighter one.
%
%  4. EXTRACT    For each particle at row R, fit a Gaussian to its
%                average first-order spectrum to find mu (peak
%                pixel) and sigma (spectral width).  Extract the
%                kymograph over mu +/- 40*sigma columns, widening
%                the row band by +/-VERTICAL_ALLOWANCE to tolerate
%                optical distortion between the two orders.
%                Apply column-wise percentile background
%                subtraction.  Convert pixel axis to nm using
%                the lab's calibration:  lambda = px * 1.2 - 40.
%
%  5. SAVE       For each particle: kymograph figure (raw, BG-
%                corrected, average spectrum, blinking trace).
%                Runs multi-Gaussian decomposition to track ZPL and
%                phonon sidebands over time.
%                A summary table (.txt) is written with the
%                per-particle zeroth-order (x,y) coordinate.
%
% ── KEY DESIGN DECISIONS ────────────────────────────────────────
%
%  • Detection uses the MIP (not mean_frame) for the spatial
%    profile and Gaussian block, so blinking does not hide bright
%    particles.
%  • The zeroth-order search is a user-defined column RANGE rather
%    than a single line, allowing for wide-field illumination where
%    particles are spread across multiple x columns.
%  • Each detected particle records its own x coordinate (the
%    brightest column in the range at its row in the MIP).
%  • spec_start is anchored to ZORDER_COL_MAX so the first-order
%    region is consistent for all particles.
%  • mean_frame is used for all downstream extraction and display.
%
% ── DEPENDENCIES ────────────────────────────────────────────────
%  Base MATLAB only — no toolboxes required.
%  (findpeaks is in Signal Processing Toolbox; a fallback is
%   included for systems without it.)
%
% ============================================================

clear; clc; close all;

%% ── USER PARAMETERS ─────────────────────────────────────────────────────────

% Full path to your TIFF movie stack.
TIFF_PATH = '/Users/aaditya/Desktop/SMS Lab/IMP_473_CsPbBrI2_100ms_8.4watt-cm2_0OD_1_MMStack_Default.ome-1.tif';

% ── Zeroth-order detection range ─────────────────────────────────────────────
%  All columns between ZORDER_COL_MIN and ZORDER_COL_MAX (inclusive) are
%  searched for zeroth-order particle spots.  Any particle whose PSF
%  peak falls within this x range will be detected.
%
%  For IMP_473 / wide-field geometry:
%    ZORDER_COL_MIN = 13   (leftmost column to search)
%    ZORDER_COL_MAX = 64   (rightmost column to search)
%
%  The first-order (spectral) region automatically starts at
%  ZORDER_COL_MAX + ZORDER_GAP, keeping the spectrum window
%  consistent for every particle regardless of its x position.

ZORDER_COL_MIN = 13;   % px — left bound of zeroth-order search range
ZORDER_COL_MAX = 64;   % px — right bound of zeroth-order search range

% Detection sensitivity.
% A row is a candidate when its max intensity exceeds:
%     median(background) + THRESHOLD * std(background)
% Recommended starting range: 5–10.
% Raise to reduce false positives; lower to catch very dim emitters.
THRESHOLD = 8;

% Minimum row-distance (pixels) between two detected particles.
% findpeaks uses this as MinPeakDistance.
% Set to at least the expected PSF diameter (~3–5 px for 60x/1.49 NA).
MIN_SEPARATION = 5;

% Spatial half-width for first-order kymograph row integration (px).
HALF_WIDTH = 2;

% Extra vertical pixels added to HALF_WIDTH on both sides when reading
% the first-order region, to tolerate vertical optical distortion.
VERTICAL_ALLOWANCE = 2;

% Frame acquisition rate in Hz.  Set [] to label axes in frames.
FRAME_RATE = 20;

% Colormap for kymograph images.
COLORMAP = 'hot';

% Percentile used for column-wise background subtraction of kymographs.
BG_PERCENTILE = 5;

% Spectral stability threshold (nm).
STATIC_THRESHOLD_NM = 5.0;

% Cap on number of particles to process ([] = process all).
MAX_PARTICLES = [];

% Save kymograph arrays to a .mat file alongside the figures?
SAVE_MAT = true;

%% ── WAVELENGTH CALIBRATION ───────────────────────────────────────────────────
%  lambda (nm) = pixel_index * 1.2 - 40
%  pixel_index 1 = first column of the first-order region (spec_start).

px2nm = @(px) px .* 1.2 - 40;

%% ══════════════════════════════════════════════════════════════════════════════
%  STEP 1 — LOAD
%% ══════════════════════════════════════════════════════════════════════════════

fprintf('\n%s\n Chowdhury Lab — Spectral Kymograph Generator\n%s\n\n', ...
        repmat('=',1,60), repmat('=',1,60));

if ~isfile(TIFF_PATH)
    error('File not found: %s\nCheck TIFF_PATH in the USER PARAMETERS section.', ...
          TIFF_PATH);
end

info     = imfinfo(TIFF_PATH);
n_frames = numel(info);
height   = info(1).Height;
width    = info(1).Width;

fprintf('[1/5] Loading: %s\n', TIFF_PATH);
fprintf('      %d frames  x  %d rows  x  %d cols\n', n_frames, height, width);

% Read all frames — single precision halves memory vs double.
movie = zeros(height, width, n_frames, 'single');
tobj  = Tiff(TIFF_PATH, 'r');
for f = 1:n_frames
    tobj.setDirectory(f);
    movie(:,:,f) = single(tobj.read());
end
tobj.close();

% Time-averaged frame — used for all downstream steps (Steps 4 & 5).
mean_frame = mean(movie, 3);

% Maximum-intensity projection — used ONLY for particle detection (Step 3).
% Captures every QD at its all-time peak brightness, preventing
% dim-but-persistent emitters from masking bright-but-blinking ones.
max_frame  = max(movie, [], 3);

%% ══════════════════════════════════════════════════════════════════════════════
%  STEP 2 — DEFINE ZEROTH-ORDER RANGE AND FIRST-ORDER REGION
%% ══════════════════════════════════════════════════════════════════════════════

fprintf('\n[2/5] Zeroth-order range and first-order region\n');

% Validate range against actual frame dimensions.
if ZORDER_COL_MIN < 1 || ZORDER_COL_MAX > width || ZORDER_COL_MIN >= ZORDER_COL_MAX
    error(['ZORDER_COL_MIN / ZORDER_COL_MAX are out of bounds or inverted.\n' ...
           'Frame width = %d px.  Check USER PARAMETERS.'], width);
end

% The first-order region starts ZORDER_GAP columns past the right edge of
% the zeroth-order range.  Anchoring to ZORDER_COL_MAX (rather than to
% any individual particle's x) keeps spec_start uniform for all particles.
ZORDER_GAP  = 15;
spec_start  = ZORDER_COL_MAX + ZORDER_GAP;
spec_end    = width;
n_spec_cols = spec_end - spec_start + 1;

if n_spec_cols < 10
    error(['spec_start (%d) leaves fewer than 10 spectral columns.\n' ...
           'Decrease ZORDER_COL_MAX or ZORDER_GAP, or check frame width (%d).'], ...
          spec_start, width);
end

fprintf('      Zeroth-order range  : cols %d – %d  (%d px wide)\n', ...
        ZORDER_COL_MIN, ZORDER_COL_MAX, ZORDER_COL_MAX - ZORDER_COL_MIN + 1);
fprintf('      First-order region  : cols %d – %d  (%d px)\n', ...
        spec_start, spec_end, n_spec_cols);

%% ══════════════════════════════════════════════════════════════════════════════
%  STEP 3 — PARTICLE DETECTION
%% ══════════════════════════════════════════════════════════════════════════════

fprintf('\n[3/5] Detecting particles\n');

% ── 3a. Build 1-D spatial profile across the full zeroth-order range ─────────
%
%  For each row take the maximum intensity across ALL columns in
%  [ZORDER_COL_MIN, ZORDER_COL_MAX] in the MIP.  This makes the
%  profile sensitive to every particle in the range regardless of
%  which x column it occupies, and uses the MIP so blinking QDs
%  are not suppressed by time-averaging.

spatial_profile = max( ...
    max_frame(:, ZORDER_COL_MIN:ZORDER_COL_MAX), [], 2);   % [height x 1]

% ── 3b. Robust noise threshold ───────────────────────────────────────────────
sorted_prof = sort(spatial_profile);
bg_vals     = sorted_prof(1 : floor(end/2));
bg_median   = median(double(bg_vals));
bg_std      = std(double(bg_vals));
bg_std      = max(bg_std, 1e-6);
thresh_val  = bg_median + THRESHOLD * bg_std;

fprintf('      Background:  median = %.1f,  std = %.2f\n', bg_median, bg_std);
fprintf('      Threshold  = %.1f counts  (%.0f-sigma above bg)\n', ...
        thresh_val, THRESHOLD);

% ── 3c. Find local maxima with findpeaks (built-in NMS) ──────────────────────
try
    [peak_vals, peak_rows] = findpeaks(double(spatial_profile), ...
        'MinPeakHeight',   thresh_val, ...
        'MinPeakDistance', MIN_SEPARATION);
catch
    % Fallback — no Signal Processing Toolbox.
    peak_rows = []; peak_vals = [];
    above = find(spatial_profile > thresh_val);
    if ~isempty(above)
        vals_above = spatial_profile(above);
        [~, sort_ord]  = sort(vals_above, 'descend');
        sorted_idx = above(sort_ord);
        kept = true(size(sorted_idx));
        for ii = 1:numel(sorted_idx)
            if ~kept(ii), continue; end
            for jj = ii+1:numel(sorted_idx)
                if abs(sorted_idx(ii) - sorted_idx(jj)) < MIN_SEPARATION
                    kept(jj) = false;
                end
            end
        end
        peak_rows = sorted_idx(kept);
        peak_vals = spatial_profile(peak_rows);
    end
end

fprintf('      %d local maxima above threshold\n', numel(peak_rows));

% ── 3d. Per-particle x location + Gaussian shape confirmation ────────────────
%
%  For each candidate row r:
%    1. Find the particle's x coordinate — the brightest column within
%       [ZORDER_COL_MIN, ZORDER_COL_MAX] at row r in the MIP.
%       A +/-WIN row neighbourhood is summed first to reduce noise.
%    2. Build a (2*WIN+1) x (2*WIN+1) block centred on (r, px_col) in the MIP.
%    3. Fit a 1-D Gaussian along the spatial (row) direction.
%    4. Accept if: |mu| <= WIN  AND  SIGMA_MIN <= sigma <= SIGMA_MAX.
%
%  The amplitude is NOT thresholded here — that was enforced by findpeaks.
%  This step tests PSF shape only.

WIN       = 2;      % half-width for block and row neighbourhood (px)
SIGMA_MIN = 0.5;    % px
SIGMA_MAX = 4.0;    % px

confirmed_rows = [];
confirmed_cols = [];   % per-particle zeroth-order x coordinate
confirmed_vals = [];

for k = 1:numel(peak_rows)
    r = peak_rows(k);

    % ── Find per-particle x ──────────────────────────────────────────
    %  Sum a +/-WIN row neighbourhood in the MIP across the range to
    %  reduce single-row noise spikes, then locate the peak column.
    r0_nb    = max(1,      r - WIN);
    r1_nb    = min(height, r + WIN);
    col_sums = sum(max_frame(r0_nb:r1_nb, ZORDER_COL_MIN:ZORDER_COL_MAX), 1);
    [~, best_rel] = max(col_sums);
    px_col   = ZORDER_COL_MIN + best_rel - 1;   % absolute column index

    % ── Gaussian shape test centred on (r, px_col) ───────────────────
    r0 = max(1,      r      - WIN);
    r1 = min(height, r      + WIN);
    c0 = max(1,      px_col - WIN);
    c1 = min(width,  px_col + WIN);

    block = max_frame(r0:r1, c0:c1);
    sp5   = sum(block, 2);
    x_rel = double((r0:r1)' - r);

    [~, mu_g, sigma_g, ~] = fit_gaussian1d(x_rel, sp5);

    if isnan(mu_g)
        fprintf('      Row %4d (x=%3d): Gaussian fit failed — rejected\n', r, px_col);
        continue
    end
    if abs(mu_g) > WIN
        fprintf('      Row %4d (x=%3d): centre off-window (mu=%.2f) — rejected\n', ...
                r, px_col, mu_g);
        continue
    end
    if sigma_g < SIGMA_MIN || sigma_g > SIGMA_MAX
        fprintf('      Row %4d (x=%3d): sigma=%.2f out of [%.1f,%.1f] — rejected\n', ...
                r, px_col, sigma_g, SIGMA_MIN, SIGMA_MAX);
        continue
    end

    confirmed_rows(end+1) = r;            %#ok<AGROW>
    confirmed_cols(end+1) = px_col;       %#ok<AGROW>
    confirmed_vals(end+1) = peak_vals(k); %#ok<AGROW>
    fprintf('      Row %4d (x=%3d): CONFIRMED  (mu=%.2f px, sigma=%.2f px)\n', ...
            r, px_col, mu_g, sigma_g);
end

fprintf('      %d particle(s) confirmed after shape test\n', numel(confirmed_rows));

% ── 3e. Keep-brighter rule for adjacent confirmed particles ───────────────────
%
%  When two confirmed peaks are within MIN_SEPARATION rows, discard the
%  dimmer one.  This replaces the original "reject-both" rule that caused
%  the correct bright particle to be lost whenever a faint neighbour
%  passed the shape test.

keep2 = true(size(confirmed_rows));
for ii = 1:numel(confirmed_rows)
    for jj = ii+1:numel(confirmed_rows)
        if abs(confirmed_rows(ii) - confirmed_rows(jj)) <= MIN_SEPARATION
            if confirmed_vals(ii) >= confirmed_vals(jj)
                keep2(jj) = false;
                fprintf(['      Rows %d & %d within %d px: ' ...
                         'kept row %d (x=%d, %.0f cnts), ' ...
                         'dropped row %d (x=%d, %.0f cnts)\n'], ...
                    confirmed_rows(ii), confirmed_rows(jj), MIN_SEPARATION, ...
                    confirmed_rows(ii), confirmed_cols(ii), confirmed_vals(ii), ...
                    confirmed_rows(jj), confirmed_cols(jj), confirmed_vals(jj));
            else
                keep2(ii) = false;
                fprintf(['      Rows %d & %d within %d px: ' ...
                         'kept row %d (x=%d, %.0f cnts), ' ...
                         'dropped row %d (x=%d, %.0f cnts)\n'], ...
                    confirmed_rows(ii), confirmed_rows(jj), MIN_SEPARATION, ...
                    confirmed_rows(jj), confirmed_cols(jj), confirmed_vals(jj), ...
                    confirmed_rows(ii), confirmed_cols(ii), confirmed_vals(ii));
            end
        end
    end
end

n_adj_dropped = sum(~keep2);
if n_adj_dropped > 0
    fprintf('      %d particle(s) removed (dimmer of adjacent pair)\n', n_adj_dropped);
end

confirmed_rows = confirmed_rows(keep2);
confirmed_cols = confirmed_cols(keep2);
confirmed_vals = confirmed_vals(keep2);

% Apply particle cap (keep brightest)
if ~isempty(MAX_PARTICLES) && numel(confirmed_rows) > MAX_PARTICLES
    [~, bright_ord] = sort(confirmed_vals, 'descend');
    confirmed_rows  = confirmed_rows(bright_ord(1:MAX_PARTICLES));
    confirmed_cols  = confirmed_cols(bright_ord(1:MAX_PARTICLES));
    confirmed_vals  = confirmed_vals(bright_ord(1:MAX_PARTICLES));
    fprintf('      Capped at %d brightest\n', MAX_PARTICLES);
end

% Sort top-to-bottom for consistent numbering
[confirmed_rows, row_ord] = sort(confirmed_rows, 'ascend');
confirmed_cols = confirmed_cols(row_ord);
confirmed_vals = confirmed_vals(row_ord);

n_particles = numel(confirmed_rows);

% Pack into matrix:
%   particles(:,1) = row (y)
%   particles(:,2) = amplitude (profile max from MIP)
%   particles(:,3) = col (x) — individual zeroth-order column per particle
particles = [confirmed_rows(:), confirmed_vals(:), confirmed_cols(:)];

fprintf('\n      --> %d particle(s) accepted\n', n_particles);
fprintf('          %-10s  %-10s  %-10s  %-14s\n', ...
        'Particle', 'x (col)', 'y (row)', 'Profile max');
fprintf('          %s\n', repmat('-',1,46));
for i = 1:n_particles
    fprintf('          %-10d  %-10d  %-10d  %-14.1f\n', ...
            i, particles(i,3), particles(i,1), particles(i,2));
end
fprintf('          %s\n', repmat('-',1,46));
fprintf('          x is the per-particle brightest column in the ZO range\n');
fprintf('          y is the detected row in the original TIFF\n');

if n_particles == 0
    warning('No particles detected. Try lowering THRESHOLD or widening the zeroth-order range.');
    return
end

%% ── DETECTION OVERVIEW FIGURE ────────────────────────────────────────────────

[fdir, fstem, ~] = fileparts(TIFF_PATH);
out_dir = fullfile(fdir, [fstem '_kymographs']);
if ~isfolder(out_dir), mkdir(out_dir); end

fig_ov = figure('Name','Detection Overview', ...
                'Position',[50 50 1300 520], 'Visible','off');

% Left: log-scaled mean frame with zeroth-order range and detections
subplot(1,2,1);
imagesc(log1p(mean_frame)); colormap(gca, gray); axis image; hold on;
% Draw the zeroth-order search range as a dashed yellow rectangle
rectangle('Position',[ZORDER_COL_MIN-0.5, 0.5, ...
                       ZORDER_COL_MAX-ZORDER_COL_MIN+1, height], ...
          'EdgeColor','y','LineWidth',1.4,'LineStyle','--');
xline(spec_start, 'c--', 'LineWidth',0.9, 'Alpha',0.7, ...
      'Label','spec start','FontSize',7,'LabelVerticalAlignment','top');
for i = 1:n_particles
    yline(particles(i,1), 'r-', 'LineWidth',0.9, 'Alpha',0.9);
    text(particles(i,3)+4, particles(i,1)-4, num2str(i), ...
         'Color','yellow','FontSize',8,'FontWeight','bold');
    % Cyan box centred on the per-particle (x, y)
    rectangle('Position',[particles(i,3)-WIN-0.5, particles(i,1)-WIN-0.5, ...
                           2*WIN+1, 2*WIN+1], 'EdgeColor','cyan','LineWidth',1);
end
xlabel('Column (x)'); ylabel('Row (y)');
title(sprintf('Time-averaged frame  (log scale)\nYellow box = ZO search range cols %d–%d', ...
              ZORDER_COL_MIN, ZORDER_COL_MAX), 'FontSize',9);
colorbar;

% Right: MIP-based spatial profile
subplot(1,2,2);
plot(double(spatial_profile), 1:height, 'b-','LineWidth',0.8); hold on;
xline(thresh_val,'r--','LineWidth',1.4, ...
      'DisplayName',sprintf('Threshold (%.0f)', thresh_val));
for i = 1:n_particles
    scatter(particles(i,2), particles(i,1), 70,'r','filled');
    text(particles(i,2)*1.02, particles(i,1)-2, num2str(i), ...
         'Color','r','FontSize',8,'FontWeight','bold');
end
set(gca,'YDir','reverse');
xlabel(sprintf('Max intensity across cols %d–%d (MIP)  (counts)', ...
               ZORDER_COL_MIN, ZORDER_COL_MAX));
ylabel('Row (spatial pixel)');
title(sprintf('%d particle(s) detected', n_particles),'FontSize',10);
legend('Spatial profile (MIP)','Threshold','Location','southeast','FontSize',8);

sgtitle(sprintf('Particle Detection Overview  |  ZO range cols %d–%d', ...
                ZORDER_COL_MIN, ZORDER_COL_MAX), ...
        'FontSize',12,'FontWeight','bold');
exportgraphics(fig_ov, fullfile(out_dir,'00_detection_overview.png'),'Resolution',150);
close(fig_ov);
fprintf('  Detection overview saved.\n');

%% ══════════════════════════════════════════════════════════════════════════════
%  STEP 4 — EXTRACT KYMOGRAPHS
%  mean_frame is used here — time-averaged photon counts are the correct
%  quantity for spectrum fitting and kymograph display.
%% ══════════════════════════════════════════════════════════════════════════════

fprintf('\n[4/5] Extracting kymographs\n');

kymographs_raw  = cell(n_particles,1);
kymographs_corr = cell(n_particles,1);
wavelength_axes = cell(n_particles,1);
gauss_params    = zeros(n_particles,4);   % [A, mu, sigma, C] per particle

for i = 1:n_particles

    row_c = particles(i,1);

    r0 = max(1,      row_c - HALF_WIDTH - VERTICAL_ALLOWANCE);
    r1 = min(height, row_c + HALF_WIDTH + VERTICAL_ALLOWANCE);

    % ── 4a. Average spectrum -> Gaussian fit to locate emission peak ──
    avg_spec = double(sum(mean_frame(r0:r1, spec_start:spec_end), 1));
    px_local = (1:n_spec_cols)';

    [A_fit, mu_fit, sigma_fit, C_fit] = fit_gaussian1d(px_local, avg_spec);

    if isnan(mu_fit) || sigma_fit <= 0
        fprintf('      Particle %d: first-order Gaussian failed, using fallback\n',i);
        [~, mu_fit]  = max(avg_spec);
        sigma_fit    = 10;
        A_fit        = avg_spec(mu_fit) - min(avg_spec);
        C_fit        = min(avg_spec);
    end

    gauss_params(i,:) = [A_fit, mu_fit, sigma_fit, C_fit];

    fprintf(['      Particle %d: first-order  mu = %.1f px (%.1f nm),' ...
             '  sigma = %.2f px\n'], i, mu_fit, px2nm(mu_fit), sigma_fit);

    % ── 4b. Extraction window = mu +/- 40*sigma ──────────────────────
    win_lo = max(1,           round(mu_fit - 40*sigma_fit));
    win_hi = min(n_spec_cols, round(mu_fit + 40*sigma_fit));
    n_win  = win_hi - win_lo + 1;

    abs_lo = spec_start + win_lo - 1;
    abs_hi = spec_start + win_hi - 1;

    % ── 4c. Build kymograph frame-by-frame ───────────────────────────
    kymo = zeros(n_frames, n_win, 'single');
    for f = 1:n_frames
        kymo(f,:) = sum(movie(r0:r1, abs_lo:abs_hi, f), 1);
    end
    kymographs_raw{i} = kymo;

    % ── 4d. Column-wise background subtraction ───────────────────────
    bg        = prctile(kymo, BG_PERCENTILE, 1);
    kymo_corr = max(kymo - bg, 0);
    kymographs_corr{i} = kymo_corr;

    % ── 4e. Wavelength axis ──────────────────────────────────────────
    wavelength_axes{i} = px2nm((win_lo:win_hi)');
end

%% ── SAVE MAT FILE ────────────────────────────────────────────────────────────

if SAVE_MAT
    save(fullfile(out_dir,'kymographs.mat'), ...
         'kymographs_raw','kymographs_corr','wavelength_axes', ...
         'particles','gauss_params','n_frames','FRAME_RATE', ...
         'ZORDER_COL_MIN','ZORDER_COL_MAX','spec_start','spec_end','px2nm','-v7.3');
    fprintf('      Arrays saved -> kymographs.mat\n');
end

%% ══════════════════════════════════════════════════════════════════════════════
%  STEP 5 — FIGURES & SUMMARY
%% ══════════════════════════════════════════════════════════════════════════════

fprintf('\n[5/5] Saving figures\n');

if ~isempty(FRAME_RATE) && FRAME_RATE > 0
    t_axis  = (0:n_frames-1) / FRAME_RATE;
    t_label = 'Time  (s)';
else
    t_axis  = 0:n_frames-1;
    t_label = 'Time  (frames)';
end

mean_wl_static = nan(n_particles,1);
std_wl_all     = nan(n_particles,1);

% Per-particle zeroth-order coordinates: x from detected column, y = row
zorder_coords = [particles(:,3), particles(:,1)];   % [x, y]

for i = 1:n_particles

    wl        = wavelength_axes{i};
    kymo_raw  = kymographs_raw{i};
    kymo_corr = kymographs_corr{i};

    clim_r = [prctile(kymo_raw(:),2),  prctile(kymo_raw(:),99.5)];
    clim_c = [prctile(kymo_corr(:),2), prctile(kymo_corr(:),99.5)];
    if clim_r(1)==clim_r(2), clim_r(2)=clim_r(1)+1; end
    if clim_c(1)==clim_c(2), clim_c(2)=clim_c(1)+1; end

    avg_spec  = mean(kymo_corr, 1);
    int_trace = sum(kymo_corr, 2);
    wl_col    = wl(:)';

    % ── Spectral stability: per-frame peak wavelength ─────────────────
    frame_totals  = sum(kymo_corr, 2);
    on_thresh_loc = prctile(frame_totals, 10);
    on_idx        = frame_totals > on_thresh_loc;

    if sum(on_idx) >= 2
        on_spec      = kymo_corr(on_idx,:);
        row_tot      = sum(on_spec, 2);
        row_tot(row_tot==0) = 1;
        per_frame_wl = (on_spec * wl_col') ./ row_tot;
        peak_wl_mean = mean(per_frame_wl);
        peak_wl_std  = std(per_frame_wl);
    else
        peak_wl_mean = px2nm(gauss_params(i,2));
        peak_wl_std  = Inf;
    end

    std_wl_all(i) = peak_wl_std;
    if peak_wl_std <= STATIC_THRESHOLD_NM
        mean_wl_static(i) = peak_wl_mean;
    end

    % ── Build figure ──────────────────────────────────────────────────
    fig = figure('Name',sprintf('Particle %d',i), ...
                 'Position',[50 50 1350 680],'Visible','off');

    ax1 = subplot(2,3,1);
    imagesc(ax1, wl, t_axis, kymo_raw, clim_r);
    colormap(ax1,COLORMAP); colorbar(ax1);
    set(ax1,'YDir','normal');
    xlabel(ax1,'Wavelength  (nm)'); ylabel(ax1,t_label);
    title(ax1,'Raw kymograph','FontSize',10);

    ax2 = subplot(2,3,2);
    imagesc(ax2, wl, t_axis, kymo_corr, clim_c);
    colormap(ax2,COLORMAP); colorbar(ax2);
    set(ax2,'YDir','normal');
    xlabel(ax2,'Wavelength  (nm)'); ylabel(ax2,t_label);
    title(ax2,'BG-subtracted kymograph','FontSize',10);

    if ~isnan(mean_wl_static(i))
        xline(ax1, mean_wl_static(i),'w--','LineWidth',1.0,'Alpha',0.7);
        xline(ax2, mean_wl_static(i),'w--','LineWidth',1.2,'Alpha',0.9, ...
              'Label',sprintf('%.1f nm',mean_wl_static(i)), ...
              'LabelVerticalAlignment','bottom','FontSize',7,'Color','w');
    end

    ax3 = subplot(2,3,3);
    plot(ax3, wl, avg_spec,'b-','LineWidth',1.3); hold(ax3,'on');
    mu_f  = gauss_params(i,2);
    s_f   = gauss_params(i,3);
    A_f   = gauss_params(i,1);
    wlo_i = max(1, round(mu_f-40*s_f));
    whi_i = min(n_spec_cols, round(mu_f+40*s_f));
    px_ov = (wlo_i:whi_i)';
    plot(ax3, px2nm(px_ov), A_f*exp(-0.5*((px_ov-mu_f)/s_f).^2), ...
         'r--','LineWidth',1.0,'DisplayName','Gaussian fit');
    xline(ax3, px2nm(mu_f),'g:','LineWidth',1.2,'DisplayName','\mu');
    xlabel(ax3,'Wavelength  (nm)'); ylabel(ax3,'Mean counts');
    title(ax3,'Time-averaged spectrum','FontSize',10);
    legend(ax3,'Spectrum','Gaussian','\mu','Location','best','FontSize',7);
    grid(ax3,'on');

    ax4 = subplot(2,3,[4 5 6]);
    plot(ax4, t_axis, int_trace,'Color',[0.2 0.4 0.8],'LineWidth',0.9);
    hold(ax4,'on');
    fill(ax4,[t_axis,fliplr(t_axis)],[int_trace',zeros(1,n_frames)], ...
         [0.2 0.4 0.8],'FaceAlpha',0.2,'EdgeColor','none');
    xlabel(ax4,t_label); ylabel(ax4,'Integrated PL  (counts)');
    title(ax4,'PL intensity trajectory','FontSize',10);
    grid(ax4,'on'); xlim(ax4,[t_axis(1),t_axis(end)]);

    coord_str = sprintf('ZO: (x=%d, y=%d)', zorder_coords(i,1), zorder_coords(i,2));
    if ~isnan(mean_wl_static(i))
        wl_str = sprintf(['Mean \\lambda = %.1f nm  ' ...
                          '(stable, \\sigma_{\\lambda} = %.1f nm < %.0f nm)'], ...
                         mean_wl_static(i), std_wl_all(i), STATIC_THRESHOLD_NM);
    else
        wl_str = sprintf('Spectral diffusion  (\\sigma_{\\lambda} = %.1f nm > %.0f nm)', ...
                         std_wl_all(i), STATIC_THRESHOLD_NM);
    end
    sgtitle(sprintf('Particle %d  |  %s\n%s', i, coord_str, wl_str), ...
            'FontSize',10,'FontWeight','bold');

    exportgraphics(fig, fullfile(out_dir,sprintf('particle_%03d_kymograph.png',i)), ...
                   'Resolution',150);
    close(fig);
    fprintf('    Saved -> particle_%03d_kymograph.png\n', i);

    % ── Multi-Gaussian decomposition ─────────────────────────────────────
    fprintf('    Particle %d: running spectral decomposition...\n', i);

    dec = decompose_spectrum_kymograph( ...
              kymographs_corr{i}, wavelength_axes{i}, FRAME_RATE);

    % Only use frames where the dot was on (mu1 not NaN)
    on_mask = ~isnan(dec.mu1);

    % ── Figure 1: ZPL trajectory + Huang-Rhys over time ─────────────
    fig_dec1 = figure('Position',[50 50 1200 500],'Visible','off');

    ax_mu1 = subplot(1,2,1);
    plot(ax_mu1, dec.t_axis(on_mask), dec.mu1(on_mask), ...
         'Color',[0.8 0.2 0.2],'LineWidth',0.8);
    xlabel(ax_mu1, t_label);
    ylabel(ax_mu1, 'ZPL centre  \lambda_1  (nm)');
    title(ax_mu1, 'Spectral diffusion trajectory','FontSize',10);
    grid(ax_mu1,'on');

    ax_hr = subplot(1,2,2);
    plot(ax_hr, dec.t_axis(on_mask), dec.HR_factor(on_mask), ...
         'Color',[0.2 0.6 0.2],'LineWidth',0.8);
    yline(ax_hr, mean(dec.HR_factor(on_mask),'omitnan'), ...
          'k--','LineWidth',1.2, ...
          'Label',sprintf('Mean S = %.2f', mean(dec.HR_factor(on_mask),'omitnan')));
    xlabel(ax_hr, t_label);
    ylabel(ax_hr, 'Huang-Rhys proxy  S = A_2/A_1');
    title(ax_hr, 'Electron-phonon coupling','FontSize',10);
    grid(ax_hr,'on');

    sgtitle(sprintf('Particle %d — Spectral Decomposition', i), ...
            'FontSize',11,'FontWeight','bold');
    exportgraphics(fig_dec1, ...
        fullfile(out_dir, sprintf('particle_%03d_decomposition_trajectories.png',i)), ...
        'Resolution',150);
    close(fig_dec1);

    % ── Figure 2: Amplitude time series for all N components ─────────
    fig_dec2 = figure('Position',[50 50 1200 400],'Visible','off');
    colors_comp = lines(dec.N);
    comp_labels = arrayfun(@(k) sprintf('Component %d (\\Delta=%.1f nm)',k, ...
                           dec.delta_k(k)), 1:dec.N, 'UniformOutput',false);

    hold on;
    for k = 1:dec.N
        plot(dec.t_axis(on_mask), dec.amplitudes(on_mask,k), ...
             'Color',colors_comp(k,:),'LineWidth',0.8, ...
             'DisplayName',comp_labels{k});
    end
    xlabel(t_label);
    ylabel('Amplitude  A_k  (counts)');
    title(sprintf('Particle %d — Component amplitude trajectories',i), ...
          'FontSize',10,'FontWeight','bold');
    legend('Location','best','FontSize',8);
    grid on;
    exportgraphics(fig_dec2, ...
        fullfile(out_dir, sprintf('particle_%03d_decomposition_amplitudes.png',i)), ...
        'Resolution',150);
    close(fig_dec2);

    % ── Figure 3: Example decomposed spectrum (mean of on-frames) ────
    fig_dec3 = figure('Position',[50 50 700 400],'Visible','off');
    wl_i     = wavelength_axes{i};
    mean_on  = mean(kymographs_corr{i}(on_mask,:), 1)';

    % Reconstruct components at mean mu1
    mu1_mean  = mean(dec.mu1(on_mask), 'omitnan');
    A_mean    = mean(dec.amplitudes(on_mask,:), 1);
    wl_fine   = linspace(min(wl_i), max(wl_i), 500)';
    I_total   = zeros(size(wl_fine));

    hold on;
    plot(wl_i, mean_on, 'k-', 'LineWidth',1.4,'DisplayName','Observed');
    for k = 1:dec.N
        mu_k   = mu1_mean + dec.delta_k(k);
        I_comp = A_mean(k) * exp(-0.5*((wl_fine-mu_k)/dec.sigma_k(k)).^2);
        I_total = I_total + I_comp;
        patch([wl_fine; flipud(wl_fine)], [I_comp; zeros(size(I_comp))], ...
              colors_comp(k,:), 'FaceAlpha',0.35,'EdgeColor','none', ...
              'DisplayName',comp_labels{k});
    end
    plot(wl_fine, I_total,'r--','LineWidth',1.2,'DisplayName','Total fit');

    xlabel('Wavelength  (nm)'); ylabel('Mean counts');
    title(sprintf('Particle %d — Decomposed mean spectrum',i),'FontSize',10);
    legend('Location','best','FontSize',8);
    grid on;
    exportgraphics(fig_dec3, ...
        fullfile(out_dir, sprintf('particle_%03d_decomposition_spectrum.png',i)), ...
        'Resolution',150);
    close(fig_dec3);

    % ── Save decomposition to .mat ────────────────────────────────────
    dec_file = fullfile(out_dir, sprintf('particle_%03d_decomposition.mat',i));
    save(dec_file, 'dec', '-v7.3');
    fprintf('    Decomposition saved -> particle_%03d_decomposition*.png\n', i);

end

% ── Mosaic ────────────────────────────────────────────────────────────────────

ncols_m = min(4, n_particles);
nrows_m = ceil(n_particles/ncols_m);
fig_mos = figure('Name','Mosaic', ...
                 'Position',[50 50 ncols_m*320 nrows_m*280],'Visible','off');
for i = 1:n_particles
    ax     = subplot(nrows_m, ncols_m, i);
    kymo_c = kymographs_corr{i};
    wl     = wavelength_axes{i};
    clim_c = [prctile(kymo_c(:),2), prctile(kymo_c(:),99.5)];
    if clim_c(1)==clim_c(2), clim_c(2)=clim_c(1)+1; end
    imagesc(ax, wl, t_axis, kymo_c, clim_c);
    colormap(ax,COLORMAP); set(ax,'YDir','normal');
    if ~isnan(mean_wl_static(i))
        ttl = sprintf('P%d (%d,%d) %.0fnm*', i, zorder_coords(i,1), ...
                      zorder_coords(i,2), mean_wl_static(i));
    else
        ttl = sprintf('P%d (%d,%d) diff~', i, zorder_coords(i,1), ...
                      zorder_coords(i,2));
    end
    title(ax,ttl,'FontSize',7);
    xlabel(ax,'\lambda (nm)','FontSize',7);
    ylabel(ax,t_label,'FontSize',7);
    ax.FontSize = 7;
end
sgtitle('All Particles — Kymograph Mosaic','FontSize',11,'FontWeight','bold');
exportgraphics(fig_mos, fullfile(out_dir,'00_all_kymographs_mosaic.png'),'Resolution',120);
close(fig_mos);
fprintf('  Mosaic saved.\n');

% ── Summary table ─────────────────────────────────────────────────────────────

fprintf('\n%s\n PARTICLE SUMMARY\n%s\n', repmat('-',1,70), repmat('-',1,70));
fprintf(' %-5s  %-14s  %-12s  %-10s  %-22s\n', ...
        'No.','ZO coord (x,y)','Mean lambda','SD (nm)','Stability');
fprintf('%s\n', repmat('-',1,70));

fid = fopen(fullfile(out_dir,'00_particle_summary.txt'),'w');
fprintf(fid,'No.\tZO_x\tZO_y\tMean_wl_nm\tSD_nm\tStability\n');

for i = 1:n_particles
    cx = zorder_coords(i,1);  cy = zorder_coords(i,2);
    if ~isnan(mean_wl_static(i))
        wl_d   = sprintf('%.1f nm', mean_wl_static(i));
        sd_d   = sprintf('%.2f nm', std_wl_all(i));
        st_d   = sprintf('STATIC (<= %.0f nm)', STATIC_THRESHOLD_NM);
        fid_wl = sprintf('%.2f', mean_wl_static(i));
    else
        wl_d   = 'N/A';
        sd_d   = sprintf('%.2f nm', std_wl_all(i));
        st_d   = sprintf('DIFFUSING (> %.0f nm)', STATIC_THRESHOLD_NM);
        fid_wl = 'NaN';
    end
    fprintf(' %-5d  (%-5d, %-5d)  %-12s  %-10s  %s\n', ...
            i, cx, cy, wl_d, sd_d, st_d);
    fprintf(fid,'%d\t%d\t%d\t%s\t%.2f\t%s\n', ...
            i, cx, cy, fid_wl, std_wl_all(i), st_d);
end
fprintf('%s\n', repmat('-',1,70));
fclose(fid);
fprintf('  Summary saved -> 00_particle_summary.txt\n');

fprintf('\n%s\n Done!  All outputs -> %s\n%s\n\n', ...
        repmat('=',1,60), out_dir, repmat('=',1,60));

%% ══════════════════════════════════════════════════════════════════════════════
%  LOCAL HELPER: fit_gaussian1d
%  Fits  y = A*exp(-0.5*((x-mu)/sigma)^2) + C  via weighted moments.
%  No iterative solver, no toolbox.
%  Returns NaN for all outputs if the fit is degenerate.
%% ══════════════════════════════════════════════════════════════════════════════

function [A, mu, sigma, C] = fit_gaussian1d(x, y)
    x = double(x(:));  y = double(y(:));
    C  = min(y);
    yc = y - C;
    S  = sum(yc);
    if S <= 0
        A = NaN; mu = NaN; sigma = NaN; return
    end
    mu    = sum(x .* yc) / S;
    var_g = sum((x - mu).^2 .* yc) / S;
    if var_g <= 0
        A = NaN; mu = NaN; sigma = NaN; return
    end
    sigma = sqrt(var_g);
    A     = max(yc);
    if isnan(mu) || isnan(sigma) || sigma == 0
        A = NaN; mu = NaN; sigma = NaN;
    end
end

%% ══════════════════════════════════════════════════════════════════════════════
%  LOCAL HELPER: decompose_spectrum_kymograph
%  Fits a sum of N Gaussian components to every frame of a kymograph.
%% ══════════════════════════════════════════════════════════════════════════════

function results = decompose_spectrum_kymograph(kymo_corr, wl, FRAME_RATE)
%DECOMPOSE_SPECTRUM_KYMOGRAPH
%  Fits a sum of N Gaussian components to every frame of a kymograph.
%
%  Each component is centred at mu_1 + delta_k where delta_k are fixed
%  phonon offsets. Amplitudes are solved via linear least squares for each
%  candidate mu_1. The outer search over mu_1 uses fminbnd (no toolbox).
%
%  INPUTS
%    kymo_corr  [n_frames x n_wl]  background-subtracted kymograph
%    wl         [n_wl x 1]         wavelength axis in nm
%    FRAME_RATE scalar Hz (or [] for frame-number axis)
%
%  OUTPUT  results struct with fields:
%    mu1        [n_frames x 1]  ZPL centre per frame (nm)
%    amplitudes [n_frames x N]  amplitude of each component per frame
%    residuals  [n_frames x 1]  fit residual (goodness of fit)
%    HR_factor  [n_frames x 1]  Huang-Rhys proxy = A2/A1
%    t_axis     [n_frames x 1]  time axis (s or frames)

    % ── Model parameters — adjust these to your system ──────────────────
    delta_k_nm  = [0,  -5.8,  -11.6];   % nm offsets — N = 3 components
    sigma_k_nm  = [4.0,  5.0,   6.0];   % fixed widths per component (nm)
    N           = numel(delta_k_nm);

    % ── Setup ────────────────────────────────────────────────────────────
    n_frames = size(kymo_corr, 1);
    wl       = wl(:);           % ensure column vector
    n_wl     = numel(wl);

    % Search range for mu_1: centre +/- 20 nm from the mean-spectrum peak
    mean_spec   = mean(kymo_corr, 1)';
    [~, pk_idx] = max(mean_spec);
    mu1_center  = wl(pk_idx);
    mu1_lo      = mu1_center - 20;
    mu1_hi      = mu1_center + 20;

    % Pre-allocate outputs
    mu1_out  = zeros(n_frames, 1);
    amp_out  = zeros(n_frames, N);
    res_out  = zeros(n_frames, 1);

    % ── Frame-by-frame fitting ────────────────────────────────────────────
    for f = 1:n_frames

        I_obs = double(kymo_corr(f, :))';   % [n_wl x 1]

        % Skip near-zero frames (dot is off)
        if max(I_obs) < 1e-3 * max(mean_spec)
            mu1_out(f) = NaN;
            amp_out(f,:) = 0;
            res_out(f) = NaN;
            continue
        end

        % ── Outer 1D search over mu1 using fminbnd ───────────────────
        % fminbnd is base MATLAB — no toolbox needed
        [mu1_best, ssr_best] = fminbnd(@fit_given_mu1, mu1_lo, mu1_hi);

        % Recover amplitudes at the best mu1
        [~, A_best] = fit_given_mu1(mu1_best);

        mu1_out(f)   = mu1_best;
        amp_out(f,:) = A_best';
        res_out(f)   = ssr_best;
    end

    % ── Build time axis ──────────────────────────────────────────────────
    if ~isempty(FRAME_RATE) && FRAME_RATE > 0
        t_axis = (0:n_frames-1)' / FRAME_RATE;
    else
        t_axis = (0:n_frames-1)';
    end

    % ── Huang-Rhys proxy ─────────────────────────────────────────────────
    % S ≈ A2/A1 (ratio of 1-phonon to ZPL amplitude)
    % Protect against divide-by-zero on blinked-off frames
    HR = amp_out(:,2) ./ max(amp_out(:,1), 1e-9);
    HR(isnan(mu1_out)) = NaN;

    % ── Pack results ─────────────────────────────────────────────────────
    results.mu1        = mu1_out;
    results.amplitudes = amp_out;
    results.residuals  = res_out;
    results.HR_factor  = HR;
    results.t_axis     = t_axis;
    results.delta_k    = delta_k_nm;
    results.sigma_k    = sigma_k_nm;
    results.N          = N;

    % ── NESTED FUNCTION (Moved outside the loop) ─────────────────────────
    % Because this is nested inside decompose_spectrum_kymograph, it 
    % automatically shares the parent's variables (like I_obs and wl).
    function [ssr, A_fit] = fit_given_mu1(mu1_candidate)
        % Build the design matrix Phi [n_wl x N]
        Phi = zeros(n_wl, N);
        for k = 1:N
            mu_k     = mu1_candidate + delta_k_nm(k);
            Phi(:,k) = exp(-0.5 * ((wl - mu_k) / sigma_k_nm(k)).^2);
        end

        % Constrained linear LS: amplitudes >= 0
        % Use lsqnonneg (base MATLAB, no toolbox)
        A_fit = lsqnonneg(Phi, I_obs);

        % Sum of squared residuals
        residual = I_obs - Phi * A_fit;
        ssr      = sum(residual.^2);
    end

end