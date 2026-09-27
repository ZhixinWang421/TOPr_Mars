function [topr_series,mc] = TOPr_Mars(data,Nmc,cfg)
%TOPR_MARS Monte Carlo parameter sensitivity for Mars orbital power.
%
%   TOPR_SERIES = TOPR_MARS(DATA)
%   [TOPR_SERIES,MC] = TOPR_MARS(DATA,NMC)
%   [TOPR_SERIES,MC] = TOPR_MARS(DATA,NMC,CFG)
%
% DATA must be an evenly sampled N-by-2 matrix:
%   column 1  time in kyr, strictly increasing
%   column 2  observed value
%
% For each moving window, the total orbital power ratio (TOPr) is
%
%       power in the union of the target orbital-frequency bands
%       --------------------------------------------------------- .
%                  power between FtMin and FtMax
%
% The Monte Carlo ensemble samples plausible analysis parameters. By
% default it varies the running-window length and one shared bandwidth
% factor for all nine target periods. It does not add noise to DATA,
% recenter the ensemble, clip the results, or transform TOPr to 1-TOPr.
% The percentile envelopes therefore describe conditional sensitivity to
% the configured method choices; they are not p values, measurement-error
% intervals, or age-model uncertainty intervals.
%
% TOPR_SERIES contains [time_kyr, median_topr] on the common-support
% interval. MC contains the full percentile summaries, parameter draws,
% fixed-parameter reference, realized band limits, configuration, and (by
% default) all Monte Carlo realizations.
%
% Default target periods (kyr):
%   [125 94.6 54.4 1260 116.4 65.8 118.4 52.6 37.2]
%
% Default analysis settings:
%   reference window             1260 kyr
%   Monte Carlo window range     945--1575 kyr
%   reference and MC NW          2
%   reference bandwidth factor   1.2
%   MC bandwidth-factor range    0.9--1.2
%   total-power frequency range  0--0.1 kyr^-1
%   NFFT                          4096
%   output percentiles            5, 25, 50, 75, 95
%   random seed                   20260912
%
% Example:
%   data = readmatrix('mars_series.txt');
%   [topr,mc] = TOPr_Mars(data,1000);
%
% Example enabling sampling-interval and NW sensitivity:
%   cfg = struct('SampleRangeKyr',[4.5 5.5], ...
%                'NWValues',[2 2.5 3], ...
%                'Seed',20260912, ...
%                'PlotBaseline',true);
%   [topr,mc] = TOPr_Mars(data,1000,cfg);
%
% For downstream spectral uncertainty, analyze every column of
% MC.topr_simulations and summarize the resulting spectra. The spectrum of
% a pointwise median curve need not equal the median spectrum.
%
% Requirement: MATLAB Signal Processing Toolbox (PMTM).
% Version: 1.0.0 (2026-09-28).

if nargin<2 || isempty(Nmc)
    Nmc = 1000;
end
if nargin<3 || isempty(cfg)
    cfg = struct();
end

validateattributes(Nmc,{'numeric'}, ...
    {'scalar','integer','positive','finite'},mfilename,'Nmc');
data = validate_input_data(data);
cfg = resolve_config(data,cfg);
stream = RandStream('mt19937ar','Seed',cfg.Seed);

%% Fixed-parameter reference and common output grid.
baselineDt = median(diff(data(:,1)));
baselineBands = target_bands(cfg.TargetPeriodsKyr, ...
    cfg.ReferenceWindowKyr,cfg.ReferenceNW, ...
    cfg.ReferenceBandWidthFactor,cfg.FtMin,cfg.FtMax, ...
    1/(2*baselineDt));
pow0 = compute_topr(data,baselineBands,cfg.ReferenceWindowKyr, ...
    cfg.ReferenceNW,cfg.FtMin,cfg.FtMax,cfg.StepSamples,cfg.Nfft);
baselineTime = pow0(:,1);
baselineTOPr = pow0(:,2);

%% Reproducible parameter draws.
if isempty(cfg.SampleRangeKyr)
    sampleDraw = repmat(baselineDt,Nmc,1);
else
    sampleDraw = uniform_draw(stream,cfg.SampleRangeKyr,Nmc);
end
windowDraw = uniform_draw(stream,cfg.WindowRangeKyr,Nmc);
bandFactorDraw = uniform_draw(stream,cfg.BandWidthFactorRange,Nmc);
nwIndex = randi(stream,numel(cfg.NWValues),Nmc,1);
nwDraw = reshape(cfg.NWValues(nwIndex),[],1);

rawTOPr = nan(numel(baselineTime),Nmc);
lowCutoff = nan(Nmc,numel(cfg.TargetPeriodsKyr));
highCutoff = nan(Nmc,numel(cfg.TargetPeriodsKyr));

progressEvery = cfg.ProgressEvery;
if progressEvery<=0
    progressEvery = max(1,round(Nmc/20));
end
if cfg.ShowProgress
    fprintf(['[TOPr_Mars] N=%d, seed=%d, window=%.6g--%.6g kyr, ', ...
        'band factor=%.6g--%.6g\n'],Nmc,cfg.Seed, ...
        cfg.WindowRangeKyr(1),cfg.WindowRangeKyr(2), ...
        cfg.BandWidthFactorRange(1),cfg.BandWidthFactorRange(2));
end

for i = 1:Nmc
    sample = sampleDraw(i);
    realizationData = resample_if_needed(data,sample,baselineDt);

    [bands,lo,hi] = target_bands(cfg.TargetPeriodsKyr,windowDraw(i), ...
        nwDraw(i),bandFactorDraw(i),cfg.FtMin,cfg.FtMax,1/(2*sample));
    lowCutoff(i,:) = lo;
    highCutoff(i,:) = hi;

    powi = compute_topr(realizationData,bands,windowDraw(i),nwDraw(i), ...
        cfg.FtMin,cfg.FtMax,cfg.StepSamples,cfg.Nfft);
    rawTOPr(:,i) = interp1(powi(:,1),powi(:,2),baselineTime, ...
        'linear',NaN);

    if cfg.ShowProgress && (mod(i,progressEvery)==0 || i==Nmc)
        fprintf('  completed %d / %d\n',i,Nmc);
    end
end

%% Pointwise summaries and representative curve.
percentiles = cfg.Percentiles(:).';
[toprPercentiles,validCount] = row_percentiles(rawTOPr,percentiles);
minimumValid = max(1,ceil(cfg.MinValidFraction*Nmc));
medianIndex = find_percentile(percentiles,50);
commonMask = validCount>=minimumValid & ...
    isfinite(toprPercentiles(:,medianIndex));
if ~any(commonMask)
    error('TOPr_Mars:NoCommonOutput', ...
        ['No output time has the requested Monte Carlo coverage. Narrow ', ...
         'the parameter ranges or reduce CFG.MinValidFraction.']);
end

medianTOPr = toprPercentiles(:,medianIndex);
topr_series = [baselineTime(commonMask),medianTOPr(commonMask)];

parameterMatrix = [sampleDraw,windowDraw,nwDraw,bandFactorDraw];
parameterNames = {'sample_interval_kyr','window_kyr','NW', ...
    'shared_band_width_factor'};
parameterSpan = max(parameterMatrix,[],1)-min(parameterMatrix,[],1);
parameterScale = max(abs(parameterMatrix),[],1);
varyingParameter = parameterSpan>100*eps(max(parameterScale,1));

mc = struct();
mc.function_name = mfilename;
mc.version = '1.0.0';
mc.time_kyr = baselineTime;
mc.percentiles = percentiles;
mc.topr_percentiles = percentile_table( ...
    baselineTime,toprPercentiles,percentiles);
mc.median_topr = medianTOPr;
mc.representative_series = array2table(topr_series, ...
    'VariableNames',{'time_kyr','median_topr'});
mc.common_support_mask = commonMask;
mc.valid_realizations_per_time = validCount;
mc.deterministic_baseline = table(baselineTime,baselineTOPr, ...
    'VariableNames',{'time_kyr','topr'});
mc.parameter_draws = array2table(parameterMatrix, ...
    'VariableNames',parameterNames);
mc.varying_parameters = parameterNames(varyingParameter);
mc.fixed_parameters = parameterNames(~varyingParameter);
mc.realized_lower_cutoff_kyr_inv = array2table(lowCutoff, ...
    'VariableNames',band_variable_names(cfg.TargetPeriodsKyr,'low'));
mc.realized_upper_cutoff_kyr_inv = array2table(highCutoff, ...
    'VariableNames',band_variable_names(cfg.TargetPeriodsKyr,'high'));
mc.config = cfg;
mc.definition = ['Raw total orbital power ratio evaluated over the union ', ...
    'of the configured target-frequency bands and normalized by total ', ...
    'power from FtMin to FtMax.'];
mc.note = ['The percentile envelopes describe sensitivity to the sampled ', ...
    'analysis parameters. For spectral uncertainty, analyze each Monte ', ...
    'Carlo realization rather than only the pointwise median curve.'];

if cfg.ReturnSimulations
    mc.topr_simulations = rawTOPr;
else
    mc.topr_simulations = [];
end

if cfg.MakePlot
    plot_uncertainty(baselineTime(commonMask), ...
        toprPercentiles(commonMask,:),baselineTOPr(commonMask), ...
        percentiles,Nmc,cfg.PlotBaseline);
end
end

function data = validate_input_data(data)
if ~(isnumeric(data) && isreal(data) && ismatrix(data) && ...
        size(data,2)==2 && size(data,1)>=20 && all(isfinite(data(:))))
    error('TOPr_Mars:BadData', ...
        'DATA must be a finite real N-by-2 [time_kyr,value] matrix.');
end
data = double(data);
increments = diff(data(:,1));
if any(increments<=0)
    error('TOPr_Mars:BadTime', ...
        'DATA time coordinates must be strictly increasing.');
end
referenceStep = median(increments);
if max(abs(increments-referenceStep))>1e-6*max(1,referenceStep)
    error('TOPr_Mars:UnevenTime', ...
        ['DATA must be evenly sampled. Resample the series before calling ', ...
         'TOPr_Mars, or use CFG.SampleRangeKyr only for sensitivity tests ', ...
         'after providing an evenly sampled baseline series.']);
end
end

function cfg = resolve_config(data,user)
defaults = struct();
defaults.TargetPeriodsKyr = ...
    [125 94.6 54.4 1260 116.4 65.8 118.4 52.6 37.2];
defaults.ReferenceWindowKyr = 1260;
defaults.ReferenceNW = 2;
defaults.ReferenceBandWidthFactor = 1.2;
defaults.WindowRangeKyr = [945,1575];
defaults.NWValues = 2;
defaults.BandWidthFactorRange = [0.9,1.2];
defaults.SampleRangeKyr = [];
defaults.FtMin = 0;
defaults.FtMax = 0.1;
defaults.StepSamples = 10;
defaults.Nfft = 4096;
defaults.Seed = 20260912;
defaults.Percentiles = [5,25,50,75,95];
defaults.MinValidFraction = 1;
defaults.ReturnSimulations = true;
defaults.MakePlot = true;
defaults.PlotBaseline = false;
defaults.ShowProgress = true;
defaults.ProgressEvery = 0;

if ~isstruct(user) || ~isscalar(user)
    error('TOPr_Mars:BadConfig','CFG must be a scalar structure.');
end
cfg = defaults;
names = fieldnames(user);
for i = 1:numel(names)
    if ~isfield(defaults,names{i})
        error('TOPr_Mars:UnknownConfig', ...
            'Unknown configuration field: %s',names{i});
    end
    cfg.(names{i}) = user.(names{i});
end

validateattributes(cfg.TargetPeriodsKyr,{'numeric'}, ...
    {'vector','nonempty','finite','positive'});
validateattributes(cfg.ReferenceWindowKyr,{'numeric'}, ...
    {'scalar','finite','positive'});
validateattributes(cfg.ReferenceNW,{'numeric'}, ...
    {'scalar','finite','>=',0.5});
validateattributes(cfg.ReferenceBandWidthFactor,{'numeric'}, ...
    {'scalar','finite','positive'});
validate_two_value_range(cfg.WindowRangeKyr,'WindowRangeKyr');
validateattributes(cfg.NWValues,{'numeric'}, ...
    {'vector','nonempty','finite','>=',0.5});
validate_two_value_range(cfg.BandWidthFactorRange,'BandWidthFactorRange');
if ~isempty(cfg.SampleRangeKyr)
    validate_two_value_range(cfg.SampleRangeKyr,'SampleRangeKyr');
end
validateattributes(cfg.FtMin,{'numeric'}, ...
    {'scalar','finite','nonnegative'});
validateattributes(cfg.FtMax,{'numeric'}, ...
    {'scalar','finite','positive'});
validateattributes(cfg.StepSamples,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(cfg.Nfft,{'numeric'}, ...
    {'scalar','integer','>=',2});
validateattributes(cfg.Seed,{'numeric'}, ...
    {'scalar','integer','nonnegative','<=',2^32-1});
validateattributes(cfg.Percentiles,{'numeric'},{'vector','finite'});
validateattributes(cfg.MinValidFraction,{'numeric'}, ...
    {'scalar','finite','>',0,'<=',1});
validateattributes(cfg.ProgressEvery,{'numeric'}, ...
    {'scalar','integer','nonnegative'});

logicalFields = {'ReturnSimulations','MakePlot','PlotBaseline','ShowProgress'};
for i = 1:numel(logicalFields)
    cfg.(logicalFields{i}) = validate_logical_scalar( ...
        cfg.(logicalFields{i}),logicalFields{i});
end

if cfg.FtMax<=cfg.FtMin
    error('TOPr_Mars:BadFrequencyRange', ...
        'CFG.FtMax must be larger than CFG.FtMin.');
end
if any(cfg.Percentiles<0 | cfg.Percentiles>100) || ...
        ~any(abs(cfg.Percentiles-50)<1e-12)
    error('TOPr_Mars:BadPercentiles', ...
        'CFG.Percentiles must lie in [0,100] and include 50.');
end
if cfg.MakePlot && ...
        ~all(arrayfun(@(p) any(abs(cfg.Percentiles-p)<1e-12), ...
        [5,25,50,75,95]))
    error('TOPr_Mars:PlotPercentiles', ...
        'Plotting requires CFG.Percentiles to include [5 25 50 75 95].');
end

recordLength = data(end,1)-data(1,1);
if cfg.WindowRangeKyr(2)>=recordLength || ...
        cfg.ReferenceWindowKyr>=recordLength
    error('TOPr_Mars:WindowTooLong', ...
        'Reference and Monte Carlo windows must be shorter than the record.');
end

baselineDt = median(diff(data(:,1)));
if isempty(cfg.SampleRangeKyr)
    minimumSample = baselineDt;
    maximumSample = baselineDt;
else
    minimumSample = min(baselineDt,cfg.SampleRangeKyr(1));
    maximumSample = max(baselineDt,cfg.SampleRangeKyr(2));
end
minimumNyquist = 1/(2*maximumSample);
if cfg.FtMin>=min(cfg.FtMax,minimumNyquist)
    error('TOPr_Mars:FrequencyAboveNyquist', ...
        ['CFG.FtMin leaves no total-power interval for at least one ', ...
         'configured sampling interval.']);
end
minimumWindowSamples = fix(cfg.WindowRangeKyr(1)/maximumSample);
referenceWindowSamples = fix(cfg.ReferenceWindowKyr/baselineDt);
if minimumWindowSamples<2 || referenceWindowSamples<2
    error('TOPr_Mars:WindowTooShort', ...
        'A configured window contains fewer than two samples.');
end
maximumWindowSamples = max(referenceWindowSamples, ...
    fix(cfg.WindowRangeKyr(2)/minimumSample));
if cfg.Nfft<maximumWindowSamples
    error('TOPr_Mars:NfftTooSmall', ...
        ['CFG.Nfft must be at least the largest possible window sample ', ...
         'count (%d).'],maximumWindowSamples);
end

cfg.TargetPeriodsKyr = double(cfg.TargetPeriodsKyr(:).');
cfg.NWValues = unique(double(cfg.NWValues(:).'));
cfg.WindowRangeKyr = double(cfg.WindowRangeKyr(:).');
cfg.BandWidthFactorRange = double(cfg.BandWidthFactorRange(:).');
if ~isempty(cfg.SampleRangeKyr)
    cfg.SampleRangeKyr = double(cfg.SampleRangeKyr(:).');
end
cfg.Percentiles = unique(double(cfg.Percentiles(:).'));
end

function value = validate_logical_scalar(value,name)
if ~(isscalar(value) && (islogical(value) || isnumeric(value)) && ...
        isreal(value) && isfinite(double(value)) && ...
        ismember(double(value),[0,1]))
    error('TOPr_Mars:BadLogicalConfig', ...
        'CFG.%s must be a scalar logical value.',name);
end
value = logical(value);
end

function validate_two_value_range(value,name)
if ~(isnumeric(value) && isreal(value) && numel(value)==2 && ...
        all(isfinite(value)))
    error('TOPr_Mars:BadRange', ...
        'CFG.%s must contain two finite numeric values.',name);
end
value = value(:).';
if value(2)<value(1) || any(value<=0)
    error('TOPr_Mars:BadRange', ...
        'CFG.%s must be an ordered positive [minimum maximum] range.',name);
end
end

function data = resample_if_needed(data,sample,baselineDt)
if abs(sample-baselineDt)<=100*eps(max(sample,baselineDt))
    return
end
coordinates = (data(1,1):sample:data(end,1)).';
values = interp1(data(:,1),data(:,2),coordinates,'pchip');
data = [coordinates,values];
end

function values = uniform_draw(stream,limits,n)
limits = double(limits(:).');
if limits(1)==limits(2)
    values = repmat(limits(1),n,1);
else
    values = limits(1)+diff(limits).*rand(stream,n,1);
end
end

function [bands,low,high] = target_bands(periods,window,nw,factor, ...
        ftmin,ftmax,nyquist)
center = 1./periods;
halfWidth = factor*nw/window;
upperLimit = min(ftmax,nyquist);
low = max(center-halfWidth,ftmin);
high = min(center+halfWidth,upperLimit);
if any(high<=low)
    first = find(high<=low,1);
    error('TOPr_Mars:EmptyBand', ...
        ['The configured parameters produce an empty target band for ', ...
         'period %.6g kyr.'],periods(first));
end
bands = reshape([low;high],1,[]);
end

function pow = compute_topr(data,bands,window,nw,ftmin,ftmax,step,nfft)
% Vectorized moving-window multitaper power decomposition.
dt = median(diff(data(:,1)));
fs = 1/dt;
nyquist = fs/2;
nrow = size(data,1);
npts = fix(window/dt);
if npts<2 || npts>nrow
    error('TOPr_Mars:BadWindowSamples', ...
        'A window produces an invalid number of samples.');
end
if nw<0.5 || nw>=npts/2
    error('TOPr_Mars:BadNW', ...
        'NW must be at least 0.5 and smaller than half the window size.');
end
if nfft<npts
    error('TOPr_Mars:NfftTooSmall', ...
        'NFFT must be at least the number of samples in every window.');
end
if step>=nrow/2
    error('TOPr_Mars:StepTooLarge', ...
        'CFG.StepSamples is too large for a realization.');
end

upperTotal = min(ftmax,nyquist);
if upperTotal<=ftmin
    error('TOPr_Mars:EmptyTotalBand', ...
        'No total-power frequency interval remains below Nyquist.');
end

windowCount = fix((nrow-npts)/step)+1;
starts = 1+step*(0:windowCount-1);
indices = (0:npts-1).'+starts;
values = data(:,2);
segments = values(indices);
segments = segments-mean(segments,1);
[powerByFrequency,frequency] = pmtm(segments,nw,nfft,fs);
frequency = frequency(:);

totalMask = frequency>=ftmin & frequency<=upperTotal;
targetMask = false(size(frequency));
for b = 1:numel(bands)/2
    low = min(bands(2*b-1),bands(2*b));
    high = max(bands(2*b-1),bands(2*b));
    targetMask = targetMask | (frequency>=low & frequency<=high);
end
targetMask = targetMask & totalMask;
if ~any(totalMask)
    error('TOPr_Mars:EmptyTotalBand', ...
        'No PMTM frequency bin lies in the total-power interval.');
end
if ~any(targetMask)
    error('TOPr_Mars:EmptyTargetBand', ...
        'No PMTM frequency bin lies in the target-band union.');
end

totalPower = sum(powerByFrequency(totalMask,:),1).';
targetPower = sum(powerByFrequency(targetMask,:),1).';
if any(~isfinite(totalPower) | totalPower<=0)
    error('TOPr_Mars:BadTotalPower', ...
        'A window produced nonpositive or nonfinite total power.');
end

time = linspace(data(1,1)+window/2, ...
    data(end,1)-window/2,windowCount).';
pow = [time,targetPower./totalPower,targetPower,totalPower];
end

function [P,validCount] = row_percentiles(X,percentiles)
P = nan(size(X,1),numel(percentiles));
validCount = sum(isfinite(X),2);
for row = 1:size(X,1)
    values = sort(X(row,isfinite(X(row,:))));
    n = numel(values);
    if n==0
        continue
    end
    for j = 1:numel(percentiles)
        position = 1+(n-1)*percentiles(j)/100;
        lo = floor(position);
        hi = ceil(position);
        if lo==hi
            P(row,j) = values(lo);
        else
            P(row,j) = values(lo)+ ...
                (position-lo)*(values(hi)-values(lo));
        end
    end
end
end

function T = percentile_table(time,values,percentiles)
names = cell(1,numel(percentiles));
for i = 1:numel(percentiles)
    label = strrep(sprintf('p_%g',percentiles(i)),'.','p');
    names{i} = matlab.lang.makeValidName(label);
end
T = array2table(values,'VariableNames',names);
T = addvars(T,time,'Before',1,'NewVariableNames','time_kyr');
end

function names = band_variable_names(periods,suffix)
names = cell(1,numel(periods));
for i = 1:numel(periods)
    periodLabel = strrep(sprintf('%g',periods(i)),'.','p');
    names{i} = matlab.lang.makeValidName( ...
        sprintf('P%s_%s',periodLabel,suffix));
end
end

function plot_uncertainty(time,P,baseline,percentiles,Nmc,plotBaseline)
i05 = find_percentile(percentiles,5);
i25 = find_percentile(percentiles,25);
i50 = find_percentile(percentiles,50);
i75 = find_percentile(percentiles,75);
i95 = find_percentile(percentiles,95);

figure('Color','w','Name','Mars TOPr parameter sensitivity');
hold on;
h90 = fill([time;flipud(time)], ...
    [P(:,i95);flipud(P(:,i05))],[0.82 0.91 0.84], ...
    'LineStyle','none');
h50 = fill([time;flipud(time)], ...
    [P(:,i75);flipud(P(:,i25))],[0.55 0.79 0.61], ...
    'LineStyle','none');
hMedian = plot(time,P(:,i50),'-','Color',[0 0.40 0], ...
    'LineWidth',1.6);
handles = [h90,h50,hMedian];
labels = {'5--95% interval','25--75% interval','MC median'};
if plotBaseline
    hBaseline = plot(time,baseline,':','Color',[0.25 0.25 0.25], ...
        'LineWidth',1.1);
    handles(end+1) = hBaseline;
    labels{end+1} = 'Fixed-parameter reference';
end
hold off;
xlim([time(1),time(end)]);
ylim([0,1]);
xlabel('Time (kyr)');
ylabel('Total orbital power ratio (TOPr)');
title(sprintf('Mars TOPr parameter sensitivity, N=%d',Nmc));
legend(handles,labels,'Location','best');
set(gca,'XMinorTick','on','YMinorTick','on','TickDir','out');
box on;
end

function index = find_percentile(percentiles,target)
index = find(abs(percentiles-target)<=100*eps(max(1,target)),1);
if isempty(index)
    error('TOPr_Mars:MissingPercentile', ...
        'Percentile %g is required.',target);
end
end
