

function points = locate_cell_bodies_NaJi(file_name,options)

% Finds 
% 1. Transform image (sqrt) 
% 2. High-pass filter for each frame
% 3. Dynamic thresholding for each frame (select threshold with greatest number of components after strel cleanup)
% 4. Morphological opening on each frame
% 5. Convolve thresholded image with 3D gaussian -> 3D Watershed
% 6. Peaks of 3D watershed --> Identified points.
% 
% Plausible improvements include local thresholding and 
% 
% 
% Alex Vaughan, 2017

% Load tsStack
if nargin == 0
    
    %% Na Ji's data.
    [file_dir,~] = fileparts(which(mfilename));

    %file = '/data/Dataset1/20161017_11_L30_D200_300pixelms_x1um_Gaussian3d.tif';
    %file = '/data/Dataset2/20160211_03_D150_L10_X0.8um_gaussian3D.tif';
    file = '/data/Dataset3/20160212_02_D250_L20_X1um_Gaussian3D.tif';
    
    file_name = fullfile(file_dir,file);
    
    options.flip_dimensions = '';
    %Paramters - for now all are required.
    options.frame_channel_order     = 'cf'; % outer channels,inner frames
    options.dapi       = 1;
    options.dii        = NaN;
    options.gcamp      = NaN;
    options.n_channels = 1;
    options.do_plots   = 1;
    % Variance scaling - sqrt original image?
    options.do_sqrt = 1;
    % High-pass filter
    options.sigma_hp   = 50; % Background subtraction
    % Strel size - sets size of discovered components
    options.strel_size = 5; 
    % Gaussian convolution size - also sets size of discovered components
    options.gaussian_sigma = 5;

    % Range of dynamic thresholds (constant for each frame)
    options.thresholds = 0.1:0.05:0.5;
    % Drop a percentage of watershed outliers
    options.watershed_outlier_fraction = 0.01;

end

%% Default options

do_plots   = 1;

%Image paramters - for now all are required.
flip_dimensions = '';
frame_channel_order     = 'cf'; % outer channels,inner frames
dapi       = 1;
n_channels = 1;
z_scale    = 1;

% Variance scaling - sqrt original image?
do_sqrt = 1;
% High-pass filter
sigma_hp   = 50; % Background subtraction
% Strel size - sets approximate size of components
strel_size = 4;
% Gaussian convolution size - also sets size of discovered components
options.gaussian_sigma = 4;

% Range of dynamic thresholds (constant for each frame)
thresholds = 0.1:0.05:0.5;
% Drop a percentage of watershed outliers
watershed_outlier_fraction = 0.01;

if exist('options','var')
    % Pass options into main inputs.
    input_fields = fieldnames(options);
    for f = 1:length(input_fields)
        eval(sprintf('%s = options.%s;',input_fields{f},input_fields{f}));
    end
end

%% Load data with TiffStack 

% Adjust tiff /montage warnings
warning off MATLAB:imagesci:tiffmexutils:libtiffWarning
warning off MATLAB:imagesci:tifftagsread:expectedTagDataFormat
warning off MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning
warning off images:imshow:magnificationMustBeFitForDockedFigure

% Load data as tsStack using TIFFStack library
tsStack  = TIFFStack(file_name);
[xpix,ypix,n_pseudoframes] = size(tsStack);
n_frames = n_pseudoframes / n_channels;
assert(rem(n_frames ,1) == 0,'# of frames or channels is wrong')

% Helper: frame lookup accomodates flipping a slice by reversing order in z
% Lookup frames using tsFrames(x,y,z)
flipx = ~isempty(strfind(flip_dimensions,'x'));
flipy = ~isempty(strfind(flip_dimensions,'y'));
flipz = ~isempty(strfind(flip_dimensions,'z'));
if flipz
    % Inverted z lookup
    switch frame_channel_order
        case 'fc'
            frame_lookup = @(frame,channel) 3*((n_frames-frame+1)-1)+channel;
        case 'cf'
            frame_lookup = @(frame,channel) z_scale*(channel-1)+n_frames-frame+1;
    end

else
    % Normal z lookup
    switch frame_channel_order
        case 'fc'
            frame_lookup = @(frame,channel) 3*(frame-1)+channel;
        case 'cf'
            frame_lookup = @(frame,channel) z_scale*(channel-1)+frame;
    end
end
if flipx && flipy
    tsFrames = @(frame,channel) double(tsStack(end:-1:1,end:-1:1,frame_lookup(frame,channel)));
elseif flipy
    tsFrames = @(frame,channel) double(tsStack(end:-1:1,:,frame_lookup(frame,channel)));
elseif flipx
    tsFrames = @(frame,channel) double(tsStack(:,end:-1:1,frame_lookup(frame,channel)));
else
    tsFrames = @(frame,channel) double(tsStack(:,:,frame_lookup(frame,channel)));
end

% Helper : namedFigure with black background.  isempty() is just so that it passes a value.
namedFigure_k = @(xpix) isempty(namedFigure(xpix)) & set(gcf, 'color', [0 0 0]); % Black-background named figure

% Helper : vline with quantile of data, labeled.
labeled_quantile = @(data,q) set(vline(quantile(data,q),'k--'),'LineWidth',2) & isempty(text(quantile(data,q)*1.1,max(ylim)*0.8,sprintf('%.0f%%',q*100)));


%% Sanity check figure.

% namedFigure('Sanity check - frames.');
% for c = 1:n_channels
%     subplot(3,1,c)
%     imagesc(tsFrames(10,c));
%     axis image
%     frame_lookup(10,c)
% end


%% BEGIN ALGORITHM %%
% Main processing pipeline.
% 1. Frame-wise low-pass filter to equalize local variance.
% 2. find best threshold (maximize components that survive morphological opening)
% 3. Using best threshold, perform 2d + 3d morphological opening
% 4. Trim components based on area (toss biggest/smallest 5%)
% 5. Done.

[stack_raw,stack_filtered,stack_threshold] = deal(zeros(xpix,ypix,n_frames));
[best_threshold_f,best_ind_f,n_components_f] = deal(zeros(n_frames,1));
for i = 1:n_frames
    % Sqrt transform original stack?
    if do_sqrt
        stack_raw(:,:,i) = sqrt(tsFrames(i,dapi));
    else
        stack_raw(:,:,i) = tsFrames(i,dapi);
    end
    % High-pass filter.
    stack_filtered(:,:,i) = stack_raw(:,:,i) - imgaussfilt(stack_raw(:,:,i),sigma_hp);
    % Find best threshold for filtered images
    [best_threshold_f(i),best_ind_f(i),n_components_f(i)] = find_threshold(stack_filtered(:,:,i,dapi),thresholds,strel_size);
    fprintf('Slice %.0f :: best threshold_ind %.0f :: best threshold %.2f :: %.0f\n',i,best_ind_f(i),best_threshold_f(i),n_components_f(i))
    % Filter by best threshold
    stack_threshold(:,:,i) = stack_filtered(:,:,i) >= quantile(ravel(stack_filtered(:,:,i)),best_threshold_f(i));
    % Open image
    stack_threshold(:,:,i) = imopen(stack_threshold(:,:,i),strel('disk',strel_size));
end

    
%% Watershed

fprintf('')
% Make a guassian psf for convolution for watershedding
n = 2*gaussian_sigma; 
nset = -n:n; 
gpsf = zeros(2*n+1,2*n+1,2*n+1);
for i = 1:2*n+1, for j = 1:2*n+1, for k = 1:2*n+1
    gpsf(i,j,k) =  1/(gaussian_sigma^3*sqrt((2*pi)^3)) * exp( - (nset(i)^2 + nset(j)^2 + nset(k)^2 )/(gaussian_sigma^2));
end, end, end

% Convolve
stack_convolved = convn(stack_threshold,gpsf,'same');
stack_watershed = watershed(-stack_convolved); 
components = regionprops(stack_watershed,'Area', 'BoundingBox', 'Centroid', 'FilledArea', 'FilledImage', 'Image','PixelIdxList','PixelList')


%% Drop outliers.
components = components([components.Area] > quantile([components.Area],watershed_outlier_fraction/2));
components = components([components.Area] < quantile([components.Area],1-watershed_outlier_fraction/2));

%% Clean up components
lower_quantile = quantile([components.Area],0.05);
upper_quantile = quantile([components.Area],0.95);
components = components([components.Area] > lower_quantile);
components = components([components.Area] < upper_quantile);
watershed_centroids    = reshape([components.Centroid]',3,length(components))';
assert(all(watershed_centroids(end,:) == components(end).Centroid),'Centroid matrix reshaping is incorrect.')

%% Find peak in each centroid.

watershed_peaks = zeros(size(components,1),3);
for ww = 1:size(components,1)
   % This is a dumb way to index; there has got to be a better way
   this_pix = [components(ww).PixelList];
   this_pix_ind = sub2ind(size(stack_watershed),this_pix(:,2),this_pix(:,1),this_pix(:,3)); 
   [~,max_ind] = max(stack_convolved(this_pix_ind));
   watershed_peaks(ww,:) = this_pix(max_ind,:);
end

% namedFigure('Peaks vs. centroids')
% subplot(1,3,1); plot(watershed_centroids(:,1),watershed_peaks(:,1),'.')
% subplot(1,3,2); plot(watershed_centroids(:,2),watershed_peaks(:,2),'.')
% subplot(1,3,3); plot(watershed_centroids(:,3),watershed_peaks(:,3),'.')

% Outputs
points = watershed_peaks;


%%%%%% DONE WITH ALGORITHM %%%%%


%%%%%% DONE WITH ALGORITHM %%%%%

fprintf('Found %.0f centroids,',length(points))


if ~do_plots
    return
end



%% %%%% BEGIN PLOTTING %%%%%%%%%%

%% 
namedFigure('Watershed outlines')
f = 5 + ceil(40 * rand);
%w = watershed(-stack_convolved(:,:,f)); 
w = stack_watershed(:,:,f);
outline_watersheds = @(mat)  max( max(ravel(double(mat)))*double(w==0),double(mat));
% Plot centroids in a range around the slice of interest.
z_range = 10;
centroid_ids = find( watershed_peaks(:,3) < f+z_range &...
                     watershed_peaks(:,3) > f-z_range);
ax = [];
ax(end+1) = subplot_tight(3,2,1); hold on
imagesc(stack_raw(:,:,f))
scatter(watershed_peaks(centroid_ids,1),watershed_peaks(centroid_ids,2),10,'r','filled')
%imagesc(w); 
axis xy; axis image; %colorbar
title(sprintf('Original DAPI image (%.0f)',f))

ax(end+1) = subplot_tight(3,2,2);
imagesc((stack_filtered(:,:,f)))
%imagesc(w); 
axis image; axis xy;  %colorbar
title('Filtered')

ax(end+1) = subplot_tight(3,2,3);
imagesc((stack_threshold(:,:,f)))
%imagesc(w); 
axis xy; axis image; %colorbar
title('Thresholded')


ax(end+1) = subplot_tight(3,2,4);
imagesc(stack_convolved(:,:,f))
%image( tsStack(f,dapi))
axis xy;axis image; %colorbar
title('Convolved')

                 
ax(end+1) = subplot_tight(3,2,5); cla; hold on
imagesc(outline_watersheds(stack_convolved(:,:,f)))
%imagesc((w~=0) .* (conv_out(:,:,f)));
scatter(watershed_peaks(centroid_ids,1),watershed_peaks(centroid_ids,2),10,'r','filled')
axis xy; axis image
title('Watershed detection')

ax(end+1) = subplot_tight(3,2,6); cla; hold on;
imagesc(outline_watersheds(stack_filtered(:,:,f)))
scatter(watershed_peaks(centroid_ids,1),watershed_peaks(centroid_ids,2),10,'r','filled')
axis xy; axis image;
%drawnow
linkaxes(ax)
title('Alignment with DAPI')



%% Histogram of component areas
namedFigure('Component areas');  clf
hx = logspace(log10(min([components.Area])),log10(max([components.Area])),100);
bar(hx,histc([components.Area],hx))

set(gca,'xscale','log'); 
set(gca,'yscale','log')
xlabel('Component Volume')
ylabel('Frequency')
labeled_quantile([components.Area],0.5)
labeled_quantile([components.Area],0.9)

set(vline(4/3*pi*9^3,'r'),'LineWidth',2)
text(4/3*0.65*pi*9^3,850,{'9µm diameter','50% packing'},'Color','r')

set(gca,'XTick',[1000 2000 5000 10000])
axis tight


%% Plot all cell bodies
namedFigure('Identified components'); clf; 
hold on
scatter3(watershed_peaks(:,1),watershed_peaks(:,2),watershed_peaks(:,3),10,'r','filled')
imagesc([0,xpix],[0,ypix],sqrt(sum(tsFrames(1:n_frames,dapi),3)))
axis image;  axis xy
title('Sum stack vs. identified centers')



%% Plot overlay with points in a small volume

z_range = 10;  % 10µm above and below.
f = z_range + ceil(rand*(n_frames-2*z_range)); % Random frame
namedFigure_k('DAPI overlay( small volume )'); clf

ax = subplot(1,2,1); hold on
imagesc([0 xpix],[0 ypix],sum(stack_raw(:,:,f-z_range:f+z_range),3))
centroid_ids = find(watershed_peaks(:,3) < f+z_range & watershed_peaks(:,3) > f-z_range);
scatter3(watershed_peaks(centroid_ids,1),watershed_peaks(centroid_ids,2),watershed_peaks(centroid_ids,3) - min(watershed_peaks(centroid_ids,3)),10,'r','filled')
axis tight; axis xy; axis image; view(2)
title(sprintf('DAPI cell bodies :: raw image :: frame %.0f',f),'Color','w')

zscore_fn = @(x) (x - mean(x(:)))/var(x(:));
ax(2) = subplot(1,2,2); hold on
imagesc([0 xpix],[0 ypix],zscore_fn(sum(stack_raw(:,:,f-z_range:f+z_range),3)))
scatter3(watershed_peaks(centroid_ids,1),watershed_peaks(centroid_ids,2),watershed_peaks(centroid_ids,3) - min(watershed_peaks(centroid_ids,3)),10,'r','filled')
axis tight;  axis image; axis xy;
axis equal; view(2)
title(sprintf('DAPI cell bodies :: Z-scored image :: frame %.0f',f),'Color','w')
linkaxes(ax)

%% Montage - linked axes
ns = 5;
frames = ceil(linspace(1,n_frames,ns^2));
all_axes = [];
z_range = 5;
namedFigure_k('Dapi Montage'); clf
for i = 1:ns^2
    all_axes(end+1) = subplot_tight(ns,ns,i,[0.01,0.01]); hold on
    imagesc(tsFrames(frames(i),dapi))
    centroid_ids = find(watershed_peaks(:,3) < frames(i)+z_range*2 & watershed_peaks(:,3) > frames(i)-z_range*2);
    scatter3(watershed_peaks(centroid_ids,1),watershed_peaks(centroid_ids,2),watershed_peaks(centroid_ids,3) - min(watershed_peaks(centroid_ids,3)),10,'r','filled')
    axis off; axis xy; axis equal
    text(10,100,num2str(frames(i)),'Color','r')
end
namedFigure_k('Dapi Threshold'); clf
for i = 1:ns^2
    all_axes(end+1) = subplot_tight(ns,ns,i,[0.01,0.01]); cla; hold on
    imagesc(stack_threshold(:,:,frames(i)))
        centroid_ids = find(watershed_peaks(:,3) < frames(i)+z_range*2 & watershed_peaks(:,3) > frames(i)-z_range*2);
    scatter3(watershed_peaks(centroid_ids,1),watershed_peaks(centroid_ids,2),watershed_peaks(centroid_ids,3) - min(watershed_peaks(centroid_ids,3)),10,'r','filled')

    axis off; axis xy; axis image; view(2)
    text(10,100,num2str(frames(i)),'Color','r')
end
linkaxes(all_axes)

%% Histogram of DAPI values.
namedFigure('Dapi Intensities + Best Thresholds')
ns = 3;
frames = ceil(linspace(1,n_frames,ns^2));
%hx = linspace(0,4096,100);
for i = 1:ns^2
    subplot_tight(ns,ns,i,[0.05,0.05])
    %bar(hx,histc(ravel(tsFrames(frames(i),dapi)),hx))
    hist(ravel(tsFrames(frames(i),dapi)),100)%,'Normalization','pdf')
    labeled_quantile(ravel(tsFrames(frames(i),dapi)),best_threshold_f(frames(i)));
    %set(gca,'yscale','log'); 
   % set(gca,'xscale','log'); 
    title(sprintf('Frame %.0f',frames(i)))
end


keyboard

function [best_threshold,best_ind,n_components] = find_threshold(img,thresholds,strel_size)
% Given a scalar image 'img', find the binarization threshold such that we 
% recover the maximum number of components after cleanup with a disk-shaped strel 
% of size strel_size.

for t = 1:length(thresholds)
    CC(t) = bwconncomp(imopen(img > quantile(img(:),thresholds(t)),strel('disk',strel_size)));
end

[n_components,best_ind] = max([CC.NumObjects]);
best_threshold = thresholds(best_ind);
