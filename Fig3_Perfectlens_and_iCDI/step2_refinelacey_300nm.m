%{
======== step 2: lacey refinement ======== 

% 1: try 3-4 times HIO reconstruction with the 1st dp of 20 
% do "circshift(u,[y,x])" or "-conj(recs)" within the FOV of loose support for each round
% 2: GPS shrinkwrap
% 3: test FRC
% repeat step1-3 to refine lacey reconstruction

%}
 
clear all
close all
addpath data\
addpath model_and_support_functions_iCDI\
reduce_ratio = 0.3;%300nm

runid = 7;

%% data loading
load(['simdp_3.5e',num2str(runid-1),'_',num2str(reduce_ratio*1000),'nm.mat']);

%% loose support 
loosesize = 4;% loose support, for high dose, we can reconstruct tight support
smoothratio= 0.5;
sz=1024; 
supp_radius=128;shift_xycell = 124;shift_xylacey = round(sz*0.23);
distanceRatioFromCenter=0.23;sz=1024;
shift_xy=164;

supp_lacey = circshift(makeSquareMask(loosesize+supp_radius,sz),...
    [-shift_xylacey -shift_xylacey]);%...
supp_cell = circshift(makeSquareMask(loosesize+supp_radius,sz), [shift_xycell shift_xycell] );

supp_lacey = supp_lacey~=0;
supp_cell = supp_cell~=0;

%% round1: try multiple times HIO reconstruction with 1 slice of 20
support = supp_lacey | supp_cell;
non_support_ref = ~supp_lacey;
no_supp = ~support;
dp = diff_pats(:,:,1);
[N1,N2] = size(dp);
stopper = dp==-1;
norm_dp = norm ( dp(~stopper), 'fro' );

iterations=10000;
u = 0.1*randn(N1) .*supp_lacey;

%% %%%%%%%%%%%%%%%%%%%%%% take 2 or more times step1-step 3 %%%%%%%%%%%%%%%%%%
% !!! after step 3 FRC check, back here to choose one of them as initial of next round HIO
%{
u = (recs);
u = -conj(recs);
%}

%% step 1: few rounds HIO reconstruction: shift reconstruction and try again

% flip the image if needed
%{
lac = u(supp_lacey);
lac = lac(end:-1:1);
u(supp_lacey) = lac;
%}


% !!!!!! [0,0] for 1st round
% move around to remove artifact (by monitoring error
shiftu_right = -8; % 0 for 1st round, then 8/-8/0 as you like
shiftu_down = -8; % 0 for 1st round, then 8/-8/0 as you like


u = circshift(u,[shiftu_down,shiftu_right]); 
thresh=0;
err_best=1;
ker = fspecial('disk', 5);

for t=1:iterations
    Fu = fftshift(fft2(u));        
    z_stopper = Fu(stopper);

    Fu = dp .* exp(1i*angle(Fu));
    Fu(stopper) = z_stopper;

    u_new = ifft2(ifftshift(Fu));
    u_new(no_supp) = u(no_supp) - 1*u_new(no_supp);
    u = u_new;

    if t>0
        u_abs = abs(u);
        u_abs = conv2(u_abs, ker, 'same');
        support_refined = u_abs.*supp_lacey;

        support_refined = support_refined/max(support_refined(:)); % make in range [0,1]
        support_refined = support_refined > thresh;
        support_refined1 = support_refined | non_support_ref;

        u = u.*support_refined1;
        if mod(t,50000)==0
            figure(201); img(support_refined1,'support'); drawnow;
        end
        support = support_refined | supp_cell;
        non_support_ref = ~support_refined;
        no_supp = ~support;
    end

    % error
    Fu = fftshift(fft2(u.*support));
    diff_z = abs(Fu(~stopper)) - dp(~stopper);
    err = norm(diff_z, 'fro') / norm_dp;
    if err<err_best
        u_best = u;
        err_best = err;
    end

    if mod(t,20)==0
        figure(4);
        rec = u_new.*support;
        img( rec, ['mag ',num2str(t)], 'colormap','jet')
        drawnow
        
        fprintf('%d. error = %f\n',t,err);
%         if err>err_best+0.001;break;end

    end
end
% figure(20); img( real(recs.*supp_cell),'real', imag(recs.*supp_cell),'imag','colormap','gray')
pause(1)
%% step 2: GPS_shrinkwrap reconstruction
clear GPS_input;

GPS_input.initial = circshift(u,[4,4]); % shift HIO output to the center

GPS_input.support_ref = supp_lacey;
GPS_input.support_obj = supp_cell;
GPS_input.probe = 1; 
GPS_input.iterations = 1000;
GPS_input.sigma=0.0;
GPS_input.thresh=0.12;
GPS_input.diffpat = diff_pats(:,:,1);

[recs,probe_new] = GPS_shrinkwrap2(GPS_input);
figure(33);img(real(recs.*supp_cell))

%%  step 3: check FRC and back to step 1 if necessary
% Fourier shell correlation: show result of one case
sz_half = ceil( size(diff_pats,1)/2 ) ;
x_cen = sz_half+shift_xycell;
y_cen = sz_half+shift_xycell;

size_crop=250;
size_crop2=250;

params.size_crop = size_crop;
params.size_crop2 = size_crop2;

ind = 1;ind2 = ind;
smoothmask = 1;%smooth3D(makeSquareMask(FOVradius,size_crop),smoothratioFSC);

recsi = real(recs(:,:,ind));
model_truth_framesi = real(model_truth_frames{ind2});
per_image_noisei = real(per_image_noise{ind2});

model_i = crop_roi((model_truth_framesi), size_crop, y_cen,x_cen);
model_i = model_i/mean(model_i(:));
minv = min(model_i(:));maxv=max(model_i(:));
rec_i        = crop_roi((recsi)/max((recsi(:))), size_crop, y_cen,x_cen); 
rec_i_shift = align2D( model_i, rescale((rec_i),minv,maxv) );

perfecLens_i = crop_roi((per_image_noisei),size_crop, y_cen,x_cen);  
perfecLens_i = rescale((perfecLens_i),minv,maxv);

model_i        = croppedOut(model_i, size_crop2); 
rec_i_shift        = croppedOut(rec_i_shift, size_crop2); 
perfecLens_i        = croppedOut(perfecLens_i, size_crop2); 

nn=20;
[correlation_rec_i(:,ind), freq]  = FourierShellCorrelate(model_i, rec_i_shift,nn);
[correlation_PL_i(:,ind) , freq]  = FourierShellCorrelate(model_i, perfecLens_i,nn);

figure(70); hold on
correlation_PL_im =  mean(correlation_PL_i,2);
correlation_rec_im =  mean(correlation_rec_i,2);
plot(freq(1:end),correlation_PL_im(1:end),'b-',...
freq(1:end),mean(correlation_rec_im(1:end),2),'r-', 'LineWidth',1.5 );   
pause(1)
hold on
plot(freq,exp(-1).*ones(20,1),'k-', 'LineWidth',1 )

xlim([0,1]);
ylim([0,1.01]);

%% save lacey
%{
filepath='./data/';
load('rec_lacey_300_QE08.mat');
rec_lacey{runid-5}=recs.*supp_lacey;save([filepath,'rec_lacey_300_QE08'],'rec_lacey')
%}
