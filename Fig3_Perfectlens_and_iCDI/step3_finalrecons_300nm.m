%{
======== step 3: final reconstruction of cell and lacey ======== 

use refined lacey as input from step 2
Novel part: update lacey for each slice of dynamic sample

%}

clear all
close all
addpath data\
addpath model_and_support_functions_iCDI\

loaddata = 0;
%%
reduce_ratio = 0.3;%300nm
% load(['rec_ov4_complex_',num2str(1000*reduce_ratio),'nm_5to9_0328.mat']); % rec of all 
load(['rec_lacey_300_QE08.mat'])% rec of lacey from step 2

%%
loosesize = 4;% loose support, for high dose, we can reconstruct tight support
smoothratio= 1;
iid=1:20;% frames id used in the GPS; 1 is static; 1:20 is dynamic
params.loosesize = loosesize;
params.smoothratio = smoothratio;

%%
sz=1024; supp_radius=128;
shift_xycell = 124;%round(sz*distanceRatioFromCenter)-50;
shift_xylacey = round(sz*0.23); 

lacey_supp = circshift(makeSquareMask(loosesize+1*supp_radius,sz),[-shift_xylacey -shift_xylacey]);%...
cell_supp = circshift(makeSquareMask(loosesize+1*supp_radius,sz),[+shift_xycell +shift_xycell]);%...

%%
runind=0;indpattern=0;
for runid = 6:9

    runind = runid-5;
    load(['simdp_3.5e',num2str(runid-1),'_300nm.mat']);

    if ~loaddata 
        clear GPS_input;
        GPS_input.support_ref = lacey_supp > 0;
        GPS_input.support_obj = cell_supp > 0;
        GPS_input.probe = 1; 
        GPS_input.iterations = 300;
        GPS_input.sigma=0.0;
        params.iterations = GPS_input.iterations;
        % use reconstructed lacey as initial input for other dose
        rec_laceyi = rec_lacey{runid-5};
        rec_input = smooth3D(rec_laceyi,smoothratio).*lacey_supp;%+cell_supp.* rand(size(initial)) ;%better initial
        GPS_input.thresh = 0.25;
        GPS_input.diffpats = diff_pats(:,:,iid);
        GPS_input.initial = circshift(rec_input,[0,0]);
    
        [recs,probe_new] = ePIE_insitu2(GPS_input);
    
        figure(13); img( recs(:,:,1).*lacey_supp,'', angle(recs(:,:,1).*lacey_supp),'' )
        recs5{runind} = recs;
        params.thresh = GPS_input.thresh;
    else
        % load data
        recs = recs5{runind};
    end
    %% show result of one case
    sz_half = ceil( size(diff_pats,1)/2 ) ;
    x_cen = sz_half+shift_xycell;
    y_cen = sz_half+shift_xycell;
    
    size_crop=250;
    size_crop2=250;
    
    params.size_crop = size_crop;
    params.size_crop2 = size_crop2;
    
    indpattern = indpattern+1;
    for ind = 1:20;
        ind2 = ind;
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
                
        if ind==4;
            figure(51);hold on
            subplot(3,4,indpattern);
            imagesc(model_i);colormap jet;axis image off
            subplot(3,4,indpattern+4);
            imagesc(perfecLens_i);title('PL');colormap jet;axis image off
            subplot(3,4,indpattern+8);
            imagesc(rec_i_shift);title('CDI');colormap jet;axis image off
        end
        % FRC
        nn=20;
        [correlation_rec_i(:,ind), freq]  = FourierShellCorrelate(model_i, rec_i_shift,nn);
        [correlation_PL_i(:,ind) , freq]  = FourierShellCorrelate(model_i, perfecLens_i,nn);
    end
    
    figure(70); hold on
    correlation_PL_im =  mean(correlation_PL_i,2);
    correlation_rec_im =  mean(correlation_rec_i,2);
    plot(freq(1:end),correlation_PL_im(1:end),'b-',...
        freq(1:end),mean(correlation_rec_im(1:end),2),'r-', 'LineWidth',1.5 );   
end
hold on
plot(freq,exp(-1).*ones(nn,1),'k-', 'LineWidth',1 )
legend('PL 3.5e5', 'iCDI 3.5e5','PL 3.5e6', 'iCDI 3.5e6',...
    'PL 3.5e7', 'iCDI 3.5e7','PL 3.5e8', 'iCDI 3.5e8', '1/e');
xlim([0,1]);ylim([0,1.01]);
title([num2str(1000*reduce_ratio),'nm'])

%%
%{

p.e.FOV=FOV;
p.e.loosesize = loosesize;
p.e.smoothratio = smoothratio;
save(['rec_ov4_complex_',num2str(reduce_ratio*1000),'nm_5to9_0606.mat'],'recs5','p','params','-v7.3')

%}
