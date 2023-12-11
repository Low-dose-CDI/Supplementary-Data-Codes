%{
======== step 3: test FRC ======== 

% change line 24 runid=1:4 (including ptychography and Lo-CDI)
% if figures flip, change size in line 17

%}

clear all
close all

addpath model_and_support_functions_iCDI\
addpath dataoutput\

%% parameters
sz=1024;
size_crop =274;
size_crop2 = 240;
shift_xycell = 184;
shift_xylacey = 256;
indpattern = 0;

%%
for runid = 1:4;

    load(['ptycho10x10_300_ref',num2str(0),'id',num2str(runid),'.mat']);
    big_obj_crop1 = epieoutput.object;
    aperture_crop1 = epieoutput.probe;
    load(['ptycho10x10_300_ref',num2str(1),'id',num2str(runid+5),'.mat']);
    big_obj_crop2 = epieoutput.object;

    %% FRC
    sz_half = ceil( size(aperture_crop1,1)/2 ) ;
    x_cen = sz_half+shift_xycell;
    y_cen = sz_half+shift_xycell;
    
    model_truth_framesi = abs(epieoutput.model_truth_framesi);
    per_image_noisei = abs(epieoutput.per_image_noisei);
    
    model_i0 = crop_roi((model_truth_framesi), size_crop, y_cen,x_cen);
    model_i0 = model_i0/mean(model_i0(:));
    model_i = crop_roi(model_i0,size_crop2);
    minv = min(model_i(:));maxv=max(model_i(:));

    recsi1 = abs(crop_roi(big_obj_crop1,[sz sz]));
    rec_i1       = crop_roi((recsi1)/max((recsi1(:))), size_crop, y_cen,x_cen);
    rec_i_shift1 = align2D( model_i0, rec_i1 );
    rec_i_shift1 = crop_roi(rec_i_shift1,size_crop2);
    rec_i_shift1 = rescale((rec_i_shift1),minv,maxv);
    
    recsi2 = abs(crop_roi(big_obj_crop2,[sz sz]));
    rec_i2       = crop_roi((recsi2)/max((recsi2(:))), size_crop, y_cen,x_cen);
    rec_i_shift2 = align2D( model_i0, rec_i2 );
    rec_i_shift2 = crop_roi(rec_i_shift2,size_crop2);
    rec_i_shift2 = rescale((rec_i_shift2),minv,maxv);
    
    perfecLens_i = crop_roi(per_image_noisei,size_crop, y_cen,x_cen);
    perfecLens_i = crop_roi(perfecLens_i,size_crop2);
    perfecLens_i = rescale((perfecLens_i),minv,maxv);
  
    indpattern=indpattern+1;

    figure(52);
    hold on
    subplot(3,4,indpattern);imagesc(perfecLens_i);colormap jet;axis equal off;%,perfecLens_i,'')
    subplot(3,4,indpattern+4);imagesc(rec_i_shift1);colormap jet;axis  equal off;%,perfecLens_i,'')
    subplot(3,4,indpattern+8);imagesc(rec_i_shift2);colormap jet;axis  equal off;%,perfecLens_i,'')

    [correlation_PL_i0 , freq]  = FourierShellCorrelate(model_i, perfecLens_i,21);
    [correlation_rec_i1, freq]  = FourierShellCorrelate(model_i, rec_i_shift1,21);
    [correlation_rec_i2, freq]  = FourierShellCorrelate(model_i, rec_i_shift2,21);
    figure(1);ylim([-0.1,1.01]);hold on
    plot( freq,correlation_PL_i0,'b--',freq,correlation_rec_i2, 'r-','LineWidth',1 );hold on;
    figure(2);ylim([-0.1,1.01]);hold on
    plot( freq,correlation_rec_i1,'b--',freq,correlation_rec_i2, 'r-','LineWidth',1 );hold on;
   
%     figure(71); img(perfecLens_i,'PL', rec_i_shift1,'PTY',rec_i_shift2,'INSITU PTY','abs','off','colormap','jet')
    pause(1)
        
end
figure(2);hold on
xlim([0.0,1]);ylim([0,1.01]);
hold on;plot(freq,exp(-1).*ones(21,1),'k--', 'LineWidth',1 )
legend('Ptycho 3.5E5','LoCDI 3.5E5','Ptycho 3.5E6','LoCDI 3.5E6','Ptycho 3.5E7','LoCDI 3.5E7','Ptycho 3.5E8','LoCDI 3.5E8','1/e')
figure(1);hold on
xlim([0.0,1]);ylim([0,1.01]);
hold on;plot(freq,exp(-1).*ones(21,1),'k--', 'LineWidth',1 )
legend('PL 3.5E5','LoCDI 3.5E5','PL 3.5E6','LoCDI 3.5E6','PL 3.5E7','LoCDI 3.5E7','PL 3.5E8','LoCDI 3.5E8','1/e')
