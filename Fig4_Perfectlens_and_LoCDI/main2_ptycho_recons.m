clear all
close all

addpath model_and_support_functions_iCDI\

%% could be loose
loosesize=4;%% loose support
zz = 12;% real distance is 8um
supp_radius=106/2;
scanFOV=256/2;
shift_xycell = 184;
shift_xylacey = 256;

%% 
sz=1024;NN=1290;

for runid = [2,4]; 
    if runid<=5;ref=0;else;ref=1;end
    
    load(['ptycho10x10_300dp_ref',num2str(ref),'_3.5e',num2str(5-ref*5+runid-1),'.mat']);

    %%%  rpie    
    probe_cell0 = circshift(makeCircleMask(-loosesize+1*supp_radius,sz),[+shift_xycell +shift_xycell]);%...
    cell_supp_p = circshift(makeCircleMask(0+1*supp_radius,sz),[+shift_xycell +shift_xycell]);%...
    cell_supp_O = circshift(makeSquareMask(10+scanFOV,sz),[+shift_xycell +shift_xycell]);%...
    support_O = My_paddzero( cell_supp_O ,[NN,NN]);
    
    clear inputs
    
    if ref==1
        probe_lacey0 = circshift(makeCircleMask(-loosesize+1*supp_radius,sz),[-shift_xylacey -shift_xylacey]);%...
        lacey_supp_p = circshift(makeCircleMask(30+1*supp_radius,sz),[-shift_xylacey -shift_xylacey]);%...
        lacey_supp_O = circshift(makeSquareMask(10+scanFOV,sz),[-shift_xylacey -shift_xylacey]);%...
        inputs.InitialObj = My_paddzero(lacey_supp_O+0,[NN,NN]);
        inputs.InitialAp  = fresnel_advance(probe_cell0./(10^(2))+probe_lacey0,p.e.px_size,p.e.px_size,zz,p.e.lambda);
    else
        inputs.InitialObj = rand(NN).*support_O;
        inputs.InitialAp  = fresnel_advance(probe_cell0,p.e.px_size,p.e.px_size,zz,p.e.lambda);
    end
    inputs.lacey_supp = 0;
    %%%
    inputs.Iterations = 200;
    inputs.FileName  = 'iCDI';
    inputs.Patterns  = dp_noisy;
    inputs.Positions = [centery' centerx'];
    inputs.showim     = 1;
    inputs.ApRadius  = -1;
    inputs.do_posi = 0;
    inputs.PixelSize = 1;
    inputs.GpuFlag=1;
    inputs.bin = 1;
    inputs.obj_scale = 1;
    inputs.updateAp  = 1;
    beta_obj = 0.1; beta_ap = 0.1;
    ratio = 1;
    
    for round0 = 1:4  %:6
   
%         ePIE
%         if round0>2; beta_obj = 0.1; beta_ap = 0.5;inputs.Iterations = 200;end
%         [big_obj_crop1,aperture_crop1,fourier_error,~,~] = ePIE(inputs,beta_obj,beta_ap);
%          rPIE
        if round0>2; beta_obj = 0.9; beta_ap = 0.5;inputs.Iterations = 200;end
        [big_obj_crop1,aperture_crop1,fourier_error,~,~] = rPIE(inputs,beta_obj,beta_ap);
    
        %%% prepare for next epie
        if ref==1
            support_refined = big_obj_crop1;
            support_refined1 = support_refined/max(support_refined(:)); % make in range [0,1]
            inputs.InitialObj  = (support_refined1)>0.6;
            inputs.InitialAp  = aperture_crop1.*(lacey_supp_p+cell_supp_p);
        else
            inputs.InitialObj  = support_O.*(big_obj_crop1) ;
            inputs.InitialAp  = aperture_crop1.*(cell_supp_p);
        end
        figure(22);img(big_obj_crop1,'',aperture_crop1,'',inputs.InitialObj,'',inputs.InitialAp,'')

        %% FRC
        size_crop =268;
        size_crop2 = 250;
        sz_half = ceil( size(aperture_crop1,1)/2 ) ;
        x_cen = sz_half+shift_xycell;
        y_cen = sz_half+shift_xycell;
        
        model_truth_framesi = real(modeltruth );
        per_image_noisei = real(per_image_noise );
        
        model_i0 = crop_roi((model_truth_framesi), size_crop, y_cen,x_cen);
        model_i0 = model_i0/mean(model_i0(:));
        model_i = crop_roi(model_i0,size_crop2);
        minv = min(model_i(:));maxv=max(model_i(:));
    
        recsi1 = real(crop_roi(big_obj_crop1,[sz sz]));
        rec_i1       = crop_roi((recsi1)/max((recsi1(:))), size_crop, y_cen,x_cen);
        rec_i_shift1 = align2D( model_i0, rec_i1 );
        rec_i_shift1 = crop_roi(rec_i_shift1,size_crop2);
        rec_i_shift1 = rescale((rec_i_shift1),minv,maxv);
        
        perfecLens_i = crop_roi(per_image_noisei,size_crop, y_cen,x_cen);
        perfecLens_i = crop_roi(perfecLens_i,size_crop2);
        perfecLens_i = rescale((perfecLens_i),minv,maxv);
          
        figure(72);ylim([-0.1,1.01]);hold on
        [correlation_PL_i0 , freq]  = FourierShellCorrelate(model_i, perfecLens_i,21);
        [correlation_rec_i1, freq]  = FourierShellCorrelate(model_i, rec_i_shift1,21);

        plot( freq,correlation_PL_i0,'b--',freq,correlation_rec_i1, 'r-','LineWidth',1 )
        xlim([0.0,1]);ylim([0,1.01]);hold on;
        plot(freq,exp(-1).*ones(21,1),'k--', 'LineWidth',1 )
        
        figure(71); img(model_i,'',perfecLens_i,'PL', rec_i_shift1,'INSITU PTY','abs','off','colormap','jet')
        pause(1)
        %%
    end
    
    epieoutput.object=big_obj_crop1;
    epieoutput.probe=aperture_crop1;
    epieoutput.model_truth_framesi= (modeltruth );
    epieoutput.per_image_noisei= (per_image_noise );
    save(['ptycho10x10_300_ref',num2str(ref),'id',num2str(runid),'.mat'],'epieoutput','p')
    
end