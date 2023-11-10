
function [big_obj,aperture,fourier_error,initial_obj,initial_aperture] = rPIE(ePIE_inputs,varargin)
%varargin = {beta_obj, beta_ap}
optional_args = {0.05 .1}; %default values for varargin parameters
rng('shuffle','twister');
%% setup working and save directories

[~,jobID] = system('echo $JOB_ID');
jobID = jobID(~isspace(jobID));
if isempty(jobID)
    time_str = char(datetime('now','Format','uuuuMMdd_HHmm'));
    jobID = sprintf('local_%s',time_str);
end

%% Load inputs from struct
diffpats = ePIE_inputs(1).Patterns;
positions = ePIE_inputs(1).Positions;
% filename = ePIE_inputs(1).FileName;
% filename = strcat(filename, 'ePIE_mono_roundShiftObj');
pixel_size = ePIE_inputs(1).PixelSize;
big_obj = ePIE_inputs(1).InitialObj;
aperture_radius = ePIE_inputs(1).ApRadius;
aperture = ePIE_inputs(1).InitialAp;
iterations = ePIE_inputs(1).Iterations;
filename = ePIE_inputs(1).FileName;
filename = strcat('rPIE_',filename,'_',jobID);
filename = strrep(filename,'__','_');
%% parameter inputs
if isfield(ePIE_inputs, 'updateAp')
    update_aperture = ePIE_inputs.updateAp;
else
    update_aperture = 1;
end

lacey_supp = ePIE_inputs.lacey_supp;


if isfield(ePIE_inputs, 'GpuFlag')
    gpu = ePIE_inputs.GpuFlag;
else
    gpu = 0;
end

probe_mask_flag = 0;
if isfield(ePIE_inputs, 'R_probe_mask') %real space probe mask
    R_probe_mask = single(ePIE_inputs.R_probe_mask);
    %params.R_probe_mask = R_probe_mask;
    probe_mask_flag = 1;
end

if isfield(ePIE_inputs, 'F_probe_mask') %k-space probe mask (don't fftshift)
    F_probe_mask = single(ePIE_inputs.F_probe_mask);
    %params.F_probe_mask = F_probe_mask;
    probe_mask_flag = 2;
end

if isfield(ePIE_inputs, 'miscNotes')
    miscNotes = ePIE_inputs.miscNotes;
else
    miscNotes = 'None';
end

if isfield(ePIE_inputs, 'dp_mask')
    dp_mask = ePIE_inputs(1).dp_mask;
    dp_mask = fftshift(dp_mask);
    dp_mask_flag = 1;
else
    dp_mask_flag=0;
end


if isfield(ePIE_inputs, 'showim')
    showim = ePIE_inputs(1).showim;
else
    showim = 0;
end

if isfield(ePIE_inputs, 'obj_scale'), obj_scale = ePIE_inputs.obj_scale;
else obj_scale = 1;
end

if isfield(ePIE_inputs, 'do_posi')
    do_posi = ePIE_inputs.do_posi;
else
    do_posi = 0;
end

if isfield(ePIE_inputs, 'save_path'), save_string = ePIE_inputs.save_path;
else save_string = [ pwd '/Results_ptychography/']; % Place to save results
end

if isfield(ePIE_inputs, 'save_intermediate')
    save_intermediate = ePIE_inputs.save_intermediate;
else save_intermediate = 0;
end

if isfield(ePIE_inputs, 'save_inter_freq')
    save_inter_freq = ePIE_inputs.save_inter_freq;
else save_inter_freq = 100;
end

clear ePIE_inputs

%% === Reconstruction parameters frequently changed === %%
nva = length(varargin);
optional_args(1:nva) = varargin;
[beta_obj, beta_ap] = optional_args{:};
freeze_aperture = Inf;
%% print parameters
fprintf('iterations = %d\n', iterations);
fprintf('beta probe = %0.1f\n', beta_ap);
fprintf('beta obj = %0.1f\n', beta_obj);
fprintf('gpu flag = %d\n', gpu);
fprintf('updating probe = %d\n', update_aperture);
fprintf('positivity = %d\n', do_posi);
fprintf('misc notes: %s\n', miscNotes);
%% Define parameters from data and for reconstruction
[N1,N2,nApert] = size(diffpats); % Size of diffraction patterns
norm_dp_avg = ( sum(diffpats(:).^2) ) / nApert;
for ii = 1:size(diffpats,3)
    diffpats(:,:,ii) = fftshift((diffpats(:,:,ii)));
end
diffpats = single(diffpats);


%% Get centre positions for cropping (should be a 2 by n vector)
%positions = convert_to_pixel_positions(positions,pixel_size);
[pixelPositions, bigx, bigy] = ...
    convert_to_pixel_positions_testing5(positions,pixel_size,N1);
centrey = round(pixelPositions(:,2));
centrex = round(pixelPositions(:,1));
centBig = round((bigx+1)/2);
Y1 = centrey - floor(N1/2); Y2 = Y1+N1-1;
X1 = centrex - floor(N1/2); X2 = X1+N1-1;
%% create initial aperture?and object guesses
if aperture == 0
    aperture = single(((2*makeCircleMask(round(aperture_radius./pixel_size),N1))));
    initial_aperture = aperture;
else
    aperture = single(aperture);
    initial_aperture = aperture;
end
aperture = complex( single( aperture/norm(aperture,'fro') * sqrt(norm_dp_avg) / N1 ));


if big_obj == 0
    %big_obj = single(rand(bigx,bigy)).*exp(1i*(rand(bigx,bigy)));
    big_obj = obj_scale*ones(bigx,bigy,'single') ;
    initial_obj = big_obj;
else
    big_obj = single(big_obj);
    initial_obj = big_obj;
end
fourier_error = zeros(iterations,nApert);

if gpu==1
    disp('========rPIE reconstructing with GPU========')
    diffpats = gpuArray(diffpats);
    fourier_error = gpuArray(fourier_error);
    big_obj = gpuArray(big_obj);
    aperture = gpuArray(aperture);
else
    disp('========rPIE reconstructing with CPU========')
end

sum_dp = zeros(nApert,1);
for aper=1:nApert
    current_dp = diffpats(:,:,aper);
    missing_data = current_dp == -1;
    sum_dp(aper) = sum(sum(current_dp(~missing_data)));
end

best_err = 100; % check to make sure saving reconstruction with best error

x_center = ceil(150); y_center = ceil(150);
[KK,JJ] = meshgrid(1:300,1:300);
R2 = ( (KK - x_center).^2 + (JJ - y_center).^2 );
Kfilter = exp( -R2 ./ (0.01* 300^2) );
Kfilter = ifftshift(Kfilter);

%% Main ePIE itteration loop
disp('========beginning reconstruction=======');
for itt = 1:iterations
    if itt>2.*iterations/3
        beta_obj = 0.95;
%         beta_ap = 0.95;
    end
    for aper = randperm(nApert)
        %         bigObjShifted = circshift(big_obj, [-1*(centrey(aper) - centBig) -1*(centrex(aper) - centBig)]);
        %         rspace = croppedOut(bigObjShifted,y_kspace);
        O_old = big_obj(Y1(aper):Y2(aper), X1(aper):X2(aper));
        O_max = max(abs(O_old(:))).^2;
        P_max = max(abs(aperture(:))).^2;

        %% Create new exitwave
        Pu_old = O_old.*(aperture);
        temp_dp = fft2(Pu_old);
        check_dp = abs(temp_dp);
        current_dp = diffpats(:,:,aper);

        if dp_mask_flag && itt>30,  missing_data = dp_mask;
        else missing_data = current_dp == -1;
        end

        k_fill = temp_dp(missing_data);
        temp_dp = current_dp.*exp(1i*angle(temp_dp)); % Replace reconstructed magnitudes with measured magnitudes
        temp_dp(missing_data) = k_fill;

        %% Update the object

        Pu_new = ifft2(temp_dp);
        diff_PO = Pu_new-Pu_old;
        %update_factor_ob = beta_obj/probe_max;
        %new_rspace = buffer_rspace + update_factor_ob.*conj(aperture).*(diff_exit_wave);
        O_new = O_old + conj(aperture).*diff_PO ./ ...
            ((1-beta_obj).*abs(aperture).^2 + beta_obj.*P_max);
        if do_posi == 1
            O_new = max(0, real(O_new));
        end
        big_obj(Y1(aper):Y2(aper), X1(aper):X2(aper)) = O_new;
        %% Update the probe

        if update_aperture
%             if itt > iterations - freeze_aperture
%                 new_beta_ap = beta_ap*sqrt((iterations-itt)/iterations);
%                 ds_P = new_beta_ap./O_max;
%             else
                ds_P = beta_ap./O_max;
%             end
            aperture = aperture + ds_P*conj(O_old).*(diff_PO);
            %aperture = real(aperture);
%             ap_lacey = aperture(lacey_supp);
%             ap_lacey = aperture(lacey_supp);
%             index_hi = abs(ap_lacey)>0.01*max(abs(ap_lacey(:)));
%             ap_lacey_mean = mean(ap_lacey(index_hi));
%             ap_lacey(index_hi) = ap_lacey_mean;
% %             aperture(lacey_supp) = ap_lacey;
%             aperture = ap_lacey;
        end

        switch probe_mask_flag
            case 0
            case 1 %real space
                aperture = aperture .* R_probe_mask;
            case 2 %k-space
                aperture = fftshift(ifft2(fft2(ifftshift(aperture)) .* F_probe_mask));
            otherwise
                error('probe mask flag can only be 0,1 or 2');
        end

        if sum_dp(aper)~=0
            fourier_error(itt,aper) = sum(abs( current_dp(~missing_data) - check_dp(~missing_data) )) ./ sum_dp(aper);
        end
    end
    %scale = max(max(abs(aperture))); aperture = aperture/scale;
    errors(itt) = sum(fourier_error(itt,:),2)/nApert;

    x_c = 678;
    y_c = 306;
    if itt<0
        rec_crop = big_obj(x_c+(-150:149),y_c+(-150:149));
        rec_crop = ifft2( fft2(rec_crop) .* Kfilter );
        big_obj(x_c+(-150:149),y_c+(-150:149)) = rec_crop;
    end

    if  mod(itt,2) == 0 && showim
        fprintf('%d. Error = %f\n',itt, errors(itt));

        %rec_crop = crop_roi(big_obj, 230, 296,668);
        rec_crop = big_obj(x_c+(-150:149),y_c+(-150:149));
        rec_crop = rec_crop(end:-1:1, end:-1:1);

        figure(33); img((abs(big_obj)),num2str(itt),aperture,'aperture', 'colormap', 'jet');
        %hsv_ap = make_hsv(aperture,1);
        %figure(34); imagesc(hsv_ap); title(num2str(itt)); colormap gray; axis image;
        %figure(35); img( abs(rec_crop),'rec mag',  angle(rec_crop), 'rec angle', 'abs','off')
        drawnow

    end

    %     toc
    %     tic

    mean_err = sum(fourier_error(itt,:),2)/nApert;
    if best_err > mean_err
        best_obj = big_obj;
        best_err = mean_err;
    end

    if save_intermediate == 1 && mod(itt,save_inter_freq) == 0
        save(sprintf('%s%s_itt%04d.mat',save_string,filename,itt),...%[save_string 'inter_output_rPIE_ID_',jobID,'_itt_',num2str(itt),'.mat'],...
            'best_obj','big_obj','aperture','fourier_error','-v7.3');
    end

end
disp('======reconstruction finished=======')


drawnow

if strcmp(gpu,'GPU')
    fourier_error = gather(fourier_error);
    best_obj = gather(best_obj);
    aperture = gather(aperture);
end

% if saveOutput == 1
%     save([save_string 'best_obj_' filename '.mat'],'best_obj','aperture','initial_obj','initial_aperture','fourier_error');
% end

%% Function for converting positions from experimental geometry to pixel geometry

%     function [positions] = convert_to_pixel_positions(positions,pixel_size)
%         positions = positions./pixel_size;
%         positions(:,1) = (positions(:,1)-min(positions(:,1)));
%         positions(:,2) = (positions(:,2)-min(positions(:,2)));
%         positions(:,1) = (positions(:,1)-round(max(positions(:,1))/2));
%         positions(:,2) = (positions(:,2)-round(max(positions(:,2))/2));
%         positions = round(positions);
%         bigx =little_area + max(positions(:))*2+10; % Field of view for full object
%         bigy = little_area + max(positions(:))*2+10;
%         big_cent = floor(bigx/2)+1;
%         positions = positions+big_cent;
%
%
%     end
    function [shiftDistancesX, shiftDistancesY, truePositions, positions, bigx, bigy] = ...
            convert_to_pixel_positions_testing3(positions,pixel_size,little_area)


        positions = positions./pixel_size;
        positions(:,1) = (positions(:,1)-min(positions(:,1)));
        positions(:,2) = (positions(:,2)-min(positions(:,2)));
        truePositions(:,1) = (positions(:,1) - max(positions(:,1))/2);
        truePositions(:,2) = (positions(:,2) - max(positions(:,2))/2);
        positions(:,1) = (positions(:,1)-round(max(positions(:,1))/2));
        positions(:,2) = (positions(:,2)-round(max(positions(:,2))/2));

        positions = round(positions);
        bigx =little_area + max(positions(:))*2+10; % Field of view for full object
        bigy = little_area + max(positions(:))*2+10;
        %         bigx = little_area;
        %         bigy = little_area;
        big_cent = floor(bigx/2)+1;
        positions = positions+big_cent;
        truePositions = truePositions + big_cent;

        shiftDistancesX = truePositions(:,1) - positions(:,1);
        shiftDistancesY = truePositions(:,2) - positions(:,2);



    end

    function [shiftDistancesX, shiftDistancesY, truePositions, positions, bigx, bigy] = ...
            convert_to_pixel_positions_testing4(positions,pixel_size,little_area)


        positions = positions./pixel_size;
        positions(:,1) = (positions(:,1)-min(positions(:,1)));
        positions(:,2) = (positions(:,2)-min(positions(:,2)));
        truePositions(:,1) = (positions(:,1) - round(max(positions(:,1))/2));
        truePositions(:,2) = (positions(:,2) - round(max(positions(:,2))/2));

        positions = round(truePositions);
        bigx =little_area + max(positions(:))*2+10; % Field of view for full object
        bigy = little_area + max(positions(:))*2+10;
        %         bigx = little_area;
        %         bigy = little_area;
        big_cent = floor(bigx/2)+1;
        positions = positions+big_cent;
        truePositions = truePositions + big_cent;

        shiftDistancesX = truePositions(:,1) - positions(:,1);
        shiftDistancesY = truePositions(:,2) - positions(:,2);



    end

    function [pixelPositions, bigx, bigy] = ...
            convert_to_pixel_positions_testing5(positions,pixel_size,little_area)


        pixelPositions = positions./pixel_size;
        pixelPositions(:,1) = (pixelPositions(:,1)-min(pixelPositions(:,1))); %x goes from 0 to max
        pixelPositions(:,2) = (pixelPositions(:,2)-min(pixelPositions(:,2))); %y goes from 0 to max
        pixelPositions(:,1) = (pixelPositions(:,1) - round(max(pixelPositions(:,1))/2)); %x is centrosymmetric around 0
        pixelPositions(:,2) = (pixelPositions(:,2) - round(max(pixelPositions(:,2))/2)); %y is centrosymmetric around 0

        bigx = little_area + round(max(pixelPositions(:)))*2+10; % Field of view for full object
        bigy = little_area + round(max(pixelPositions(:)))*2+10;

        big_cent = floor(bigx/2)+1;

        pixelPositions = pixelPositions + big_cent;


    end

%% 2D gaussian smoothing of an image

    function [smoothImg,cutoffRad]= smooth2d(img,resolutionCutoff)

        Rsize = size(img,1);
        Csize = size(img,2);
        Rcenter = round((Rsize+1)/2);
        Ccenter = round((Csize+1)/2);
        a=1:1:Rsize;
        b=1:1:Csize;
        [bb,aa]=meshgrid(b,a);
        sigma=(Rsize*resolutionCutoff)/(2*sqrt(2));
        kfilter=exp( -( ( ((sqrt((aa-Rcenter).^2+(bb-Ccenter).^2)).^2) ) ./ (2* sigma.^2) ));
        kfilter=kfilter/max(max(kfilter));
        kbinned = my_fft(img);

        kbinned = kbinned.*kfilter;
        smoothImg = my_ifft(kbinned);

        [Y, X] = ind2sub(size(img),find(kfilter<(exp(-1))));

        Y = Y-(size(img,2)/2);
        X = X-(size(img,2)/2);
        R = sqrt(Y.^2+X.^2);
        cutoffRad = ceil(min(abs(R)));
    end

%% Fresnel propogation
    function U = fresnel_advance (U0, dx, dy, z, lambda)
        % The function receives a field U0 at wavelength lambda
        % and returns the field U after distance z, using the Fresnel
        % approximation. dx, dy, are spatial resolution.

        k=2*pi/lambda;
        [ny, nx] = size(U0);

        Lx = dx * nx;
        Ly = dy * ny;

        dfx = 1./Lx;
        dfy = 1./Ly;

        u = ones(nx,1)*((1:nx)-nx/2)*dfx;
        v = ((1:ny)-ny/2)'*ones(1,ny)*dfy;

        O = my_fft(U0);

        H = exp(1i*k*z).*exp(-1i*pi*lambda*z*(u.^2+v.^2));

        U = my_ifft(O.*H);
    end

%% Make a circle of defined radius

    function out = makeCircleMask(radius,imgSize)


        nc = imgSize/2+1;
        n2 = nc-1;
        [xx, yy] = meshgrid(-n2:n2-1,-n2:n2-1);
        R = sqrt(xx.^2 + yy.^2);
        out = R<=radius;
    end

%% Function for creating HSV display objects for showing phase and magnitude
%  of a reconstruction simaultaneously

    function [hsv_obj] = make_hsv(initial_obj, factor)

        [sizey,sizex] = size(initial_obj);
        hue = angle(initial_obj);

        value = abs(initial_obj);
        hue = hue - min(hue(:));
        hue = (hue./max(hue(:)));
        value = (value./max(value(:))).*factor;
        hsv_obj(:,:,1) = hue;
        hsv_obj(:,:,3) = value;
        hsv_obj(:,:,2) = ones(sizey,sizex);
        hsv_obj = hsv2rgb(hsv_obj);
    end
%% Function for defining a specific region of an image

    function [roi] = get_roi(image, centrex,centrey,crop_size)

        bigy = size(image,1);
        bigx = size(image,2);

        half_crop_size = floor(crop_size/2);
        if mod(crop_size,2) == 0
            roi = {centrex - half_crop_size:centrex + (half_crop_size - 1);...
                centrey - half_crop_size:centrey + (half_crop_size - 1)};

        else
            roi = {centrex - half_crop_size:centrex + (half_crop_size);...
                centrey - half_crop_size:centrey + (half_crop_size)};

        end
    end

%% Fast Fourier transform function
    function kout = my_fft(img)
        kout = fftshift(fftn((ifftshift(img))));
    end
%% Inverse Fast Fourier transform function
    function realout = my_ifft(k)
        realout =fftshift((ifftn(ifftshift(k))));
        %realout =ifftshift((ifftn(fftshift(k))));

    end
end



