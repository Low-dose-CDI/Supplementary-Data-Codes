function [O2D_best, probe] = ePIE_insitu2(inputs)
% Performs phase retrieval on a series of in situ CDI diffraction patterns
%   Inputs:
%       diff_pats - stack of in situ CDI diffraction patterns
%       probe - pre-determined probe function
%       ref_mask - mask for cropping out static structure
%       iterations - # of iterations to run
%       showprogress - flag for displaying reconstruction progress
%   Outputs:
%       R - reconstructed series of objects
%       r_factor - Fourier R-factor
%       params - data-specific parameters. Used for display purposes

diffpatts  = ifftshift(ifftshift(inputs.diffpats,1),2);
probe      = inputs.probe;
support_ref= inputs.support_ref;
support_obj= inputs.support_obj;
support_all = support_ref|support_obj;
non_support_all = ~support_all;
non_support_ref = ~support_ref;

iterations = inputs.iterations;

unknown = diffpatts(:,:,1)==-1;
known = ~unknown;


[N1, N2, num_frames] = size(diffpatts);
recs = rand(N2, N1, num_frames,'single');
Ys   = zeros(N2,N1,num_frames,'single');
ref_area = single(support_ref);

if isfield(inputs,'initial')
    recs = repmat(inputs.initial,[1,1,num_frames]);
end

if isfield(inputs,'sigma')
    sigma = inputs.sigma;
else
    sigma=0.1;
end
if isfield(inputs,'thresh')
    thresh = inputs.thresh; 
else  
    thresh=0.3; 
end

%norm_dp_avg = sum(diffpatts(:)) / num_frames;
%probe = probe/norm(probe,'fro') * sqrt(norm_dp_avg) / N1;

%d_probe = 0.5*abs(probe) / max(abs(probe(:))) .* conj(probe) ./ (eps + abs(probe).^2);
%d_probe =  conj(probe) ./ max (abs(probe(:)).^2);

sum_dp = zeros(num_frames,1);
for i=1:num_frames
    dp_i      = diffpatts(:,:,i);
    sum_dp(i) = sum(dp_i(known));
end

x_center = ceil(N1/2); 
y_center = ceil(N2/2);
X = 1:iterations; filtercount = 10;
sub_iter = iterations/filtercount;
FX=(filtercount + 1-ceil(X*filtercount/iterations))*ceil(iterations/(1*filtercount));
FX=((FX-ceil(iterations/filtercount))*(2*N1)/max(FX))+(2*N1/10);
[KK,JJ] = meshgrid(1:N1,1:N2);
R2 = (KK - x_center).^2 + (JJ - y_center).^2 ;

if true
    probe = gpuArray(probe);
    support_ref = gpuArray(support_ref);
    support_obj = gpuArray(support_obj);
    recs = gpuArray(recs);
    ref_area = gpuArray(ref_area);
end
%d_probe =  conj(probe) ./ (0.9*abs(probe).^2 + 0.1*max (abs(probe(:)).^2) );

O2D_best = recs;
Y2D_best = zeros(N1,N2,num_frames);
errorF_best = ones(num_frames,1);
dt=1;
ds=1;
h = fspecial('disk', 5);
% thresh = 0.1;
%% loops
for t = 1:iterations
    errors=0;
    if t==round(0.1*iterations), sigma=0.2; end

    for i = randperm(num_frames)
        dp_i = diffpatts(:, :, i);
        O_i = recs(:,:,i);
        Y_i = Ys(:,:,i);

        % the novel part
        if t<3000
            O_i(support_ref) = 0.5.*ref_area(support_ref) + 0.5 * O_i(support_ref);%change ref to lacey
            ref_area = support_ref .* O_i;
        end

        O_new = O_i - dt*Y_i;
        % Fourier projection
        PO = O_new .* probe;

        Zi = fft2(PO);        
        Z_i_missing = Zi(unknown);
        %abs_Z = abs(Zi);
        angles = angle(Zi);
        Zi = (dp_i .* exp(1j*angles) + sigma*Zi) ./ (1+sigma);
        Zi(unknown) = Z_i_missing;

        PO_new = ifft2(Zi);
        diff_PO = PO_new - PO;

        %d_probe =  conj(probe) ./ (0.9*abs(probe).^2 + 0.1*max (abs(probe(:)).^2) );
        O_new = O_new + diff_PO .* conj(probe) ./ (0.9*abs(probe).^2 + 0.1*max (abs(probe(:)).^2) );
        O_new (support_obj) =  (O_new (support_obj));
        %%
% if t<=round(iterations)
%         O_new2 = abs(O_new);
%         O_new2 = conv2(O_new2, h, 'same');
%         O_new_ref = O_new2.*support_ref;
% 
%         O_new_ref = O_new_ref/max(O_new_ref(:)); % make in range [0,1]
%         O_new_ref = (O_new_ref > thresh) | non_support_ref;
%   
%         O_new = O_new.*O_new_ref - 0*O_new;
% end
%%
        %O_new=O_new.*support_all;

%         if t>30000
%             d_O = conj(O_i) .* diff_PO ./ max( abs(O_i(:)).^2);
%             probe = probe + (sqrt( (iterations-t)/iterations ) ) *d_O;
%         end

        Y_i = Y_i + ds*(1.5*O_new - 0.5*O_i);

        % update y
%         if isreal
%             index= support_obj & real(Y_i)>0;
%             Y_i(index) = 1j*imag(Y_i(index));
%         end
        Y_i(support_all) = 0;
        %y(support_real) = 0;

        % smooth out z FX(t)^2)
        if mod(t,sub_iter)==1
            kfilter = exp( -R2 /(2* FX(t)^2 )); kfilter = kfilter/max(max(kfilter));
            kfilter_shift = ifftshift(kfilter);
        end
        % type 2: Poisson noise removal on Fourier space
        %y_temp = y.*kfilter; y(index0) = y_temp(index0);
        %Y_i = Y_i.*kfilter;
        
%         Fy = fft2(Y_i) .* kfilter_shift; 
%         Y_i = ifft2(  Fy );

        rec_i = support_all .* O_new;
        z_temp = abs(  fft2(rec_i.*probe) );
        errorF = sum( abs( z_temp(known) - dp_i(known) )) /  sum_dp(i) ;
        if errorF < errorF_best(i)
            errorF_best(i) = errorF;
            O2D_best(:,:,i) = O_new;
            Y2D_best(:,:,i) = Y_i;
        end

        if mod(t,sub_iter)==0
            Y_i   = Y2D_best(:,:,i);
            O_new = O2D_best(:,:,i);
        end

        % store updated object and probe
        Ys(:,:,i) = Y_i;
        recs(:, :, i) = O_new;
        errors = errors + errorF ;

        if mod(t, 10) == 0 && i==1
            %R_f = errors/num_frames;
            fprintf('%d. error = %.4f\n',t,errorF);
            figure(10); img(O_new.*support_ref,'ref', real(O_new).*support_obj,'mag', imag( O_new.*support_obj),'imag','colormap','gray');
            drawnow
        end % display progress

    end % frame



end % iterations



end % function










