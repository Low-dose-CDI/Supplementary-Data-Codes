function [u2D_best,y, support_refined] = GPS_shrinkwrap2(input)
diffpat = input.diffpat;
[n1, n2]      = size(diffpat)   ;

support_obj   = input.support_obj;
support_ref   = input.support_ref;
ntimes        = input.iterations; 

if isfield(input,'sigma')
    sigma = input.sigma;
else
    sigma=0.1;
end

if isfield(input,'isreal')
    isreal=true;
else
    isreal=false;
end

support_all     = support_ref | support_obj;


non_support_all = ~support_all;
non_support_ref = ~support_ref;

if isfield(input,'thresh'),  thresh = input.thresh; else  thresh=0.3; end

x_center = ceil(n1/2); y_center = ceil(n2/2);

% filter parameter
X = 1:ntimes; filtercount = 10;
sub_iter = ntimes/filtercount;
FX=(filtercount + 1-ceil(X*filtercount/ntimes))*ceil(ntimes/(1*filtercount));
FX=((FX-ceil(ntimes/filtercount))*(2*n1)/max(FX))+(2*n1/10);
[KK,JJ] = meshgrid(1:n1,1:n2);
R2 = (KK - x_center).^2 + (JJ - y_center).^2 ;

%% OSS
Z_initial = diffpat;
stopper = find(diffpat==-1);
unknown = diffpat==-1; known = ~unknown;
Z_initial (stopper) = 0;


if isfield(input,'beta')
    beta = input.beta;
else
    beta = 1;
end

if isfield(input, 'initial')
    u = input.initial;
else
    u =  support_ref.*randn(size(support_ref));
    u = u /norm(u,'fro')* norm(diffpat,'fro');
end
u2D_best = u;
Y2D_best = zeros(size(u));

fprintf('Done with initialization, starting loop\n')

y = zeros(n1,n2);
dt=1; ds=1;
errorF_best=1;

h = fspecial('disk', 5);
h2 = fspecial('gaussian',7,1);
%thresh = 0.65;

%% Loop
for t = 1:ntimes
    if t==round(0.8*ntimes), sigma=0.1; end
    
    u_new = u - dt*y;        
    z = fftshift(fft2(u_new));
    z_stopper = z(stopper);
    phase = angle(z);    
    z = ( diffpat.* exp(1i*phase) + sigma*z ) / (1+sigma);    
    z(stopper) = z_stopper;

    u_new = ifft2(ifftshift(z));
    if mod(t,5)==-1
        u_new = smooth3D(u_new,0.4);
    end    
    
    if mod(t,1)==-1, u_new( abs(u_new) <0.09 ) = 0; end

    u_hat = 2*u_new - u;
    u = u_new;

    if t<=round(ntimes)
        u_abs = abs(u_hat);
        u_abs = conv2(u_abs, h, 'same');
        support_refined = u_abs.*support_ref;

        support_refined = support_refined/max(support_refined(:)); % make in range [0,1]
        support_refined = support_refined > thresh;
        support_refined1 = support_refined | non_support_ref;
		
        u_hat = u_hat.*support_refined1;% - 0*u_hat;
        if mod(t,50)==0
            figure(201); img(support_refined1,''); drawnow;
        end
    end
    
    y = y + ds*u_hat;

    %y_hat = real(y); y_hat(y_hat<0 | index_real0) = 0;    
    %y =  y - y_hat;
    if isreal
        index= support_ref & real(y)>0;
        y(index) = 1j*imag(y(index)); 
    end    
    y(support_all) = 0;
    %y(support_real) = 0;
    
    % smooth out z FX(t)^2)
    if mod(t,sub_iter)==1
        kfilter = exp( -R2 /(2* FX(t)^2 )); kfilter = kfilter/max(max(kfilter));
    end
    % type 2: Poisson noise removal on Fourier space
    %y_temp = y.*kfilter; y(index0) = y_temp(index0);
    y = y.*kfilter;
    
    % type 1: low pass filter
    %Fy = my_fft(y) .* kfilter; y = my_ifft(Fy);
    %Fy = fftshift(fft2((y))) .* kfilter; y = (ifft2(ifftshift(Fy)));
    
    % type 3: apply type1 and type2 alternatively
    %if mod(t,2)==0, Fy = fftshift(fft2(y)) .* kfilter; y = ifft2(ifftshift(Fy)); else y = y.*kfilter; end
    
    
    if mod(t,sub_iter)==0 && t>10*ntimes
        u = u2D_best;
        y = Y2D_best;
    end
    
    if mod(t,2)==0
        %u = ifft2(ifftshift(z));
        rec = support_all .* u;
        z_temp = abs( my_fft(rec));
        errorF = sum(sum(abs(z_temp(known)-diffpat(known)))) / sum(sum(diffpat(known)));
        
        if errorF <= errorF_best
            errorF_best = errorF;
            u2D_best = u;
            Y2D_best = y;
        end
                
    end
    
    if mod(t,100)==0
        fprintf('%d. ErrorF = %f\n', t, errorF) ;
        figure(101); img(rec.*support_ref,'ref',rec.*support_obj, 'mag', imag(rec.*support_obj),'imag', 'colormap','jet'); drawnow();
    end
    
end

fprintf('All done!\n')

%% Show some results
%diff = known.*abs(diffpat - abs(fftshift(fft2(real(sample)))) );
%z_in = sample(r1:r2, c1:c2);
%figure(1); imagesc(z); axis image;
%figure(2); imagesc(sample); axis image;
%figure(3); imagesc(diff); axis image; colorbar('fontsize',20); set(gca,'xtick',[],'ytick',[]); colormap(jet) %title('Residual','fontsize',18);
%figure(4); imagesc(z_in); axis image; set(gca,'xtick',[],'ytick',[]); colormap(jet)

%%

%{
figure; p3 = bar(-log10(x),y3); ylim([0 .2]);
p3(1).FaceColor = [0, 0.4470,0.7410]; p3(2).FaceColor = [0.8500, 0.3250, 0.0980];
set(gca,'Xtick',-9:-7); %// adjust manually; values in log scale
set(gca,'Xticklabel',10.^(-get(gca,'Xtick'))); %// use labels with linear values
set(gca, 'fontsize',18);
%}


