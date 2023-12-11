%{
======== in situ CDI simulation ======== 

1. define support shape and static structure
2. define model object using scattering factors of different elements and
other factors for more realistic simulation 
3. make noisy diffraction patterns
4. reconstruct using dose_reduction_reconstruct.m

%}

clear all;close all

%%
addpath model_and_support_functions_iCDI\
load('position.mat')
load('model_ptycho.mat')
%%
%%% parameters of ptycho
scannumx = 10;% ptycho scan
sump1=0;sump2=0;sum1to2=0;

%%%  flux 0.8273 = flux inside object FOV / total flux
reference_photons_series = (1/0.8273).*1.2e12/(scannumx.^2) * [0 0 0 0 0 1 1 1 1 1 1];
object_photons_series = (1/0.8273).*3*1e6/(scannumx.^2) * [1e0 1e1 1e2 1e3 1e4 1e0 1e1 1e2 1e3 1e4]; %100

%%% sample: thickness and size of FOV
thickness_ratio = 0.3;% thickness 1 means 1um; 0.3 means 300nm
p.h2o.thickness_um = 1*thickness_ratio; % will 0.2;um. will be converted into px
p.cell.thickness_um = 1*thickness_ratio; % will x0.2 2um. note that this is the MAX thickness, mean thickness is much smaller (see p.e.mean_thickness_um)
p.lacey.thickness_um = .02; %um

%%% parameters of support
sz = 1024; % CCD array size (px)
shift_xycell = 184; % 164
shift_xylacey = 256; 

% predefined functions
sp = @(x,y,z) subtightplot(x,y,z,[.07 .07],[.15 .15],[.15 .15]);
rrfac = @(in,model) sum(sum(abs(in-model)))/sum(sum(abs(model)));
krfac = @(in,model,beamstop) sum(sum(abs(in(~beamstop)-model(~beamstop))))/sum(sum(abs(model(~beamstop))));
interpolate = @(tabulated_energies,values,desired_energies) interp1(tabulated_energies,values,desired_energies,'linear');
        
% physical dimensions
avoNum = 6.022e23; % Avogadro's number

% SACLA parameters:
% Song et al., Multiple application X-ray imaging chamber for single-shot diffraction experiments with femtosecond X-ray laser pulses (2013)
p.e.QE                  = 1;%0.8; % detector efficiency at 710 eV < .2

%% experimental parameters
p.e.zz                  = 8; % fresnel from pinhole to object : um
step                    = 26;% 18 (256-l_pixel)/scannumx; % pixel scan number = 5*5
l_pixel                 = 106;%96
ovlap = sum(sum(single((makeCircleMask(l_pixel/2,256)+circshift(makeCircleMask(l_pixel/2,256),[step 0])>1))))./sum(sum(single(makeCircleMask(l_pixel/2,256))))
p.e.detector_size       = sz;
p.e.detector_px_size    = 10*1024/sz; % um
p.e.energy              = 5.3e2; % ev
p.e.lambda              = 1.2398/p.e.energy; % um
p.e.re                  = 2.82e-15*1e6; % classical electron radius (um)
p.e.CD                  = .05e6; % um

p.e.px_size             = p.e.lambda*p.e.CD/(p.e.detector_size*p.e.detector_px_size); % 0.0114um/px
p.e.l                   = p.e.px_size.*l_pixel; % /2 for ptycho
p.e.scanarea            = p.e.px_size.*256; % /2 for ptycho
sx                      =  (256-l_pixel-0)/2+40;
sy                      = (256-l_pixel-0)/2+40;
p.e.os                  = p.e.lambda*p.e.CD/(p.e.detector_px_size*p.e.l);
p.e.pulse_duration      = 1; % sec

p.e.resolution          = p.e.px_size*2; % real space resolution (um)
p.e.k_resolution        = 1/(2*p.e.px_size); % reciprocal space resolution (1/um)
p.e.diameter            = p.e.l; % um

%% make support
supp_radius = round( p.e.diameter / p.e.px_size / 2 );
supp_radius_scan = round( p.e.scanarea/ p.e.px_size / 2 );

m.lacey_supp = circshift(makeCircleMask(1*supp_radius,sz),[-shift_xylacey -shift_xylacey]);%...
m.lacey_supp(m.lacey_supp<max(m.lacey_supp(:))*.1)=0;

m.cell_supp = circshift(makeCircleMask(1*supp_radius,sz),[+shift_xycell +shift_xycell]);%...
m.cell_supp(m.cell_supp<max(m.cell_supp(:))*.1)=0;

m.lacey_suppB = circshift(makeSquareMask(1*supp_radius_scan,sz),[-shift_xylacey -shift_xylacey]);%...
m.cell_suppB = circshift(makeSquareMask(1*supp_radius_scan,sz),[+shift_xycell +shift_xycell]);%...

m.cell_suppBb = circshift(makeSquareMask(1*supp_radius_scan+8,sz),[+shift_xycell +shift_xycell]);%...

p.e.area  = nnz(m.cell_supp).*p.e.px_size^2;
p.e.scanarea2  = nnz(m.cell_suppB).*p.e.px_size^2;

p.lacey.thickness_px    = p.lacey.thickness_um / p.e.px_size; % px
p.h2o.thickness_px      = p.h2o.thickness_um / p.e.px_size; % px
p.cell.thickness_px     = p.cell.thickness_um / p.e.px_size; % px

%% projection model in pixel
h2o_2d_projection  = thickness_ratio*model_ptycho.h2o_2d_projection ;
h2o_2d_projection = circshift( h2o_2d_projection, [shift_xycell,shift_xycell]);

cell_2d_projection = thickness_ratio*model_ptycho.cell_2d_projection ;
cell_2d_projection = circshift( cell_2d_projection, [shift_xycell,shift_xycell]);

lacey_2d_projection = model_ptycho.lacey_2d_projection;
lacey_2d_projection = p.lacey.thickness_px.*lacey_2d_projection;
h2o_2d_projection_ss = ones(p.e.detector_size) * p.h2o.thickness_px - lacey_2d_projection;

%% scattering factors for common protein elements
scattering_factors;
element_type={'h' 'c' 'n' 'o' 'au' 'h2o' 's'};
for nn=1:length(element_type)
    Type=element_type{nn};
    eval(sprintf('f1.%s=interpolate(sc.%s.f1(:,1)*1e3,sc.%s.f1(:,2),p.e.energy);',Type,Type,Type));
    eval(sprintf('f2.%s=interpolate(sc.%s.f2(:,1)*1e3,sc.%s.f2(:,2),p.e.energy);',Type,Type,Type));
end

N_cell = 50+30+9+10+1;
M_cell = 50*1 + 30*12 + 9*14 + 10*16 + 1*32;
f1.all=(50*f1.h + 30*f1.c + 9*f1.n + 10*f1.o + 1*f1.s)/N_cell; % e/atom
f2.all=(50*f2.h + 30*f2.c + 9*f2.n + 10*f2.o + 1*f2.s)/N_cell; % e/atom
f1.h2o=(2*f1.h + 1*f1.o)/3; % e/atom
f2.h2o=(2*f2.h + 1*f2.o)/3; % e/atom

%% model: atoms/pixel ==> e/pixels
% The total number of projected electrons in each object pixel.
% This is calculated via dimensional analysis. atoms/um3 * um3/vox
% * vox/px3 * px * e/atom = e/px2
% n2 is an electron density map of the object in e/px2
% n is electron density map of object in e/um2

%%%% water
p.h2o.mol_weight = 2*1 + 1*16; % g/mol
p.h2o.density = 1; % 1g/cm^3
p.h2o.Natoms = 3;
p.h2o.na = p.h2o.Natoms * p.h2o.density * avoNum / p.h2o.mol_weight / 1e12; % atoms/um^3
%%%% water in cell
p.h2o.n    = p.h2o.na * p.e.px_size^3 * h2o_2d_projection * (f1.h2o + 1i*f2.h2o);%cell
m.h2o = p.h2o.n;
%%%% water in lacey
p.h2o.n_ss = p.h2o.na * p.e.px_size^3 * h2o_2d_projection_ss * (f1.h2o + 1i*f2.h2o);
m.h2o_ss = p.h2o.n_ss;

%%%% cell
p.cell.mol_weight = 50 + 30*12 + 9*14 + 10*16 + 32; % g/mol
p.cell.density = 1.35; % g/cm^3
p.cell.N = 50+30+9+10+1;
p.cell.na = p.cell.N * p.cell.density * avoNum / p.cell.mol_weight / 1e12; % atoms/um^3
p.cell.n = p.cell.na * p.e.px_size^3 * cell_2d_projection * (f1.all + 1i*f2.all);
m.cell = p.cell.n;
%%%% lacey
p.lacey.mol_weight = 197; %g/mol
p.lacey.density = 19.3; %g/cm^3
p.cell.Nlacey = 1;
p.lacey.na = p.cell.Nlacey * p.lacey.density * avoNum / p.lacey.mol_weight / 1e12;
p.lacey.n = p.lacey.na * p.e.px_size^3 * lacey_2d_projection * (f1.au + 1i*f2.au);
m.lacey = p.lacey.n;

fprintf('cell thickness: %.2f nm \nlacey thickness: %.2f nm \nH2O thickness: %.2f nm\n',p.cell.thickness_um*1e3,1e3*p.lacey.thickness_um,1e3*p.h2o.thickness_um)

%% final model and support
ss = (m.h2o_ss+m.lacey).*m.lacey_suppB;
dd = (m.h2o + m.cell).*m.cell_suppB;
m.supp = m.cell_supp | m.lacey_supp;     

%% difference flux
for run_id = [7]  %  1:4 no lacey 6:9 with lacey
    dp_noisy = zeros(1024,1024,scannumx^2);
    if run_id<=5;laceyflag=0;ref=0;else;ref=1;laceyflag=1;end

    % flux
    ref_flux = reference_photons_series(run_id) 
    obj_flux = object_photons_series(run_id) 
    p.e.Io_cell             = obj_flux * p.e.pulse_duration; % ph. total deposited photons.
    p.e.Io_lacey            = ref_flux * p.e.pulse_duration; % ph. total deposited photons
    p.e.true_flux_ds        = p.e.Io_cell / p.e.area; % ph/um2 on the cell
    p.e.true_flux_ss        = p.e.Io_lacey / p.e.area;

    % flux_scaling--scale intensities going through each pinhole
    flux_scaling = ((p.e.re*p.e.lambda) / (p.e.os*p.e.diameter)); % scaling
    p.cell.scaling = sqrt(p.e.true_flux_ds * p.e.QE) * flux_scaling; % scale flux to cell
    p.lacey.scaling = sqrt(p.e.true_flux_ss * p.e.QE) * flux_scaling; % scale flux to static structure
   
    scanid=0;
    
    for scanx=0:scannumx-1
    for scany=0:scannumx-1
        scanid = scanid+1; 
%         rand1=floor(10-20*rand(1));rand2=floor(10-20*rand(1));
%         centerx(scanid) =  sz/2-(step*scanx-sx) +rand1;%
%         centery(scanid) = sz/2-(step*scany-sy) +rand2;%
        centerx(scanid) =  position.cx(scanid);% random scan  
        centery(scanid) = position.cy(scanid);% random scan 
        shfx = sz/2-centerx(scanid);
        shfy = sz/2-centery(scanid);
            
        %%%%%%%%%%%%%%%%%%%%%%% model and probe  %%%%%%%%%%%%%%%%%%%%%%%
        if laceyflag
        modeltruth = ss+dd;
        ssdd = circshift(modeltruth,[shfx,shfy]);
        % probe            
        probe  = fresnel_advance( p.lacey.scaling.*m.lacey_supp+p.cell.scaling.* ...
        m.cell_supp,p.e.px_size,p.e.px_size,p.e.zz,p.e.lambda);         
        else
        modeltruth = dd;
        ssdd = circshift(modeltruth,[shfx,shfy]);
        % probe            
        probe  = fresnel_advance( p.cell.scaling.* m.cell_supp,p.e.px_size,p.e.px_size,p.e.zz,p.e.lambda);
        end            
        % define model            
        m.truth = probe.*(ssdd);

        %%%%%%%%%%%%%%%%%%%%%%% CDI--- exit wave & diffraction pattern  %%%%%%%%%%%%%%%%%%%%%%%
        PO = m.truth;
        kPO = my_fft( PO );
        d.dp = abs(kPO).^2;   
        %%% check flux
        p.e.integrated_counts = sum(sum(d.dp));
        
        dc = zeros(sz);
        dc(sz/2+1,sz/2+1)=true;
        d.dp(dc == 1) = 0;

        %%% add Poisson noise
        d.dp_noisy = poissrnd(d.dp);
        d.dp_noisy(d.dp_noisy<0)=0;
        d.dp_noisy=sqrt(d.dp_noisy);
        beamstop = logical(pad2(ones(5),[sz,sz]));

        d.dp_noisy(dc == 1)=-1; % ignore DC
        d.dp_noisy(d.dp_noisy == 0) = 0; % -1, should zero valued pixels be enforced?
        dp_noisy(:,:,scanid) = d.dp_noisy;
        d.dp(dc == 1 ) = -1;
        d.rfac=krfac(d.dp_noisy, sqrt(d.dp), beamstop); % calculate R-factor of noisy pattern
        p.e.dp_rfac = krfac(d.dp, d.dp_noisy.^2, beamstop);

        %%%%%%%%%%%%%%%%%%%%%%% Perfect Lens  %%%%%%%%%%%%%%%%%%%%%%%
        if scanid==1;
        %%% total flux p.e.Io_cell2 is flux used for perfect lens
        p.e.Io_cell2 = p.e.Io_cell.*0.8273; %sump2/sump1=0.  
        
        % absorption
        na_f2_z = p.e.px_size* ( p.h2o.na.*f2.h2o.*h2o_2d_projection + p.cell.na.*f2.all.*cell_2d_projection );
        miu_z = (4*pi./p.e.lambda) *(p.e.re*p.e.lambda.^2)./(2*pi) * na_f2_z; %um
        % absorp = 1; %phase only
        absorp = sqrt(exp(-miu_z)) ;% consider absorption
        p.absorp = absorp;
        
%         fprintf('absorp: min = %.2f, max = %.2f\n',min(min(absorp(absorp~=0))),max(max(absorp(absorp~=0))))
        
        % phase
        na_f_z = p.e.px_size* ( p.h2o.na.*f1.h2o.*h2o_2d_projection + p.cell.na.*f1.all.*cell_2d_projection );
        delta_z = (p.e.re*p.e.lambda.^2)./(2*pi) * na_f_z; %um
        phi = 2*pi./p.e.lambda *delta_z ;
        phasecontrast = max(phi(:))-min( phi(phi~=0));
        p.phi = phi;

%         fprintf('phi: min=%.2f, max= %.2f, contras = %.2f.\n',min( phi(phi~=0)), max(phi(:)),phasecontrast);
        
        % PL model
        perfect_lens_image = (absorp.^2+2*absorp.*(phi)).*m.cell_suppB;
        
        % flux scale
        % apply correct scaling factor so total photon counts on the
        % detector is consistent with total flux right after the sample
        sumimage = sum(perfect_lens_image(:));
        absorp_inten_ratio = sum(sum(absorp.^2.*m.cell_suppB))./sum(sum(m.cell_suppB));
        %%%% total flux = p.e.Io_cell2*scannumx^2;
        perfect_lens_image = ( absorp_inten_ratio*p.e.QE * p.e.Io_cell2*scannumx^2/sumimage).*perfect_lens_image;

        perfect_lens_image_noise = poissrnd(perfect_lens_image);
        per_image_noise = perfect_lens_image_noise;%crop_roi(perfect_lens_image_noise,[512 512]);
        
        end
        
        %%% add Poisson noise
        im_noise = per_image_noise;
        im_noise1 = im_noise(im_noise>0);
        max_im1 = max(im_noise1(:));
        min_im1 = min(im_noise1(:));

        %%% dose
        %%% calculate dose imparted on cell
        eval(sprintf('protein_AbCoef=interpolate(sc.protein.mu(:,1)*1e3,sc.protein.mu(:,2),p.e.energy);'));
        % dose = flux * mu*E/rho, where mu is linear absorption
        % coefficient, E is photon energy, and rho is cell density. mu/rho
        % is the mass absorption coefficient. Last 3 terms are unit
        % conversion from g to kg, um to cm, and from eV to J.%             
        % miu_cell =  (4*pi./p.e.lambda) *(p.e.re*p.e.lambda.^2)./(2*pi) * ( p.cell.na.*f2.all ); %um
        % p.e.dose = p.e.true_flux_ds * p.e.energy * miu_cell / p.cell.density * 1e3 * (1e4)^3 * 1.6e-19; % not complete. see Rodriguez 2015
        p.e.dose = p.e.true_flux_ds * p.e.energy * protein_AbCoef / p.cell.density * 1e3 * (1e4)^2 * 1.6e-19; % not complete. see Rodriguez 2015

        %%% transmission
        % cross_section = 2 * p.e.re * p.e.lambda * f2.all; % source: http://henke.lbl.gov/optical_constants/intro.html
        % map = p.cell.thickness_um*cell_2d_projection;
        % p.e.mean_thickness_um = mean(mean(map(map~=0)));
        % p.e.transmission = exp(-p.cell.na * cross_section * p.e.mean_thickness_um);

        figure(2); img(m.supp,'m.supp',probe,'',log(1+1e3.*abs(ssdd.*probe)),'',(modeltruth),'m.truth',per_image_noise,'');%,log(1+d.dp_noisy),'');
        caxis([min_im1,max_im1]);colormap(jet);drawnow
        
        if scanid==1;
        fprintf('dose = %.2s, total flux on cell = %.2s, flux per um2 = %.2s\n ',...
            p.e.dose, p.e.Io_cell2*scannumx.^2,p.e.Io_cell2*scannumx.^2./(p.e.scanarea.^2));
        end
        p.e.Io_cell_um2 = p.e.Io_cell2*scannumx.^2./(p.e.scanarea.^2);
        % calculate flux in FOV and flux leak         
        probecell=fresnel_advance( m.cell_supp,p.e.px_size,p.e.px_size,p.e.zz,p.e.lambda);
        probelacey=fresnel_advance( m.lacey_supp,p.e.px_size,p.e.px_size,p.e.zz,p.e.lambda);
        m.cell_suppB2 = circshift(m.cell_suppB,[shfx,shfy]);
        probecell2=probecell.*m.cell_suppB2;
        probelacey2=probelacey.*m.cell_suppB2;
        sump1 = sump1+sum(abs(probecell(:)).^2);
        sump2 = sump2+sum(abs(probecell2(:)).^2);
        sum1to2 = sum1to2+sum(abs(probelacey2(:)).^2);

        end
    end
%     fprintf('percentage of flux in FOV = %.2f %%, percentage of flux leak to cell = %.2s %%  \n',100*sump2/sump1,100*sum1to2/sump1);
    save (['ptycho10x10_300dp_ref',num2str(ref),'_3.5e',num2str(5-5*ref+run_id-1),'.mat'],'dp_noisy','per_image_noise','p','centerx','centery','modeltruth','-v7.3')
end



