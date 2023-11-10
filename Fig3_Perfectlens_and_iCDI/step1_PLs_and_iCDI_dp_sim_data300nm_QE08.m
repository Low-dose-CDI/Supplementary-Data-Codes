%{
======== step 1: PLs and in situ CDI dps simulation ======== 

1. define support shape and static structure
2. define model object using scattering factors of different elements and
other factors for more realistic simulation 
3. make noisy diffraction patterns
4. make perfectlens images using I=a^2±2a*phi and add noise

%}

clear all
close all
addpath model_and_support_functions_iCDI\
load('model_iCDI.mat')
%%
%change thickness
thickness_ratio = 0.3;%um
p.h2o.thickness_um = 1*thickness_ratio; %  um
p.cell.thickness_um = 1*thickness_ratio; %  um
p.lacey.thickness_um = .02; %um

%flux
reference_photons_series = 1.2e12 * [0 0 0 0 0 1 1 1 1 1 1 1 1 1];
object_photons_series = 3*1e6 * [1e0 1e1 1e2 1e3 1e4 1e0 1e1 1e2 1e3 1e4 1e5 1e6]; %1.108*

%other parameters
sz = 1024; % CCD array size (px)
shift_xycell = 124; shift_xylacey = 236;sz2 = 270; % size of inputmodel
avephi0 = 0;avephimax=0;%check phi of cell

% physical dimensions
avoNum = 6.022e23; % Avogadro's number

% experimental parameters
p.e.detector_size       = sz;
p.e.detector_px_size    = 10*1024/sz; % um
p.e.energy              = 5.3e2; % ev
p.e.lambda              = 1.2398/p.e.energy; % um
p.e.re                  = 2.82e-15*1e6; % classical electron radius (um)
p.e.CD                  = .05e6; % um

p.e.px_size             = p.e.lambda*p.e.CD/(p.e.detector_size*p.e.detector_px_size); % 0.0114um/px
p.e.l                   = p.e.px_size.*256;%(sz/desired_os * p.e.px_size / 2); % um. field of view. /2 is to account for 2 pinholes
p.e.os                  = p.e.lambda*p.e.CD/(p.e.detector_px_size*p.e.l);
p.e.pulse_duration      = 1; % sec
p.e.resolution          = p.e.px_size*2; % real space resolution (um)
p.e.k_resolution        = 1/(2*p.e.px_size); % reciprocal space resolution (1/um)
p.e.diameter            = p.e.l; % um
p.lacey.thickness_px    = p.lacey.thickness_um / p.e.px_size; % px
p.h2o.thickness_px      = p.h2o.thickness_um / p.e.px_size; % px
p.cell.thickness_px     = p.cell.thickness_um / p.e.px_size; % px

%%% make support
supp_radius = round( p.e.diameter / p.e.px_size / 2 );
m.lacey_supp = circshift(makeSquareMask(1*supp_radius,sz),...
            [-shift_xylacey -shift_xylacey]);%...
m.cell_supp = circshift(makeSquareMask(1*supp_radius,sz),...
            [+shift_xycell +shift_xycell]);%...        
m.cell_supp(m.cell_supp<max(m.cell_supp(:))*.1)=0;
m.lacey_supp(m.lacey_supp<max(m.lacey_supp(:))*.1)=0;
p.e.area  = nnz(m.cell_supp).*p.e.px_size^2;

% SACLA parameters:
% Song et al., Multiple application X-ray imaging chamber for single-shot diffraction experiments with femtosecond X-ray laser pulses (2013)
p.e.QE                  = 0.8;%0.8; % detector efficiency at 710 eV < .2


% predefined functions
sp = @(x,y,z) subtightplot(x,y,z,[.07 .07],[.15 .15],[.15 .15]);
rrfac = @(in,model) sum(sum(abs(in-model)))/sum(sum(abs(model)));
krfac = @(in,model,beamstop) sum(sum(abs(in(~beamstop)-model(~beamstop))))/sum(sum(abs(model(~beamstop))));
interpolate = @(tabulated_energies,values,desired_energies) interp1(tabulated_energies,values,desired_energies,'linear');

%%% scattering factors for common protein elements
scattering_factors;
element_type={'h' 'c' 'n' 'o' 'au' 'h2o' 's'};
for nn=1:length(element_type)
    Type=element_type{nn};
    eval(sprintf('f1.%s=interpolate(sc.%s.f1(:,1)*1e3,sc.%s.f1(:,2),p.e.energy);',Type,Type,Type));
    eval(sprintf('f2.%s=interpolate(sc.%s.f2(:,1)*1e3,sc.%s.f2(:,2),p.e.energy);',Type,Type,Type));
end

% The total number of projected electrons in each object pixel.
% This is calculated via dimensional analysis. atoms/um3 * um3/vox
% * vox/px3 * px * e/atom = e/px2
N_cell = 50+30+9+10+1;
M_cell = 50*1 + 30*12 + 9*14 + 10*16 + 1*32;
f1.all=(50*f1.h + 30*f1.c + 9*f1.n + 10*f1.o + 1*f1.s)/N_cell; % e/atom
f2.all=(50*f2.h + 30*f2.c + 9*f2.n + 10*f2.o + 1*f2.s)/N_cell; % e/atom
f1.h2o=(2*f1.h + 1*f1.o)/3; % e/atom
f2.h2o=(2*f2.h + 1*f2.o)/3; % e/atom

p.h2o.mol_weight = 2*1 + 1*16; % g/mol
p.h2o.density = 1; % 1g/cm^3
p.h2o.Natoms = 3;
p.h2o.na = p.h2o.Natoms * p.h2o.density * avoNum / p.h2o.mol_weight / 1e12; % atoms/um^3

p.cell.mol_weight = 50 + 30*12 + 9*14 + 10*16 + 32; % g/mol
p.cell.density = 1.35; % g/cm^3
p.cell.N = 50+30+9+10+1;
p.cell.na = p.cell.N * p.cell.density * avoNum / p.cell.mol_weight / 1e12; % atoms/um^3

p.lacey.mol_weight = 197; %g/mol
p.lacey.density = 19.3; %g/cm^3
p.cell.Nlacey = 1;
p.lacey.na = p.cell.Nlacey * p.lacey.density * avoNum / p.lacey.mol_weight / 1e12;

%%
for run_id = 6:9
    for time_id = 1:20

        % tunables
        if reference_photons_series(run_id) == 0
            gold_lacey = 0; % use lacey gold for reference struct?
        else
            gold_lacey = 1;
        end

        %%% 
        ref_flux = reference_photons_series(run_id);
        obj_flux = object_photons_series(run_id);       
        p.e.Io_cell             = obj_flux * p.e.pulse_duration; % ph. total deposited photons.
        p.e.Io_lacey            = ref_flux * p.e.pulse_duration; % ph. total deposited photons
        p.e.true_flux_ds        = p.e.Io_cell / p.e.area; % ph/um2 on the cell
        p.e.true_flux_ss        = p.e.Io_lacey / p.e.area;

        %%% dynamic model
        h2o_2d_projection  = thickness_ratio*model_iCDI.h2o_2d_projection{time_id};
        cell_2d_projection = thickness_ratio*model_iCDI.cell_2d_projection{time_id};

        cell_2d_projection = circshift(cell_2d_projection, [shift_xycell,shift_xycell]);
        h2o_2d_projection = circshift(h2o_2d_projection, [shift_xycell,shift_xycell]);

        if gold_lacey
            lacey_2d_projection = p.lacey.thickness_px.*model_iCDI.lacey_2d_projection{time_id};
            h2o_2d_projection_ss = ones(p.e.detector_size) * p.h2o.thickness_px - lacey_2d_projection;
        else
            lacey_thickness_map=0;lacey_2d_projection = 0;h2o_2d_projection_ss=0;
            m.lacey_supp=0;p.e.Io_lacey=0;
        end        
        
        p.h2o.n    = p.h2o.na * p.e.px_size^3 * h2o_2d_projection * (f1.h2o + 1i*f2.h2o);%cell
        p.h2o.n_ss = p.h2o.na * p.e.px_size^3 * h2o_2d_projection_ss * (f1.h2o + 1i*f2.h2o);
        m.h2o = p.h2o.n;
        m.h2o_ss = p.h2o.n_ss;
        
        p.cell.n = p.cell.na * p.e.px_size^3 * cell_2d_projection * (f1.all + 1i*f2.all);
        m.cell = p.cell.n;% n is electron density map of object in e/um2

        p.lacey.n = p.lacey.na * p.e.px_size^3 * lacey_2d_projection * (f1.au + 1i*f2.au);
        m.lacey = p.lacey.n;
        

        %%% calculate dose imparted on cell
        eval(sprintf('protein_AbCoef=interpolate(sc.protein.mu(:,1)*1e3,sc.protein.mu(:,2),p.e.energy);'));
        miu_cell =  (4*pi./p.e.lambda) *(p.e.re*p.e.lambda.^2)./(2*pi) * ( p.cell.na.*f2.all ); %um

        p.e.dose = p.e.true_flux_ds * p.e.energy * miu_cell / p.cell.density * 1e3 * (1e4)^3 * 1.6e-19; % not complete. see Rodriguez 2015
        % dose = flux * mu*E/rho, where mu is linear absorption
        % coefficient, E is photon energy, and rho is cell density. mu/rho
        % is the mass absorption coefficient. Last 3 terms are unit
        % conversion from g to kg, um to cm, and from eV to J.

        %%% generate diffraction pattern
        % scale intensities going through each pinhole
        flux_scaling = ((p.e.re*p.e.lambda) / (p.e.os*p.e.diameter)); % scaling
        p.cell.scaling = sqrt(p.e.true_flux_ds * p.e.QE) * flux_scaling; % scale flux to cell
        p.lacey.scaling = sqrt(p.e.true_flux_ss * p.e.QE) * flux_scaling; % scale flux to static structure

        %%%%%%%%%%%%%%%%%%%%%%% model  %%%%%%%%%%%%%%%%%%%%%%%
        ss = (m.h2o_ss+m.lacey) .* m.lacey_supp .*p.lacey.scaling;
        dd = (m.h2o + m.cell)   .* m.cell_supp.*p.cell.scaling ;
        m.truth = ss + dd;
        m.supp = m.cell_supp | m.lacey_supp;        
        model_truth_frames{time_id} = m.truth;

        %%%%%%%%%%%%%%%%%%%%%%% CDI--- exit wave & diffraction pattern  %%%%%%%%%%%%%%%%%%%%%%%
        PO = m.truth;
        kPO = my_fft( PO );
        d.dp = abs(kPO).^2;

        p.e.integrated_counts = sum(sum(d.dp));
        % applied correct scaling factor so total photon counts on the detector is consistent with total flux
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
        dp_noisy{time_id} = d.dp_noisy;
        d.dp(dc == 1 ) = -1;

        p.e.dp_rfac = krfac(d.dp, d.dp_noisy.^2, beamstop);
        d.rfac=krfac(d.dp_noisy, sqrt(d.dp), beamstop); % calculate R-factor of noisy pattern

        %%%%%%%%%%%%%%%%%%%%%%% Perfect Lens  %%%%%%%%%%%%%%%%%%%%%%%
        na_f_z = p.e.px_size* ( p.h2o.na.*f1.h2o.*h2o_2d_projection + p.cell.na.*f1.all.*cell_2d_projection );
        delta_z = (p.e.re*p.e.lambda.^2)./(2*pi) * na_f_z; %um        
        phi = 2*pi./p.e.lambda *delta_z ;

        phi0 = min( phi(phi~=0));
        phimax = max(phi(:));
        phasecontrast(time_id) = max(phi(:))-phi0;
        avephi0 = avephi0 + phi0;
        avephimax = avephimax + phimax;
        

        % amplitude (imag)
        na_f2_z = p.e.px_size* ( p.h2o.na.*f2.h2o.*h2o_2d_projection + p.cell.na.*f2.all.*cell_2d_projection );
        miu_z = (4*pi./p.e.lambda) *(p.e.re*p.e.lambda.^2)./(2*pi) * na_f2_z; %um
       
        absorp = sqrt(exp(-miu_z)); % consider absorption
        % absorp = 1; % phase only
        
        perfect_lens_image = (absorp.^2+2*absorp.*(phi)).*m.cell_supp;
        sumimage = sum(perfect_lens_image(:));
        absorp_inten_ratio = sum(sum(absorp.^2.*m.cell_supp))./sum(sum(m.cell_supp));
        perfect_lens_image = ( absorp_inten_ratio*p.e.QE * p.e.Io_cell/sumimage).*perfect_lens_image;

        perfect_lens_image_noise = poissrnd(perfect_lens_image);
        per_image_noise{time_id} = perfect_lens_image_noise;


        im_noise = per_image_noise{time_id};
        im_noise1 = im_noise(im_noise>0);
        max_im1 = max(im_noise1(:));
        min_im1 = min(im_noise1(:));
        
        im_noise00 = perfect_lens_image;
        im_noise0 = im_noise00(im_noise00>0);
        max_im0 = max(im_noise0(:));
        min_im0 = min(im_noise0(:));

        figure(100); img(m.supp,'m.supp',m.truth,'m.truth',per_image_noise{time_id},'');%,log(1+d.dp_noisy),'');
        caxis([min_im1,max_im1]); colormap(jet); drawnow

        if time_id==4;
        fprintf('cell thickness: %.2f nm \nlacey thickness: %.2f nm \nH2O thickness: %.2f nm\n',p.cell.thickness_um*1e3,1e3*p.lacey.thickness_um,1e3*p.h2o.thickness_um)
        fprintf('phi: min=%.2f, max= %.2f, contras = %.2f.\n',min( phi(phi~=0)), max(phi(:)),phasecontrast);
        fprintf('dose = %.2s, flux in dp = %.2s, total dose = %.2s, dose per um2 = %.2s\n ',...
            p.e.dose,p.e.integrated_counts,p.e.Io_cell, p.e.Io_cell./(p.e.area));
        end

    end

    initial = m.truth/max(m.truth(:));
    diff_pats = zeros(sz,sz,20);
    for i=1:size(dp_noisy,2)
        diff_pats(:,:,i) = dp_noisy{i};
    end
    
    figure();img(log(1+diff_pats),'dp_noisy',model_truth_frames,'model_truth_frames',per_image_noise,'per_image_noise')
    p.e.avephi0 = avephi0/20;
    p.e.avephimax = avephimax/20;
    save (['simdp_3.5e',num2str(run_id-1),'_',num2str(thickness_ratio*1000),'nm.mat'],'model_truth_frames','initial'...
        ,'diff_pats','per_image_noise','p')

end
