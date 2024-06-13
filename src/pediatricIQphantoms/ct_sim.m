function res = ct_sim(phantom, patient_diameter, reference_diameter, relative_lesion_diameter, I0, nb, na, ds, sdd, sid, offset_s, down, has_bowtie, add_noise, aec_on, nx, fov, fbp_kernel, nsims)
    % Run a CT simulation for a given `phantom` object
    %
    % :param phantom: phantom object to be scanned, options include ['CCT189', 'CTP404', 'Uniform']
    % :param patient_diameter: effective diameter in mm.
    % :param reference_diameter: Optional, reference effective diameter in mm for computing automatic exposure control (aec) noise index. 
    %        For example if a 200 mm reference phantom has a noise level of 24 HU at I0=3e5, smaller phantoms will scale I0 to match that noise level.
    %   Note this only applies if 'aec_on`=True.
    % :param relative_lesion_diameter: bool | list | float | int, if False lesions scale with phantom size,
    %        if model=='CTP404' lesion_diameter > 1 interpret as absolute diameter in mm,
    %        if lesion_diameter < 1 interpret as scale relative to phantom diameter.
    %        If model=='CCT189' a list of 4 diameters must be provided for the four inserts in contrast order 14, 7, 5, 3 HU.
    %        (only applies for CCCT189 MITA and CTP404 phantoms)
    % :param I0: float, fluence at the detector in the projection data for determining quantum noise
    % The following parameters all belong to MIRT's `sino_geom`_, please see the MIRT documentation for more details
    % :param nb: int, number of detector columns
    % :param na: int, number of angular views in a rotation, 
    % :param ds: detector column size in mm
    % :param sdd: source-to-isocenter distance in mm
    % :param sid: source-to-detector distance in mm
    % :param offset_s: float, lateral shift of detector (1.25 = quarter pixel offset)
    % :param down: downsampling, defaults to 1 but can be increased for faster run times for testing purposes
    % Non-MIRT parameters
    % :param has_bowtie: whether to add a patient fitting bowtie 
    % :param aec_on: 'aec' = automatic exposure control, when `true`, it ensures constant noise levels for all
    %        `patient_diameters` (see `reference_dose_level` for more info)
    % :param nx: reconstructed matrix size in pixels (square, equal on both sides)
    % :param fov: float, reconstructed field of view (FOV) units mm
    % :param fbp_kernel: str, `hanning,xxx`, xxx = the cutoff frequency
    %        see `fbp2_window.m in MIRT <https://github.com/JeffFessler/mirt/blob/main/fbp/fbp2_window.m>`_ for details.
    %        E.g. 'hanning,2.05' approximates a sharp kernel D45 in Siemens Force and 'hanning, 0.85' approximates a smooth kernel B30.
    % :param nsims: int, number of simulations to perform with different noise instantiations
    % .. sino_geom: https://github.com/JeffFessler/mirt/blob/d5685692254247ecac42e6e7dbef328af0a812d5/fbp/sino_geom.m#L17>`_
  
    warning('off', 'all');
    if ~exist('mirt-main', 'dir')
        unzip('https://github.com/JeffFessler/mirt/archive/refs/heads/main.zip', '.');
    end
    run('mirt-main/setup.m')
 
    if iscell(reference_diameter)
        reference_diameter = cell2mat(reference_diameter);
        relative_lesion_diameter = cell2mat(relative_lesion_diameter);
        nb = cell2mat(nb);
        na = cell2mat(na);
        ds = cell2mat(ds);
        sdd = cell2mat(sdd);
        sid = cell2mat(sid);
        offset_s = cell2mat(offset_s);
        down = cell2mat(down);
        has_bowtie = cell2mat(has_bowtie);
        add_noise = cell2mat(add_noise);
        aec_on = cell2mat(aec_on);
        nx = cell2mat(nx);
        nsims = cell2mat(nsims);
    end
    dod = sdd - sid;
    sg = sino_geom('fan', 'units', 'mm', ...
    'nb', nb, 'na', na, 'ds', ds, ...
    'dsd', sdd, 'dod', dod, 'offset_s', offset_s, ...
    'down', down);

    ig = image_geom('nx', nx, 'fov', fov, 'down', down);

    mu_water = 0.2059 / 10;     % in mm-1
    aec_factor = exp(mu_water*patient_diameter)./exp(mu_water*reference_diameter);

    % relative_lesion_diameter = 0.01335; <-- TODO add switch for this in config file
    % relative_lesion_location = 0.4;
    switch lower(phantom)
        case 'uniform'
            ellipse_obj = CCT189(patient_diameter, mu_water, relative_lesion_diameter);
            ellipse_obj = ellipse_obj(1, :);
        case 'cct189'
            ellipse_obj = CCT189(patient_diameter, mu_water, relative_lesion_diameter);
        case 'ctp404'
            relative_lesion_diameter = 0.08;
            relative_lesion_location = 0.38;
            ellipse_obj = CTP404(patient_diameter, mu_water, relative_lesion_diameter);
    end
    ground_truth = ellipse_im(ig, ellipse_obj, 'oversample', 4, 'rot', 0);
    ground_truth_hu = 1000*(ground_truth - mu_water)/mu_water;

    pathlength = ellipse_sino(sg, ellipse_obj, 'oversample', 4);
    % FBP reconstruction operator
    fg = fbp2(sg, ig,'type','std:mat'); %choose 'std:mat' to be able to using different recon filter
                                        %default would be 'std:mex' but only ramp filter was implemented in it
    if aec_on
        I0 = aec_factor*I0; %accounts for different patient size
    end

    if has_bowtie
        I0_afterbowtie=apply_bowtie_filter(I0, sg, mu_water, patient_diameter);           
    else
        I0_afterbowtie=I0;            
    end

    noiseless_sinogram = I0_afterbowtie .* exp(-pathlength);

    ny = nx;
    vol = zeros(nsims, nx, ny);

    for sim_idx = 1:nsims
        disp(sprintf('%s, simulation: [%d/%d]', mfilename, sim_idx, nsims))
        if add_noise   
            disk_proj = poisson(noiseless_sinogram); %This poisson generator respond to the seed number setby rand('sate',x');
        else
            disk_proj = noiseless_sinogram;
        end
        disk_proj = replace_zeros(disk_proj);
        sinogram = -log(disk_proj ./ I0_afterbowtie);            % noisy fan-beam sinogram
        mu_image = fbp2(sinogram, fg, 'window', fbp_kernel);
        recon = 1000*(mu_image - mu_water)/mu_water;
        vol(sim_idx, :, :) = recon;
    end
    res.recon = vol;
    res.ground_truth = ground_truth_hu;
    res.sinogram_noiseless = noiseless_sinogram;
end