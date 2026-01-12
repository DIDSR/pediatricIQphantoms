"""Functions and classes for interactively running CT simulations

This module provides classes and functions to interactively run CT simulations and make phantom datasets
"""

from pathlib import Path
import argparse
import tomli
import numpy as np
import os

# from oct2py import octave
import pandas as pd
import pydicom
from datetime import datetime

URL = 'https://zenodo.org/records/11267694/files/pediatricIQphantoms.zip'

class CTobj():
    """
        A class to hold CT simulation data and run simulations

        :param phantom: phantom object to be scanned, options include ['CCT189', 'CTP404', 'Uniform']
        :param patient_diameter: Optional, effective diameter in mm.
            See `AAPM TG220 <https://www.aapm.org/pubs/reports/detail.asp?docid=146>`_ for more details.
        :param reference_diameter: Optional, reference effective diameter in mm for computing automatic exposure control (aec) noise index.
            For example if a 200 mm reference phantom has a noise level of 24 HU at I0=3e5, smaller phantoms will scale I0 to match that noise level.
            Note this only applies if 'aec_on`=True.
        :param I0: Optional float, fluence at the detector in the projection data for determining quantum noise default 3e5 units of photons.
        :param ndetectors: Optional int, number of detector columns
        :param nangles: Optional int, number of views in a rotation
            (na=1160 based on `Zeng et al 2015 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4629802/>`_)
        :param detector_size: Optional float, detector column size in mm
        :param sid: Optional float, source-to-isocenter distance in mm
        :param sdd: Optional float, source-to-detector distance in mm
        :param detector_offset: Optional float, lateral shift of detector (1.25 = quarter pixel offset)
        :param has_bowtie: Optional bool, whether to add a patient fitting bowtie (TODO add different bowtie sizes.)
        :param add_noise: Optional bool, if true adds Poisson noise, noise magnitude set by `reference_dose_level`,
            noise texture set by reconstructed field of view (currently fov = 110% patient_diameter)
        :param aec_on: Optional bool, 'aec' = automatic exposure control,
            when `true`, it ensures constant noise levels for all `patient_diameters` (see `reference_dose_level` for more info)
        :param matrix_size: Optional int, reconstructed matrix size in pixels (square, equal on both sides)
        :param fov: Optional float, reconstructed field of view (FOV) units mm
        :param fbp_kernel: Optional str, `hanning,xxx`, xxx = the cutoff frequency,
            see `fbp2_window.m in MIRT <https://github.com/JeffFessler/mirt/blob/main/fbp/fbp2_window.m>`_ for details.
            E.g. 'hanning,2.05' approximates a sharp kernel D45 in Siemens Force and 'hanning, 0.85' approximates a smooth kernel B30.
        :param nsims: Optional int, number of simulations to perform with different noise instantiations
        :param lesion_diameter: Optional bool | list | float | int, if False lesions scale with phantom size,
            if model=='CTP404' lesion_diameter > 1 interpret as absolute diameter in mm,
            if lesion_diameter < 1 interpret as scale relative to phantom diameter.
            If model=='CCT189' a list of 4 diameters must be provided for the four inserts in contrast order 14, 7, 5, 3 HU.
            (only applies for CCCT189 MITA and CTP404 phantoms)
        :param framework: Optional, CT simulation framework options include `['MIRT'] <https://github.com/JeffFessler/mirt>`_
        :returns: None
    """
    def __init__(self, phantom='CCT189', patient_diameter=200, reference_diameter=200,
                 reference_fov=340, I0=3e5, ndetectors=900, nangles=580, detector_size=1, sid=595, sdd=1085.6,
                 detector_offset=1.25, down_sampling=1, has_bowtie=False, add_noise=True, aec_on=True,
                 matrix_size=512, fov=340, fbp_kernel='hanning,2.05', nsims=1, lesion_diameter=False, age=0,
                 patientname="", patientid=0, studyname="", studyid=0, seriesname="", seriesid=0, framework='MIRT') -> None:
        """Constructor method
        """
        if phantom in ['MITA-LCD', 'MITA', 'MITALCD', 'LCD']: phantom = 'CCT189'
        options = ['CTP404', 'CCT189', 'UNIFORM']
        if phantom.upper() not in options:
            raise ValueError(f'{phantom} not in available phantoms: {options}')
        self.phantom=phantom
        self.patient_diameter=patient_diameter
        self.reference_diameter=reference_diameter
        self.reference_fov=reference_fov
        self.I0=I0
        self.ndetectors=ndetectors
        self.nangles=nangles
        self.detector_size=detector_size
        self.sid=sid
        self.sdd=sdd
        self.detector_offset=detector_offset
        self.down_sampling=down_sampling
        self.has_bowtie=has_bowtie
        self.add_noise=add_noise
        self.aec=aec_on
        self.matrix_size=matrix_size
        self.fov=fov
        self.kernel=fbp_kernel
        self.nsims=nsims
        self.lesion_diameter=lesion_diameter
        self.age=age
        self.patientname=patientname or f'{patient_diameter/10} cm {phantom}'
        self.patientid=patientid
        self.studyname=studyname or self.patientname
        self.studyid=studyid
        self.seriesname=seriesname or f'{self.patientname} I0: {I0}'
        self.seriesid=seriesid
        self.framework=framework
        self.recon=None
        self.projections=None
        self.groundtruth=None

    def __repr__(self) -> str:
        repr = f'{self.__class__} {self.seriesname}'
        if self.recon is None:
            return repr
        repr += f'\nRecon: {self.recon.shape} {self.fov/10} cm FOV'
        if self.projections is None:
            return repr
        repr += f'\nProjections: {self.projections.shape}'
        return repr

    def run(self, verbose=False):
        """
            Runs the CT simulation using the stored parameters.

            :param verbose: optional boolean, if True prints out status updates, if False they are suppressed. 
        """
        phantom=self.phantom
        patient_diameter=self.patient_diameter
        reference_diameter=self.reference_diameter
        reference_fov=self.reference_fov
        lesion_diameter=self.lesion_diameter
        I0=self.I0
        nb=self.ndetectors
        na=self.nangles
        ds=self.detector_size
        sdd=self.sdd
        sid=self.sid
        offset_s=self.detector_offset
        down=self.down_sampling
        has_bowtie=self.has_bowtie
        add_noise=self.add_noise
        aec_on=self.aec
        nx=self.matrix_size
        fov=self.fov
        fbp_kernel=self.kernel
        nsims=self.nsims

        resdict = mirt_sim(phantom=phantom, patient_diameter=patient_diameter, reference_diameter=reference_diameter, reference_fov=reference_fov,
                           I0=I0, nb=nb, na=na, ds=ds, sid=sid, sdd=sdd, offset_s=offset_s, down=down, has_bowtie=has_bowtie,
                           add_noise=add_noise, aec_on=aec_on, nx=nx, fov=fov, fbp_kernel=fbp_kernel, nsims=nsims, lesion_diameter=lesion_diameter, verbose=verbose)
        self.recon = resdict['recon']
        self.projections = resdict['sinogram_noiseless']
        self.groundtruth = resdict['ground_truth']
        return self
    
    def write_to_dicom(self, fname:str|Path, groundtruth=False) -> list[Path]:
        """
            write ct data to DICOM file, returns list of written dicom file names

            :param fname: filename to save image to (preferably with '.dcm` or related extension)
            :param groundtruth: Optional, whether to save the ground truth phantom image (no noise, blurring, or other artifacts).
                If True, 'self.groundtruth` is saved, if False (default) `self.recon` which contains blurring (and noise if 'add_noise`True)
            :returns: list[Path]
        """
        fpath = pydicom.data.get_testdata_file("CT_small.dcm")
        ds = pydicom.dcmread(fpath)
        # update meta info
        ds.Manufacturer = 'Siemens (simulated)'
        ds.ManufacturerModelName = 'Definition AS+ (simulated)'
        time = datetime.now()
        ds.InstanceCreationDate = time.strftime('%Y%m%d')
        ds.InstanceCreationTime = time.strftime('%H%M%S')
        ds.InstitutionName = 'FDA/CDRH/OSEL/DIDSR'
        ds.StudyDate = ds.InstanceCreationDate
        ds.StudyTime = ds.InstanceCreationTime
        
        ds.PatientName = self.patientname
        ds.SeriesNumber = self.seriesid
        
        ds.PatientAge = f'{int(self.age):03d}Y'
        ds.PatientID = f'{int(self.patientid):03d}'
        del(ds.PatientWeight)
        del(ds.ContrastBolusRoute)
        del(ds.ContrastBolusAgent)
        ds.ImageComments = f"effctive diameter [cm]: {self.patient_diameter/10}"
        ds.ScanOptions = 'AXIAL MODE'
        ds.ReconstructionDiameter = self.fov
        ds.ConvolutionKernel ='fbp D45'
        ds.Exposure = self.I0
        
        # load image data
        ds.StudyDescription = f"{self.I0} photons " + self.seriesname + " " + ds.ConvolutionKernel
        if self.recon.ndim == 2: self.recon = self.recon[None]
        nslices, ds.Rows, ds.Columns = self.recon.shape
        assert nslices == self.nsims
        ds.SpacingBetweenSlices = ds.SliceThickness
        ds.DistanceSourceToDetector = self.sdd
        ds.DistanceSourceToPatient = self.sid
        
        ds.PixelSpacing = [self.fov/self.matrix_size, self.fov/self.matrix_size]
        ds.SliceThickness = ds.PixelSpacing[0]

        ds.KVP = 120
        ds.StudyID = str(self.studyid)
        # series instance uid unique for each series
        end = ds.SeriesInstanceUID.split('.')[-1]
        new_end = str(int(end) + self.studyid)
        ds.SeriesInstanceUID = ds.SeriesInstanceUID.replace(end, new_end)
        
        # study instance uid unique for each series
        end = ds.StudyInstanceUID.split('.')[-1]
        new_end = str(int(end) + self.studyid)
        ds.StudyInstanceUID = ds.StudyInstanceUID.replace(end, new_end)
        ds.AcquisitionNumber = self.studyid

        fname = Path(fname)
        fname.parent.mkdir(exist_ok=True, parents=True)
        # saveout slices as individual dicom files
        fnames = []
        vol = self.groundtruth if groundtruth else self.recon
        if vol.ndim == 2: vol = vol[None]
        for slice_idx, array_slice in enumerate(vol):
            ds.InstanceNumber = slice_idx + 1 # image number
            # SOP instance UID changes every slice
            end = ds.SOPInstanceUID.split('.')[-1]
            new_end = str(int(end) + slice_idx + self.studyid + self.seriesid)
            ds.SOPInstanceUID = ds.SOPInstanceUID.replace(end, new_end)
            # MediaStorageSOPInstanceUID changes every slice
            end = ds.file_meta.MediaStorageSOPInstanceUID.split('.')[-1]
            new_end = str(int(end) + slice_idx + self.studyid + self.seriesid)
            ds.file_meta.MediaStorageSOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID.replace(end, new_end)
            # slice location and image position changes every slice
            ds.SliceLocation = self.nsims//2*ds.SliceThickness + slice_idx*ds.SliceThickness
            ds.ImagePositionPatient[-1] = ds.SliceLocation
            ds.ImagePositionPatient[0] = -ds.Rows//2*ds.PixelSpacing[0]
            ds.ImagePositionPatient[1] = -ds.Columns//2*ds.PixelSpacing[1]
            ds.ImagePositionPatient[2] = ds.SliceLocation
            ds.PixelData = (array_slice.copy(order='C').astype('int16') - int(ds.RescaleIntercept)).tobytes()
            dcm_fname = fname.parent / f'{fname.stem}_{slice_idx:03d}{fname.suffix}' if nslices > 1 else fname
            fnames.append(dcm_fname)
            pydicom.dcmwrite(dcm_fname, ds)
        return fnames



def get_octave_session():
    """
    Lazy-load Octave session to avoid import-time hangs.
    """
    try:
        # Force oct2py to use the octave binary in the same environment as python
        # This prevents it from picking up a system octave which might be incompatible (e.g. 10.3.0 vs 8.4.0)
        import os
        import sys
        
        # Try to find system octave first, as the conda env one (10.x) appears broken with oct2py
        # and checking /usr/bin/octave (8.x) showed success.
        
        system_octave = "/usr/bin/octave"
        if os.path.exists(system_octave):
             os.environ['OCTAVE_EXECUTABLE'] = system_octave
             # We don't need to update PATH for oct2py if OCTAVE_EXECUTABLE is set,
             # but it might help if octave calls subprocesses.
             # However, let's trust oct2py to handle the executable.
        else:
             # Fallback to env octave if system not found
             if 'OCTAVE_EXECUTABLE' not in os.environ:
                 bin_dir = os.path.dirname(sys.executable)
                 octave_path = os.path.join(bin_dir, 'octave')
                 if os.path.exists(octave_path):
                     os.environ['OCTAVE_EXECUTABLE'] = octave_path

        
        from oct2py import octave
        # Check if engine is already ready
        if octave._engine:
            return octave
            
        return octave
    except (Exception, KeyboardInterrupt) as e:
        print("\n\033[91mCRITICAL ERROR: Failed to initialize Octave session.\033[0m")
        print("The process might have hung or been interrupted.")
        print("This often happens if there are zombie Octave processes from previous runs.")
        print("Try running the following command to clean up:")
        print("  \033[1mpkill -u $USER octave\033[0m\n")
        
        # Re-raise to ensure the program exits
        if isinstance(e, KeyboardInterrupt):
            raise e
        raise e


def mirt_sim(phantom='CCT189', patient_diameter=200, reference_diameter=200, reference_fov=340,

             I0=3e5, nb=900, na=580, ds=1, sid=595, sdd=1085.6, offset_s=1.25, down=1, has_bowtie=False,
             add_noise=True, aec_on=True, nx=512, fov=340, fbp_kernel='hanning,2.05', nsims=1, lesion_diameter=False, verbose=True) -> dict:
    """
        Python wrapper for calling Michigan Image Reconstruction Toolbox (MIRT) Octave function

        :returns: dict with keys for recon, projection data, ground truth recon.
    """
    if phantom in ['MITA-LCD', 'MITA', 'MITALCD', 'LCD']:
        phantom = 'CCT189'

    if patient_diameter == reference_diameter:
        fov = reference_fov
    else:
        fov = 1.1*patient_diameter
    if isinstance(lesion_diameter, np.ndarray): lesion_diameter = list(lesion_diameter)
    if lesion_diameter:
        if (phantom.upper() == 'CTP404') & (~isinstance(lesion_diameter, int | float)):
            raise ValueError(f'CTP lesion diameter must be length 1 int or float, entered: {lesion_diameter}')
        if (phantom == 'CCT189') & (len(lesion_diameter) != 4):
            raise ValueError(f'CCT189 phantom expects a length 4 list of ints or floats, but {len(lesion_diameter)} != 4, found: {lesion_diameter}')
    curdir = os.path.dirname(os.path.realpath(__file__))
    octave = get_octave_session()
    octave.cd(curdir)
    if verbose:
        return octave.ct_sim(phantom, patient_diameter, reference_diameter,    lesion_diameter, I0, nb, na, ds, sdd, sid, offset_s, down, has_bowtie, add_noise, aec_on, nx, fov, fbp_kernel, nsims)
    return octave.ct_sim_quiet(phantom, patient_diameter, reference_diameter,    lesion_diameter, I0, nb, na, ds, sdd, sid, offset_s, down, has_bowtie, add_noise, aec_on, nx, fov, fbp_kernel, nsims)

adult_waist_circumferences_cm = {
    # 20: 90.7,
    30: 99.9,
    40: 102.8,
    # 50: 103.3,
    60: 106.2,
    70: 106.6,
    80: 104.1
}


def round_to_decile(age): return age - age % 10


def age_to_eff_diameter(age):
    """
        Equation A-3 from TG204 <https://www.aapm.org/pubs/reports/rpt_204.pdf>
    """
    if age > 22:
        return adult_waist_circumferences_cm[round_to_decile(age)]

    x = age
    a = 18.788598
    b = 0.19486455
    c = -1.060056
    d = -7.6244784
    y = a + b*x**1.5 + c *x**0.5 + d*np.exp(-x) 
    eff_diam = y
    return eff_diam


def pediatric_subgroup(diameter):
    if diameter < age_to_eff_diameter(1):
        return 'newborn'
    elif (diameter >= age_to_eff_diameter(1)) & (diameter < age_to_eff_diameter(5)):
        return 'infant'
    elif (diameter >= age_to_eff_diameter(5)) & (diameter < age_to_eff_diameter(12)):
        return 'child'
    elif (diameter >= age_to_eff_diameter(12)) & (diameter < age_to_eff_diameter(22)):
        return 'adolescent'
    else:
        return 'adult'
    
def subgroup_to_age(group):
    if group == 'newborn': return 1/12
    if group == 'infant': return 2
    if group == 'child': return 12
    if group == 'adolescent': return 21
    if group == 'adult': return 39


def dicom_meta_to_dataframe(fname:str|Path) -> pd.DataFrame:
    """
    Takes dicom header information and exports to pandas Dataframe
    """
    fname = Path(fname)
    dcm = pydicom.dcmread(fname)
    subgroup = pediatric_subgroup(18.2)
    age_est = subgroup_to_age(subgroup)
    subgroup, age_est
    phantom = str(dcm.PatientName).split('cm')[1].strip()
    diameter = float(dcm.ImageComments.split(':')[1])
    dose = None
    repeat = 0
    if fname.stem.endswith('noisefree'):
        series = 'noise free'
    elif fname.stem.endswith('groundtruth'):
        series = 'ground truth'
    else:
        series = 'simulation'
        dose = int(str(fname).split('dose_')[1].split('/')[0])
        repeat = int(fname.stem.split('_')[1]) if len(fname.stem.split('_')) > 1 else 0

    return pd.DataFrame({'Name': [str(dcm.PatientName)], 'Patient ID': [dcm.PatientID], 'Study Name': [dcm.StudyDescription], 'Study ID': [dcm.StudyID], 'series': [series], 'effective diameter [cm]': [diameter], 'age [year]': [age_est], 'pediatric subgroup': subgroup,
     'phantom': [phantom], 'scanner': [dcm.Manufacturer.split('(simulated)')[0] + dcm.ManufacturerModelName], 'Dose [%]': [dose],
     'recon': ['fbp'], 'kernel': [dcm.ConvolutionKernel], 'FOV [cm]': [dcm.PixelSpacing[0]*dcm.Rows/10], 'repeat':[repeat], 'file': [fname]})


def run_batch_sim(image_directory: str, model=['MITA-LCD'], diameter=[200], full_dose=3e5, dose_level=[1.0], fbp_kernel='fbp', verbose=True,  **kwargs) -> pd.DataFrame:
    """
        Running simulations in batch mode

        `run_batch_sim` takes lists of parameters (phantoms, diameters, and dose levels) and iterates through all combinations,
        returns a pandas DataFrame summarizing file metadata including filename

        :param image_directory: Directory to save simulated outputs
        :type image_directory: str
        :param model: Optional list, list of phantom models to simulate (can be length 1) options include ['CTP404', 'Uniform', 'MITA-LCD']
        :type model: list[str]
        :param diameter: Optional list, list of simulated phantom diameters in mm (can be length 1)
        :param full_dose: Optional int, units of photons per pixel, multiplied by `dose_level` to determine quantum noise level in projections
        :param dose_level: Optional list, units of photons
        :param verbose: Optional bool, whether to print update information
        :returns: pd.DataFrame

        See `CTobj` for remaining keyword arguments
    """
    image_directory = Path(os.path.abspath(image_directory))
    print(image_directory)
    phantoms = model

    dose_level = full_dose * np.array(dose_level)
    recon = 'fbp'
    reference_fov = 340
    fnames = []
    for phantom_idx, phantom in enumerate(phantoms):
        print(f'{phantom} Simulation series {phantom_idx}/{len(phantoms)}')
        for patientid, patient_diameter in enumerate(diameter):
            patientname = f'{patient_diameter/10} cm {phantom}'
            for studyid, dose in enumerate(dose_level):
                rel_dose = 100*dose/dose_level.max()
                studyname = f'{int(rel_dose)}% dose'
                ct = CTobj(phantom=phantom, patient_diameter=patient_diameter, reference_fov=reference_fov,
                           I0=dose, down_sampling=1, patientname=patientname, patientid=patientid, studyname=studyname, studyid=studyid, seriesname=f'{patientname} {studyname}', **kwargs)
                ct.run(verbose=verbose)
                fname = image_directory / phantom / f'diameter{patient_diameter}mm' / f'dose_{int(rel_dose):03d}' / f"{recon} {fbp_kernel.replace(',','').replace('.','')}" / f'{ct.patientname}.dcm'
                fnames += ct.write_to_dicom(fname)
            # add noise free
            ct.nsims=1
            ct.add_noise=False
            ct.run(verbose=verbose)
            fname = image_directory / phantom / f'diameter{patient_diameter}mm' / f'{ct.patientname}_noisefree.dcm'
            fnames += ct.write_to_dicom(fname)
            # add ground truth
            fname = image_directory / phantom / f'diameter{patient_diameter}mm' / f'{ct.patientname}_groundtruth.dcm'
            fnames += ct.write_to_dicom(fname, groundtruth=True)
    metadata = pd.concat(list(map(dicom_meta_to_dataframe, fnames)), ignore_index=True)
    metadata.to_csv(image_directory / 'metadata.csv', index=False)
    return metadata


def main():
    parser = argparse.ArgumentParser(description='Make Pediatric IQ Phantoms: command line interface')
    parser.add_argument('config', nargs='?', default=None,
                        help="""input is a configuration .toml file containing simulation parameters (see configs/defaults.toml as an example)""")
    args = parser.parse_args()

    if args.config:
        with open(args.config, mode="rb") as fp:
            config = tomli.load(fp)
    else:
        with open("configs/defaults.toml", mode="rb") as fp:
            config = tomli.load(fp)
    
    sim = dict()
    for c in config['simulation']:
        sim.update(c)
        run_batch_sim(**sim)

if __name__ == '__main__':
    main()