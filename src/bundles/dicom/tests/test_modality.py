import pytest

import chimerax.dicom.modality as modality


@pytest.mark.parametrize(
    "modality, expected_values",
    [
        (modality.Autorefraction, "AR"),
        (modality.ContentAssessmentResults, "ASMT"),
        (modality.AudioECG, "AU"),
        (modality.BoneDensitometryUltrasound, "BDUS"),
        (modality.BiomagneticImaging, "BI"),
        (modality.BoneDensitometryXray, "BMD"),
        (modality.ComputedRadiography, "CR"),
        (modality.ComputedTomography, "CT"),
        (modality.CTProtocolPerformed, "CTPROTOCOL"),
        (modality.Diaphanography, "DG"),
        (modality.Document, "DOC"),
        (modality.DigitalRadiography, "DX"),
        (modality.Electrocardiography, "ECG"),
        (modality.CardiacElectrophysiology, "EPS"),
        (modality.Endoscopy, "ES"),
        (modality.Fiducials, "FID"),
        (modality.GeneralMicroscopy, "GM"),
        (modality.HardCopy, "HC"),
        (modality.HemodynamicWaveform, "HD"),
        (modality.IntraoralRadiography, "IO"),
        (modality.IntraocularLensData, "IOL"),
        (modality.IntravascularOpticalCoherenceTomography, "IVOCT"),
        (modality.IntravascularUltrasound, "IVUS"),
        (modality.Keratometry, "KER"),
        (modality.KeyObjectSelection, "KO"),
        (modality.Lensometry, "LEN"),
        (modality.LaserSurfaceScan, "LS"),
        (modality.Mammography, "MG"),
        (modality.MagneticResonance, "MR"),
        (modality.ModelFor3DManufacturing, "M3D"),
        (modality.NuclearMedicine, "NM"),
        (modality.OphthamlicAxialMeasurements, "OAM"),
        (modality.OpticalCoherenceTomographyNonOphthalmic, "OCT"),
        (modality.OphthalmicPhotography, "OP"),
        (modality.OphthalmicMapping, "OPM"),
        (modality.OphthalmicTomography, "OPT"),
        (modality.OphthalmicTomographyBScanVolumeAnalysis, "OPTSBV"),
        (modality.OphthalmicTomographyEnFace, "OPTENF"),
        (modality.OphthalmicVisualField, "OPV"),
        (modality.OpticalSurfaceScan, "OSS"),
        (modality.Other, "OT"),
        (modality.Plan, "PLAN"),
        (modality.PresentationState, "PR"),
        (modality.PositronEmissionTomography, "PT"),
        (modality.PanoramicXRay, "PX"),
        (modality.Registration, "REG"),
        (modality.RespiratoryWaveform, "RESP"),
        (modality.RadioFluoroscopy, "RF"),
        (modality.RadiographicImaging, "RG"),
        (modality.RadiotherapyDose, "RTDOSE"),
        (modality.RadiotherapyImage, "RTIMAGE"),
        (modality.RadiotherapyIntent, "RTINTENT"),
        (modality.RadiotherapyPlan, "RTPLAN"),
        (modality.RadiotherapyRadiation, "RTRAD"),
        (modality.RadiotherapyRecord, "RTRECORD"),
        (modality.RadiotherapySegmentAnnotation, "RTSEGANN"),
        (modality.RadiotherapyStructureSet, "RTSTRUCT"),
        (modality.RealWorldValue, "RWV"),
        (modality.Segmentation, "SEG"),
        (modality.SlideMicroscopy, "SM"),
        (modality.StereometricRelationship, "SMR"),
        (modality.StructuredReport, "SR"),
        (modality.SubjectiveRefraction, "SRF"),
        (modality.AutomatedSlideStainer, "STAIN"),
        (modality.TextureMap, "TEXTUREMAP"),
        (modality.Thermography, "TG"),
        (modality.Ultrasound, "US"),
        (modality.VisualAcuity, "VA"),
        (modality.XRayAngiography, "XA"),
        (modality.ExternalCameraPhotography, "XC")
        # Retired
        ,
        (modality.Angioscopy, "AS"),
        (modality.ColorFlowDoppler, "CD"),
        (modality.Cinefluorography, "CF"),
        (modality.Culposcopy, "CP"),
        (modality.Cytoscopy, "CS"),
        (modality.DuplexDoppler, "DD"),
        (modality.DigitalFluoroscopy, "DF"),
        (modality.DigitalMicroscopy, "DM"),
        (modality.DigitalSubtractionAngiography, "DS"),
        (modality.Echocardiography, "EC"),
        (modality.FluoresceinAngiography, "FA"),
        (modality.Fundoscopy, "FS"),
        (modality.Laparoscopy, "LP"),
        (modality.MagneticResonanceAngiography, "MA"),
        (modality.MagneticResonanceSpectroscopy, "MS"),
        (modality.OphthalmicRefraction, "OPR"),
        (modality.SinglePhotonEmissionComputedTomography, "ST"),
        (modality.Videofluorography, "VF"),
    ],
)
def test_modality_string_equality(modality, expected_values):
    assert modality == expected_values


@pytest.mark.parametrize(
    "modality, expected_values",
    [
        (modality.Autorefraction, False),
        (modality.ContentAssessmentResults, False),
        (modality.AudioECG, False),
        (modality.BoneDensitometryUltrasound, False),
        (modality.BiomagneticImaging, False),
        (modality.BoneDensitometryXray, False),
        (modality.ComputedRadiography, False),
        (modality.ComputedTomography, False),
        (modality.CTProtocolPerformed, False),
        (modality.Diaphanography, False),
        (modality.Document, False),
        (modality.DigitalRadiography, False),
        (modality.Electrocardiography, False),
        (modality.CardiacElectrophysiology, False),
        (modality.Endoscopy, False),
        (modality.Fiducials, False),
        (modality.GeneralMicroscopy, False),
        (modality.HardCopy, False),
        (modality.HemodynamicWaveform, False),
        (modality.IntraoralRadiography, False),
        (modality.IntraocularLensData, False),
        (modality.IntravascularOpticalCoherenceTomography, False),
        (modality.IntravascularUltrasound, False),
        (modality.Keratometry, False),
        (modality.KeyObjectSelection, False),
        (modality.Lensometry, False),
        (modality.LaserSurfaceScan, False),
        (modality.Mammography, False),
        (modality.MagneticResonance, False),
        (modality.ModelFor3DManufacturing, False),
        (modality.NuclearMedicine, False),
        (modality.OphthamlicAxialMeasurements, False),
        (modality.OpticalCoherenceTomographyNonOphthalmic, False),
        (modality.OphthalmicPhotography, False),
        (modality.OphthalmicMapping, False),
        (modality.OphthalmicTomography, False),
        (modality.OphthalmicTomographyBScanVolumeAnalysis, False),
        (modality.OphthalmicTomographyEnFace, False),
        (modality.OphthalmicVisualField, False),
        (modality.OpticalSurfaceScan, False),
        (modality.Other, False),
        (modality.Plan, False),
        (modality.PresentationState, False),
        (modality.PositronEmissionTomography, False),
        (modality.PanoramicXRay, False),
        (modality.Registration, False),
        (modality.RespiratoryWaveform, False),
        (modality.RadioFluoroscopy, False),
        (modality.RadiographicImaging, False),
        (modality.RadiotherapyDose, False),
        (modality.RadiotherapyImage, False),
        (modality.RadiotherapyIntent, False),
        (modality.RadiotherapyPlan, False),
        (modality.RadiotherapyRadiation, False),
        (modality.RadiotherapyRecord, False),
        (modality.RadiotherapySegmentAnnotation, False),
        (modality.RadiotherapyStructureSet, False),
        (modality.RealWorldValue, False),
        (modality.Segmentation, False),
        (modality.SlideMicroscopy, False),
        (modality.StereometricRelationship, False),
        (modality.StructuredReport, False),
        (modality.SubjectiveRefraction, False),
        (modality.AutomatedSlideStainer, False),
        (modality.TextureMap, False),
        (modality.Thermography, False),
        (modality.Ultrasound, False),
        (modality.VisualAcuity, False),
        (modality.XRayAngiography, False),
        (modality.ExternalCameraPhotography, False)
        # Retired
        ,
        (modality.Angioscopy, True),
        (modality.ColorFlowDoppler, True),
        (modality.Cinefluorography, True),
        (modality.Culposcopy, True),
        (modality.Cytoscopy, True),
        (modality.DuplexDoppler, True),
        (modality.DigitalFluoroscopy, True),
        (modality.DigitalMicroscopy, True),
        (modality.DigitalSubtractionAngiography, True),
        (modality.Echocardiography, True),
        (modality.FluoresceinAngiography, True),
        (modality.Fundoscopy, True),
        (modality.Laparoscopy, True),
        (modality.MagneticResonanceAngiography, True),
        (modality.MagneticResonanceSpectroscopy, True),
        (modality.OphthalmicRefraction, True),
        (modality.SinglePhotonEmissionComputedTomography, True),
        (modality.Videofluorography, True),
    ],
)
def test_modality_retired(modality, expected_values):
    assert modality.is_retired == expected_values
