# SPIRE Spectral Feature Finder Catalogue

_R. Hopwood, I. Valtchanov, N. Marchili, L.D. Spencer, J. Scott, C. Benson, N. Hładczuk, E.T. Polehamption, N. Lu, G. Makiwa, D.A. Naylor, G. Noble, M.J. Griffin_

This document: **HERSCHEL-HSC-TN-2321**

Mar 2018

## Abstract

The SPIRE Spectral Feature Finder Catalogue is a result of an automated run of the Spectral Feature Finder (FF). The FF is designed to extract significant spectral features from FTS data products. Only High Resolution (HR) sparse or mapping observations produce spectral features, while for Low Resolution (LR) sparse or mapping observations the FF only provides the best fit continuum parameters. The FF engine iteratively searches for peaks over a set of signal-to-noise ratio (SNR) thresholds, either in the HR spectra of the two central detectors (sparse mode) or in each pixel (SPIRE Long Wavelength - SLW, and SPIRE Short Wavelength - SSW) of the two hyper-spectral cubes (mapping). At the end of each iteration the FF simultaneously fits the continuum and the features found. The residual of the fit is used for the next iteration. The final FF catalogue contains emission and absorption features and their respective SNR for each observation, SNR is negative for absorption features. Line fluxes are not included as extracting reliable line flux from the FTS data is a complex process that requires careful evaluation and analysis of the associated spectra. The FTS Spectral Feature Finder Catalogue is 100% complete for features above SNR=10, and 50-70% complete down to SNR=5. The full SPIRE Automated Feature Extraction Catalogue (SAFECAT) contains 166,442 features from 822 sparse and mapping observations, 163,733 of these features are at SNR > 5 with a subset of XXX features with SNR > 10.

## Table of Contents

* [Observations](#ff_obs)
* [FF algorithm](#ff_proc)
    * [Fitting and subtracting the continuum](#step1)
    * [Iterating over SNR thresholds](#step2)
    * [Final SNR estimate and final check](#step3)
    * [FF feature flags](#step4)
    * [Source radial velocity estimate](#step5)
    * [Neutral Carbon Check](#step6)
    * [Bespoke Handling](#step7)
* [FF products](#ff_prod)
    * [Feature catalogue per observation](#cats)
    * [The SPIRE Automated Feature Extraction CATalogue: SAFECAT](#safecat)
    * [Continuum fit parameters](#cont)
    * [The Feature Finder postcards and POSTCATs](#postcats)
* [FF products access](#ff_web)
    * [Folder structure and content](#ff_folders)
    * [Individual FF catalogues and postcards](#ff_wiki)


## <a name="ff_obs"></a>Observations

The observations with FF products are listed in four coma-separated-value (.csv) files:

 - `hrSparseObservations.csv`: lists all the high resolution (HR) sparse-mode observations processed by the FF, including the HR part of H+LR observations. For each observation the file provides the following information: observation ID (`obsid`); source name (`target`); if the source is known to be featureless (`knownFeatureless`), and therefore no FF catalogue is provided; if the source has a significant spatial size (*semiExtended* or fully *extended*) or if it is *pointLike* (`sourceExt`); what data product was used for the FF (`dataUsed`), which can be the standard pipeline product *spg* from the *Herschel* Science Archive, a SPIRE Spectrometer calibration source Highly Processed Data Product *calHpdp* or data corrected for high background or foreground emission *bgs*; if a focused check of <sup>12</sup>`CO(7-6)` and `[CI](2-1)`, i.e., the neutral carbon check (ncc), resulted in one of these features being added to the associated FF catalogue (`nccApplied`); and if any `bespokeTreatment` was needed, such as special parameter settings.<br/>
 There are 868 observations in the file, but note that no FF catalogue entry is provided for cases of no spectral features found within an observation.  There are 822 observations that do have FF catalogue entries for spectral features identified.
 - `hrMappingObservations.csv`: lists the HR mapping observations included in the FF run. The file contains the `obsid`, the target name (`target`), the  spectral resolution (`resolution`), the instrument mode (`instmode`) and the spatial sampling (`mapSampling`). <br/>
 There are 180 observations and all of them have FF catalogue entries.
 - `lrSparseObservations.csv`: lists the obsids of the LR sparse observations. There are 293 observations in total. 
 - `lrMappingObservations.csv`: lists the obsids of the LR mapping observations. There are 106 observations.

## <a name="ff_proc"></a>Feature Finder Algorithm

The _Herschel_ SPIRE Spectral Feature Finder (FF) finding and fitting process is summarised by the [FF flowchart](http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/spectral_feature_catalogue/FF_1stRelease_products/doc/featureFinder_flowchart.pdf) and described in full in _Hopwood et al._ (in preparation). This document briefly describes the main steps of the FF algorithm.

The FTS has two overlapping spectral bands: SSW (191-310 &mu;m or 1568-944 GHz) and SLW (294-671 &mu;m or 1018-447 GHz). They share an overlap region in 294-310 &mu;m (1018-944 GHz). The primary input to the FF script is a single, per band spectrum, that has been extracted from one of the centre detectors (sparse mode) or a single spectrum from each of the two per band hyper-spectral cube spaxels (mapping).

Each feature found is fitted using a sinc-function (`sin(x)/x`) profile of fixed width, with the width set using the actual resolution of the input data.

The following steps are carried out:

### 1. <a name="step1"></a> Fitting and subtracting the continuum

A resampled copy of the input spectrum is shifted by one frequency bin (5 GHz) and subtracted from itself, to look for "jumps" that correspond to significant peaks. The strong peaks are masked in the input spectrum before a 3<sup>rd</sup> order polynomial is fitted (2<sup>nd</sup> order for low resolution observations). Neither the masking nor the peaks are carried forward into the main finding loop, only the polynomial model.

This is the only stage for LR observations (sparse or mapping) for which we only keep the derived continuum parameters.

### 2. <a name="step2"></a> Iterating over SNR thresholds

We use the following SNR thresholds:  +[100, 50, 30, 10, 5, 3] for emission and [-100, -50, -30] for absorption features.

For each threshold, the signal-to-noise ratio (SNR) spectrum is taken using the model subtracted residual and the spectrum dataset "error" column (for the first iteration the continuum model is subtracted).

Peaks are determined by merging all data points that sit above the SNR threshold, within a 10 GHz width per peak.

Each new peak represents a potential new feature, so a sinc function is added to the total model (polynomial + sinc) per new peak found.

A global fit is performed using the input spectrum and total model, so all found and potential features, and the continuum, are simultaneously fitted. 

For each SNR iteration, the resulting fitted sinc models for the potential features found are put through several reliability checks (e.g., checking their fitted position and looking for multiple peaks fitted to partially resolved features).

Features that are accepted as reliable have their fitted position limited to within a 2 GHz window before the global fit is repeated.

The resulting total model is carried forward to the next iteration, with the frequency either-side of each new feature masked by a SNR threshold dependent width, where no new features are permitted in any of the following finding iterations. The existing masking is updated to account for movement of already found features, noting that these have already had the position of their respective sinc models limited.

### 3. <a name="step3"></a>Final SNR estimate and final check

The final SNR is calculated using the fitted peaks and the total-model-subtracted residual spectrum as `[fitted peak]/[local standard deviation]`.

Features with absolute SNR > 5 are carried forward to the final check - a search for fitting to the sinc wings of significant features, while trying to preserve any \[CI\](2-1) detection.

### 4. <a name="step4"></a>FF feature flags

To assess the distinction between false and true spectral features and to identify the likelihood of prospective spectral features being correctly identified, a goodness of fit metric has been developed that combines the goodness of fit for each individual feature found (GoF) with the goodness of fit of the total model (Total Fit Evaluation - ToFE). GoF and ToFE are used in conjunction to assign a flag to each identified spectral feature (FF flag). The relative weighting of the Gof and ToFE in determining the spectral feature flags varies for different frequency regions within the SLW and SSW bands (e.g., the band edges and identified as suffering higher than average noise and thus are flagged using different criteria -- see below).

The Feature Finder requires a goodness of fit metric that is not sensitive to the continuum or to line flux. Such a statistical metric, `r`, is calculated using a cross-correlation function between fitted feature and total model. `r` is not sensitive to the continuum, due to the subtraction of the local mean. `r` is also not sensitive to line flux, because of division by the standard deviation of a given spectral region.

ToFE is a Bayesian method to compare the evidence (also known as model likelihood) for two concurrent models: with and without a particular feature included. It is therefore complimentary to the `r` parameter identified above. For both models the fitting engine calculates the evidence, which is a set of probability logarithms. Their difference provides the odds for the feature, the smaller the odds the better the model with the feature included. As this is a probabilistic check, no assumption about the noise is made and therefore ToFE is insensitive to systematic noise.

#### Flag definitions (metadata keyword)

- `0.0`: good fit in lower noise region (`FLAG_G_G`)
- `0.1`: good fit in noisy region (`FLAG_G_N`)
- `1.0`: poor fit in lower noise region (`FLAG_B_G`)
- `1.1`: poor fit in noisy region  (`FLAG_B_N`)

#### Flagging criteria

The flagging criteria were chosen empirically.

- GoF `r >= 0.64`: good fit
- ToFE `odds <= -6` : good fit

Both criteria must be met for a "good fit", unless there are more than 35 features found in a given spectrum, in which case only GoF is considered, as ToFE tends to return a null result.

#### Noisy regions

The following regions are considered as noisy:

- SLW: &nu; < 608.6 GHz or &nu; > 944.0 GHz
- SSW: &nu; < 1017.8 GHz or &nu; > 1504.7 GHz


### 5. <a name="step5"></a>Source radial velocity estimate

At the end of the Feature Finder (FF) process, for each high resolution (HR) observation, the radial velocity is estimated by searching for <sup>12</sup>`CO` lines and the `[NII]` 205 µm atomic line, in the respective feature catalogue with identified <sup>12</sup>`CO` taking priority over `[NII]`. In addition, a cross-correlation (XCOR) technique is applied using the feature catalogue and a template line catalogue, which includes most of the characteristic molecular and atomic lines in the far-infrared. For FF catalogues where few features have been found, XCOR includes an additional check with `[NII]` (for SSW) and <sup>12</sup>`CO(7-6)` lines. These FF based estimates are compared to radial velocities from [Simbad](http://simbad.u-strasbg.fr/simbad/) and from a collection by the HIFI team (Lisa Benamati, private communication). 

The final radial velocity estimates included in the FF catalogues (for HR observations) are selected from all the available values.

More details on the methods will be provided in *Hopwood et al* and *Hładczuk et al*, both in preparation.

#### Radial velocity metadata and flags

The following radial velocity related metadata are included in each FF catalogue:

- `RV` - the radial velocity in the local standard of rest, in km/s. For sources at high redshift we use the convention that the radial velocity is equal to `c*z`, where `c` is the speed of light and `z` is the redshift: <br/> `1+z` = &nu;<sub>emitted</sub>/&nu;<sub>observed</sub>, where '&nu;' is the frequency.
- `RV_ERR` - the error on the estimate, when available.
- `RV_FLAG` - the flag indicating the quality and the source of the radial velocity. Priority is given to the FF estimate when both XCOR and FF have comparable quality. The flags can be `FF`, `FF?`, `XCOR`, `XCOR?`, `H?`, `S?` or `nan` when unavailable (see descriptions in the list below). The question marks are used as a warning that the value is uncertain.
    * `FF` - confident estimate based on the Feature Finder <sup>12</sup>`CO` or `[NII]` checks, also in agreement with the HIFI team collected radial velocities: either the difference is less than 20 km/s or the fractional difference is within 20%.
    * `FF?` - good estimate based on the Feature Finder <sup>12</sup>`CO` or `[NII]` checks, but either the difference with the HIFI provided value is larger than 20 km/s or the fractional difference is above 20%.
    * `XCOR`- confident estimate based on the cross-correlation method and `[NII]`, <sup>12</sup>`CO(7-6)` checks, also in agreement with the HIFI team collected radial velocities: either the difference is less than 20 km/s or the fractional difference is within 20%.
    * `XCOR?` - good estimate based on the cross-correlation method and `[NII]`, <sup>12</sup>`CO(7-6)` checks, but either the difference with the HIFI provided values is larger than 20 km/s or the fractional difference is above 20%.
    * `H?` - when none of the FF or XCOR methods have reliable estimate then we use the HIFI provided value. As the target names between HIFI and SPIRE may be different, as these are provided by the _Herschel_ observers, we search the HIFI list for the nearest neighbour within 6 arcsec of the SPIRE central detector sky coordinates.
    * `S?` - when there is no HIFI provided value and there are no good estimates by the FF and XCOR we use the radial velocity from [Simbad](http://simbad.u-strasbg.fr/simbad/). We search Simbad within 6 arcsec of the FTS central detector sky coordinates and assign the radial velocity from the nearest neighbour.
    * `nan` - when neither method provide an estimate and no radial velocity information is available from Simbad and HIFI.

*Warning:*

The radial velocities are still preliminary and in many cases, where there are only a few lines or a single feature detected, they may be incorrect.

### 6. <a name="step6"></a>Neutral Carbon Check

The Neutral Carbon Check (NCC) is a focused check of the <sup>12</sup>CO(7-6) and \[CI\](2-1) spectral region using the radial velocity from the previous step and the known positions of these lines. If either one of these neighbouring features were missed by the main FF process, the NCC-identified missing feature is added to the final list of features found.

### 7. <a name="step7"></a>Bespoke Handling

####  HR mapping

For HR mapping observation the initial list of features for the subsequent FF iterations is provided by a python-based peak finder, applied on HR apodized spectra. This method provides better stability and avoids too many spurious features in noisy hyper-spectral cube pixels (spaxels). 

####  Non-standard data

_**Highly Processed Data Products**_

- Highly Processed Data Products (HPDPs) are available for the repeated observations of SPIRE Spectrometer calibration sources presented in [Hopwood et al. (2015) (arXiv:1502.05717)](http://adsabs.harvard.edu/abs/2015MNRAS.449.2274H). These HPDPs are SPIRE spectra that have been corrected for pointing offset and, where necessary, for source extent or high background emission.
- If an HPDP was used instead of the standard _Herschel_ Science Archive product, this is reported under the metadata entry `hpdp` for the individual FF feature catalogues, and in the `HPDP` column of SAFECAT.
- The `Flag` column (in [the FF HR Sparse point-source calibrated product pages](http://herschel.esac.esa.int/twiki/bin/view/Public/HrSparsePage01)) indicates whether a HPDP was used, with a `HPDP` flag.
- The HPDPs themselves can be accessed from the [SPIRE-S calibration targets legacy data page](http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/cal_targets).

_**Background Subtracted (BGS) Spectra**_

- A number of high resolution (HR) observations were visually assessed as suffering from the same "double bump" systematic noise that is corrected in LR data (where it was more extreme and prevalent in a large fraction of LR observations, see Marchili et al. 2016).
- This issues is only seen in a small number of HR observations of faint point-like targets (in total 84).
- The observations were corrected by subtracting a mean sum of smoothed off-axis detectors, as described in Hopwood et al. (2015).
- All BGS data were corrected using the HIPE Background Subtraction useful script in interactive mode.
- If BGS data were used instead of the standard _Herschel_ Science Archive product, this is reported under the metadata entry `bgs` for the individual FF feature catalogues, and in the `BGS` column of SAFECAT.
- The `Flag` column (in the FF HR Sparse point-source calibrated product pages) indicates whether BGS data was used, with a `BGS` flag.
- [The BGS spectra are available here](http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/spectral_feature_catalogue/background_subtracted/).

####  Bespoke FF settings

If, for a particular observations, non-default FF parameters were used or bespoke treatment was applied, this is indicated with a "1" in the `bespokeTreatment` column of [hrSparseObservations.csv](hrSparseObservations.csv), which lists all HR observations for which there is a set of FF products.

_**Fewer negative SNR thresholds**_

By default, the FF iterates over a number of SNR thresholds when looking for peaks: +[100, 50, 30, 10, 5, 3] followed by [-100, -50, -30, -10]. For one particularly spectral rich observation (OBSID: 1342210847), if the default set of negative SNR thresholds is used, there are so many features found that there are not enough unmasked data points available for the final SNR estimates, and thus many features are discarded. This loss of features is prevented by omitting the -30 and -10 SNR threshold iterations.

_**Final SNR estimate**_

To optimise the final SNR estimate, the FF stipulates that the number of unmasked data points in the local region used to calculate the residual standard deviation must be at least 17 (5 GHz). If this condition is not met, then the local region around the fitted peak is widened.

For the observations of two spectral rich sources (OBSIDs: 1342192834 and 1342197466) this condition is never matched for the majority of the features found; using the default FF settings, the majority of features are discarded. For these two cases, no minimum is set for the number of data points required for the the final SNR estimate to go ahead. 

_**Features added by hand**_

During the FF SNR threshold iterations, the SNR is taken using the spectral dataset "error" column. For a handful of observations this can lead to no significant SNR peak at the position of significant spectral peaks and the corresponding features are therefore never found by the FF, regardless of how low the SNR threshold drops.

A handful of missing significant features were added at the appropriate SNR threshold during the FF process for six observations: OBSIDs 1342197466, 1342248242, 1342216879, 1342197466, 1342193670, and 1342210847.

## <a name="ff_prod"></a>The FF products

### <a name="cats"></a>Feature catalogue per observation

The Feature Finder (FF) produces feature catalogues for all SPIRE Spectrometer HR sparse-and mapping mode observations, unless these are of known featureless sources (e.g. Vesta, Ceres) or the target is not included in the FF list of observations (e.g. dark sky and sources with well developed models, including Uranus, Neptune and Mars). Failed observations and calibration observations with unusual settings are also omitted.

The Feature Finder (FF) catalogues are available as FITS files with the catalogue and its metadata in the first Header Data Unit (HDU). 

The **sparse mode** catalogue table contains the following columns:

- `frequency` - the measured frequency, in GHz, of features found with SNR > 5 (noting that negative SNR correspond to absorption features)
- `frequencyError` - the error on the measured frequency, also in GHz
- `SNR` - the SNR measured using the fitted peak and the local noise in the residual spectrum (after the total fitted model has been subtracted)
- `detector` - which detector the feature was found in (SLWC3 or SSWD4 for sparse-mode observations)
- `featureFlag` - a flag to show if the fit is considered good or poor and if the feature is in a high noise or lower noise region of the spectrum. The flags are explained in more detail in section [_FF feature flags_](#step4).

**Note \#1:** the frequency axis of all SPIRE spectra, as well as those from the other _Herschel_ spectrometers, are provided in the [kinematic Local Standard of Rest (LSRk)](http://herschel.esac.esa.int/hcss-doc-15.0/load/hifi_um/html/hifi_ref_frame.html#hum_lsr_frame). Consequently the  measured frequencies in the FF are also in the same LSRk reference frame.

**Note \#2:** no consolidation of features in the overlap region (944-1018 GHz) was performed. Therefore the same feature may be present for both central detectors for the same observation.

The **mapping mode** catalogue table contains the following _additional_ columns:

- `row` - the cube row pixel where the feature is found 
- `column` - the cube column pixel
- `ra` - the RA (in degrees) for the corresponding pixel where the feature is found 
- `dec` - the Declination (in degrees) for the corresponding pixel where the feature is found 
- `array` - instead of the detector name, this column provides the FTS array name, can be `SSW` or `SLW`.

**Note:** the SSW and SLW hyper-spectral cubes have different world-coordinate-systems (WCS), their pixel size and centres are not matched, i.e., `row` and `column` are `array` dependent. This is clearly evident within the SLW and SSW maps presented in the mapping postacrds.

In addition to the catalogue table, the mapping FITS file contains 3 more HDU extensions: `velocity`, `velError` and `velFlag`. These have the same dimensions as the SSW cube and each pixel contains the derived velocity from both SSW and SLW, its error, and flag, respectively. Pixels without velocity estimate are encoded as `NaN`.

The catalogue FITS extension metadata contain the following information:

- `FLAG_*` - feature flag definitions, see section [_FF feature flags_](#step4);
- information on the observation; the observation ID (`OBS_ID`), the source name as given by the observer (`OBJECT`); the commanded source coordinates `RA_NOM` and `DEC_NOM`, the operational day (`ODNUMBER`) when the source was observed; the resolution (`COMBRES`); the bolometer detectors bias (`BIASMODE`); and the map sampling (`MAPSAMPL`);
- some FF input parameters; the minimum SNR cut applied (`MIN_SNR` equal to 5) and the region avoided at the ends of the frequency bands (`EDGE_MASK` equal to 10 GHz);
- the maximum value of the fitted continuum in SSWD4, `MAX_CONT`;
- how many features were found per detector: `N_SSW` and `N_SLW`;
- an estimate of the source radial velocity (`RV`), the error associated to that estimate (`RV_ERR`) and a radial velocity flag (`RV_FLAG`); these are described in section [_Source radial velocity estimate_](#step5).<br/>
For **mapping mode** the radial velocity, its error and flag are kept in separate extensions (see above);
- the source extent (`S_EXTENT`) as classified from assessing the quality of the spectra in comparison to any associated PACS photometer maps (_pointLike_, _semiExtended_ or _extended_). This keyword is not available for **mapping mode**;
- the calibration scheme used (`CAL_TYPE`, can be _pointSource_ or _extended_) and units of the data (`FLXUNIT` which can be Jy or W/m2/Hz/sr);
- The metadata also reports if non-standard spectra were used, i.e. not the SPG product available in the _Herschel_ Science Archive, where either a Highly Processed Data Product (`HPD_USED` = True|False) or spectra that have been corrected for high background or foreground emission (`BGS_USED` = True|False) were used. More information on the HPDPs used by the FF can be found on the [SPIRE-S calibration targets](http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/cal_targets) legacy page. The background subtracted data used for the FF can be found [in the _Herschel_ legacy area](http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/BKGS/). These are only present for **sparse mode**;
- Due to the close proximity of <sup>12</sup>`CO(7-6)` and `[CI](2-1)`, at 806.7 and 809.3 GHz respectively, a focused check is performed at the end of the Feature Finder process. If either of this pair is found to be missed by the iterative FF process, the position of the missed line and an improved estimate of its SNR is added to the catalogue and the metadata parameter `CI_CHECK` is set to 1 (True).

### <a name="safecat"></a>The SPIRE Automated Feature Extraction CATalogue: SAFECAT

The SAFECAT contains all the features found from the catalogues per observation for SPIRE Spectrometer HR observations. SAFECAT is intended as an archive mining tool that can be searched by frequency, position, etc., to provide all SPIRE Spectrometer observations with significant features that match the search criteria. 

SAFECAT is distributed as a FITS table with the catalogue and its metadata in the first HDU. The table contains the following columns:

- the name of the source observed (`object`);
- whether the source is <em>point-like</em>, <em>semi-extended</em> or <em>extended</em>, categorised as explained above (`spatialExtent`);
- the operational day the observation was taken and the observation identifier (`operationalDay`, `obsid`);
- the (row,column) pixel for mapping mode or (0,0) for sparse (`row, column`);
- the RA and Dec of the central detectors for sparse mode or the RA and Dec of the corresponding (row,column) pixel for mapping (`RA, Dec`);
- the bias mode of the observation, which can be nominal or bright (`biasMode`);
- the bolometer detector name (sparse) or array name (mapping) the feature was detected in  (`detector`). <br/>Note that because of the SLW/SSW band overlap region (944-1018 GHz), features found within the overlap may occur twice for any given observation. This is a potential duplication for **sparse mode**, but is not directly applicable for mapping mode as the spectral maps for SSW and SLW have different WCS, so the overlapping SLW/SSW pixels do not correspond to the exact same spatial coverage region of the sky;
- the maximum SLW continuum level (`maxContinuum`);
- the fitted feature position in GHz and the error on this measurement (`frequency, frequency error`), in LSRk frame.
- the signal-to-noise ratio of fitted peak to local noise in the full residual (`SNR`);
- the feature flag, as explained above (`featureFlag`);
- the source radial velocity, in km/s, estimated from the Feature Finder results or the associated HIFI estimate (which comes from published literature) or from SIMBAD;
- the radial velocity uncertainty in km/s and associated quality flag (`radVelErr, radVelFlag`);
- whether a highly processed data product or background corrected data product was used (`HPDP, BGS`);
- and whether a focused check of <sup>12</sup>`CO(7-6)` and `[CI](2-1)` resulted in one of these being added to the catalogue after the main script was run (`ciCheck`);

The metadata provides the feature flag definitions; the minimum SNR cut applied (5); the frequency range avoided at the ends of the bands (10 GHz); and lists two special calibration observations and the unique IDs assigned for the purpose of SAFECAT only, as they consist of two sparse pointings in one observations.

### <a name="cont"></a> Continuum fit parameters

The SPIRE FTS instrumental line shape is essentially a sinc-function (`sin(x)/x`). The sinc-like wings of each feature introduce ringing, which although decreasing in amplitude, does extend over the whole frequency range. Therefore, to gain a good fit to the continuum, this should be simultaneously fitted with the main spectral features.

During the Feature Finder finding and fitting process, a 3<sup>rd</sup> order polynomial (2<sup>nd</sup> order for LR) is fitted to the continuum of each _per band_ spectrum in conjunction with sinc functions for each feature found. The resulting continuum fit may be hard to precisely recreate, unless a similar procedure is carried out. Therefore the best fit polynomial parameters are provided with the other Feature Finder products.

The parameters are of the form:<br/> `p0 + p1*freq + p2*freq^2 + p3*freq^3` (or `p0 + p1*freq + p2*freq^2` for LR).

The frequency ranges for the SPIRE FTS, used in the continuum calculations, are [446.99, 1017.80] GHz for SLW and [944.05, 1567.91] GHz for SSW.

For a given observation, the associated _fittedContinuumParameters_ FITS file contains a table in the first HDU with the parameters for each detector that has been through the Feature Finder algorithm: the centre detectors for sparse observations, and each pixel for mapping observations.

For **sparse mode** the fitted polynomial is also provided as a _Herschel_ spectrum dataset, with the best fit parameters reported in the metadata. The dataset is stored in a FITS file with each detector in a separate HDU, identified by the detector name, e.g. hdu["SSWD4"] will contain the continuum spectrum for the central detector of the SSW bolometer array. The frequency grid in these spectrum datasets is the same as the original spectra, so these can be used directly to subtract the continuum.

## <a name="postcats"></a> The Feature Finder postcards and POSTCATs

### Sparse mode

The Feature Finder results are visually summarised per observation in a “postcard”. Each postcard compares the input spectra for the centre detectors, the best fit to the continuum, and includes vertical ticks representing the features found (with these arbitrarily scaled by the signal-to-noise). By default, the Feature Finder operates on point-source calibrated data, which has flux density units of Jy. For sources with some spatial extent there is also a postcard for the extended-source calibrated data, which has surface brightness units of W/m2/Hz/sr. Both postcards for such sources are included in the product pages, accessible via the tables at the end of the [Feature Finder Release Note](http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/spectral_feature_catalogue/README.html). 

If an observation of interest is of a partially or fully extended source, it is advisable to visually check both sets of FF results.

For LR observations no features are provided, as LR was primarily aimed at providing a measurement of the continuum and there are few features found in these observations. However, two postcards are provided for LR, as many targets are semi-extended in nature and simply by comparing the postcards it may be possible to gauge which calibration scheme is more appropriate (or if there is a problem that may require further processing to correct for partial extent). Visual inspection should evaluate the continuity of the SLW/SSW spectra within the spectral overlap region.

All HR and LR postcards are collected in the the [POSTcard CATalogues (POSTCATs)](http://archives.esac.esa.int/hsa/legacy/HPDP/SPIRE/SPIRE-S/spectral_feature_catalogue/FF_1stRelease_products/POSTcardCAT/). These are provided as PDF files. There are POSTCATs that compare the point-source and extended-source calibrated results. For LR there is also a mapping POSTCAT that compares the results for the two types of hyper-spectral cubes available in the <em>Herschel</em> Science Archive.

### Mapping mode

The Feature Finder mapping results are also visually summarised per observation in a “postcard”. Each postcard is comprised of a 2x3 array of figures illustrating various aspects of the observation.  The left column illustrates the integrated flux associated with each band, with SLW on top and SSW on the bottom. Note the difference in pixel sizes for the SLW and SSW maps. Each of the integrated flux maps in the left column have two pixels identified for each band with coloured box outlines. For the SSW array, these correspond to the dimmest and brightest pixels within the map.  For the SLW array these pixels correspond to the closest SLW pixels to the identified brightest and dimmest pixels for the SSW array. The central column presents the spectra corresponding to the flagged pixels in the flux maps (brightest and dimmest integrated flux pixels for SSW, closet match for the SLW array) with SLW on top and SSW on the bottom. Also shown in the central column is the spectrum corresponding to the pixel with the most spectral features identified by the FF, again with SLW on top and SSW on the bottom. The upper right figure is a map of the number of lines identified by the FF in both the SLW and SSW arrays combined. This figure also has the pixel regions associated with the most lines in the SLW and SSW arrays identified with coloured box outlines. As the SLW pixels are larger than the SSW pixels, the SLW lines may be counted in multiple SSW pixels within this FF histogram map.  The lower right figure presents the FF radial velocities provided by the radial velocity routine. The results are shown on SSW pixels, with the pixel colour indicating the velocity (colourbar on the side), and a code within each SSW pixel indicating the quality of the radial velocity. Within the radial velocity map, a marker of `A' means all of the expected CO lines contributed to the radial velocity estimate, `N' means that the [NII] nitrogen spectral feature was used, and a number indicates the number of CO lines used in the radial velocity estimate (less than all of the expected lines were identified).  The header of the mapping postcard indicates the target name, as provided by the observer, and the observation number (OBSID).  Pixels outlined in green on the velocity map highlight pixels where the feature finder did find spectral features but the radial velocity routine did not provide an accurate velocity estimate. 

# <a name="ff_web"></a>FF products access 

The FF products are available as Highly Processed Data Products (HPDP) in the *Herschel* legacy area and also in the *Herschel* Science Archive (HSA). 

The HSA is providing access to the individual FF products linked to a particular observation ID, while the HPDP gives access to all products within a browser window (web server).

All products are combined in a single tar.gz file: [FF\_v2.tar.gz](../FF_v2.tar.gz).

## <a name="ff_folders"></a>Folder structure and content

**TBW**

_the text below is obsolete_


At the top level, there are a number of folders and CSV files
([hrSparseObservations\_FF\_1stRelease.csv](http://archives.esac.esa.int
/hsa/legacy/HPDP/SPIRE/SPIRE-S/spectral_feature_catalogue/
FF_1stRelease_products/hrSparseObservations_FF_1stRelease.csv)), 

For each observation the file provides the following information:
observation ID (`obsid`); source name (`target`); if the source is known
to be featureless (`knownFeatureless`), and therefore no FF catalogue is
provided; if the source has a significant spatial size (*semiExtended*
or *fullyExtended*) or if it is *pointLike* (`sourceExt`); what data
product was used for the FF (`dataUsed`), which can be the standard
pipeline product *spg* from the *Herschel* Science Archive, a SPIRE
Spectrometer calibration source Highly Processed Data Product *calHpdp*
or data corrected for high background or foreground emission *bgs*; if a
focused check of <sup>12</sup>`CO(7-6)` and `[CI](2-1)` resulted in one
of these being added to the associated FF catalogue (`nccApplied`); and
if any `bespokeTreatment` was needed, such as special parameter
settings.

In the [*doc*](doc) folder are a number of *readMe* files that detail
the FF products and their content:

 - [FF Products](doc/readMe_products.html)
 - [FF feature flags definitions](doc/readMe_flagDefinitions.html)
 - [FF algorithm details](doc/readMe_FF_algorithm.html) and the [FF
flowchart](doc/featureFinder_flowchart.pdf)
 - [FF bespoke treatement](doc/readMe_bespokeHandling.html)
 - [Radial velocity metadata](doc/readMe_radialVelocity.html)

The [*scripts*](scripts) folder contains a script to recreate the
postcard for any given FF result. The script needs
[HIPE](http://herschel.esac.esa.int/hipe/).

There are five other folders:

 - [“HRpointProducts”](HRpointProducts), which contains the products for
all HR sparse observations processed with the FF, using the point-source
calibrated spectra.
 - [“HRextProducts”](HRextProducts), which contains FF products for
extended-source calibrated data of HR sparse observations, but only for
sources categorised as not point-like.

These product folders have sub-directories containing the FF postcards,
the best fit continuum parameters, the continuum in a spectral dataset,
and the FF feature catalogues.

 - [“LRproducts”](LRproducts) contains all low resolution (LR) FF
products for point-source and extended-source calibrated spectra and for
LR mapping observations (LRpoint, LRext and LRmap, respectively). These
also include the LR part of H+LR observations. LR FF products include
the best fit continuum parameters and FF postcards.
 - [“POSTcardCAT”](POSTcardCAT) contain four postcard catalogues
(POSTCATs), which visually summarise the FF results per observation.
There are POSTCATs that collate all HR point-source calibrated postcards
([POSTCAT\_HR\_Sparse.pdf](POSTcardCAT/POSTCAT_HR_Sparse.pdf)),
point-source and extended-source calibrated postcards for LR and HR data
(suffix “pointVsExt”) and the postcards for LR mapping
([POSTCAT\_LR\_mapping.pdf](POSTcardCAT/POSTCAT_LR_mapping_cpVsNaive.pdf
)).
 -
[“SpireAutomatedFeatureExtractionCATalogue”](
SpireAutomatedFeatureExtractionCATalogue/SAFECAT_v1.fits.gz) contains
the combined HR sparse FF significant feature catalogue: SAFECAT.

## <a name="ff_wiki"></a>Individual FF catalogues and postcards

The FF catalogues can be accessed via a number of wiki pages detailed in the tables below. Each table provides links to the FF product tables in the left-hand column. The observations included on each page is given in the middle column, with the corresponding operational days in the right-hand column.

<table align="center" border="1" cellpadding="1" cellspacing="1" style="width:500px;">
    <thead>
        <tr>
            <th scope="col">Go to page</th>
            <th scope="col">Observations covered</th>
            <th scope="col">Operational days</th>
        </tr>
    </thead>
    <tbody>
        <tr><td colspan="3" align="center">HR sparse-mode FF products</td>
</tr>
            <td style="text-align: center;">
            <a href="http://herschel.esac.esa.int/twiki/bin/view/Public/HrSparsePage01"
target="_blank"> <span style="color:#0000FF;">HR Sparse Page
1</span></a></td>
            <td style="text-align: center;">1342187893 — 1342212341</td>
            <td style="text-align: center;">209 — 602</td>
        </tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/HrSparsePage02"
target="_blank"><span style="color: rgb(0, 0, 255); text-align: center;
background-color: rgb(255, 255, 255);">HR Sparse Page 2</span></a></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1342212342 —
1342231985</span></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">602 — 908</span></td>
        </tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/HrSparsePage03"
target="_blank"><span style="color: rgb(0, 0, 255); text-align: center;
background-color: rgb(255, 255, 255);">HR Sparse Page 3</span></a></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1342231986 —
1342247763</span></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">908 — 1151</span></td>
        </tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/HrSparsePage04"
target="_blank"><span style="color: rgb(0, 0, 255); text-align: center;
background-color: rgb(255, 255, 255);">HR Sparse Page 4</span></a></td>
            <td style="text-align: center;">1342247764 — 1342258698</td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1151 — 1335</span></td>
        </tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/HrSparsePage05"
target="_blank"><span style="color: rgb(0, 0, 255); text-align: center;
background-color: rgb(255, 255, 255);">HR Sparse Page 5</span></a></td>
            <td style="text-align: center;">1342258699 — 1342270195</td>
            <td style="text-align: center;">1335 — 1434</td>
        </tr>
<tr><td colspan="3" align="center">HR mapping-mode FF products</td>
</tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/HrMappingPage01"
target="_blank"><span style="color:#800080;">HR Mapping Page 01
</span></a></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1342192173 — 1342245117
</span></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">302 — 1079</span></td>
        </tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/HrMappingPage02"
target="_blank"><span style="color:#800080;">HR Mapping Page 02
</span></a></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1342245083 — 1342270045
</span></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1080 — 1433</span></td>
        </tr>

<tr><td colspan="3" align="center">LR sparse-mode FF products</td>
</tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/LrSparsePage01"
target="_blank"><span style="color:#800080;">LR Sparse Page
1</span></a></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1342188674 —
1342248229</span></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">227 — 1160</span></td>
        </tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/LrSparsePage02"
target="_blank"><span style="color: rgb(128, 0, 128); text-align:
center; background-color: rgb(255, 255, 255);">LR Sparse Page
2</span></a></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1342248245 —
1342257934</span></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1160 — 1326</span></td>
        </tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/LrSparsePage03"
target="_blank"><span style="color: rgb(128, 0, 128); text-align:
center; background-color: rgb(255, 255, 255);">LR Sparse Page
3</span></a></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1342259570 —
1342270194</span></td>
            <td style="text-align: center;"><span style="text-align:
center; background-color: rgb(255, 255, 255);">1340 — 1434</span></td>
        </tr>
        <tr>
	<td colspan="3" align="center">LR intermediate and fully sampled mapping mode FF products</td>
</tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/LrMappingPage01
" target="_blank"><span style="color: rgb(128, 0, 128);">LR Mapping Page
1</span></a></td>
            <td style="text-align: center;">1342192179 — 1342262926</td>
            <td style="text-align: center;">302 — 1362</td>
        </tr>
        <tr>
            <td style="text-align: center;"><a
href="http://herschel.esac.esa.int/twiki/bin/view/Public/LrMappingPage02
" target="_blank"><span style="color: rgb(128, 0, 128);">LR Mapping Page
2</span></a></td>
            <td style="text-align: center;">1342262927 — 1342270038</td>
            <td style="text-align: center;">1362 — 1433</td>
        </tr>
    </tbody>
</table>

---
Ivan Valtchanov, HSC, 28 Mar 2018

