# National Scale Bathymetric Reprocessing

**Author:** Clay Caldwell  
**Contact:** wccaldwe@syr.edu  
**Organization:** Syracuse University  
**Website:** [LinkedIn](https://www.linkedin.com/in/clay-caldwell-9530011a3/)

## Project Overview

- **Problem Statement:** The U.S. Army Corps of Engineers (USACE) currently relies on Triangulated Irregular Networks (TINs) as the default method for generating bathymetric surfaces from hydrographic survey point clouds. While computationally efficient, TINs are prone to artifacts in complex environments and sparsely surveyed, and do not provide any estimate of interpolation uncertainty, potentially limiting the accuracy of dredge volume calculations.

- **Challenge Statement:** Hydrographic surveys vary widely in equipment type (single-beam vs. multi-beam echosounders), point density, spatial anisotropy, and depth variability across the USACE navigation portfolio. No national-scale comparative study has evaluated multiple automated interpolation methods across these diverse survey conditions.

- **Solution Statement:** This study automates and compares six bathymetric interpolation methods — TIN, Natural Neighbor (NN), Inverse Distance Weighting (IDW), Radial Basis Functions (RBF), Ordinary Kriging (OK), and Regression Kriging (RK) — across hundreds of hydrographic surveys from five geographically diverse USACE civil works districts, leveraging automated variogram fitting and cross-validation to eliminate manual parameter tuning.

- **Objective:** To quantify which interpolation method produces the most accurate bathymetric surfaces under varying survey characteristics (anisotropy, depth variability, and point density), and to develop a decision framework to guide USACE practitioners in selecting the most appropriate interpolation method.

- **Literature Review:** Multiple studies have compared interpolation methods for bathymetric surface generation across diverse water bodies globally. TIN, Natural Neighbor, ANUDEM, IDW, and Kriging are the most commonly compared methods (Li et al., 2023). Studies on the Mississippi River found that RBF and anisotropic Ordinary Kriging outperformed other methods in complex environments (Wu et al., 2019), while TIN and a novel Rectilinear IDW approach performed best across 158 SBES surveys in the upper Mississippi River (Andes & Cox, 2017). Similar comparisons in international water bodies have yielded mixed results, often favoring IDW or kriging variants depending on the environment. This study builds on Andes & Cox (2017) by expanding geographic scope to the national scale and addressing their limitation of static, non-optimized kriging parameters through automated variogram fitting.

- **Research Questions:**
  1. Which bathymetric interpolation method (TIN, NN, IDW, RBF, OK, RK) provides the most accurate depth predictions for USACE navigation channel surveys?
  2. How does interpolation method performance vary across different combinations of survey characteristics including levels of spatial anisotropy, depth variability, and point density?
  3. What combinations of survey characteristics can guide practitioners in selecting the most appropriate interpolation method?

## Data Sources

The dataset for this study comes from five USACE civil works districts: Buffalo (LRB), Philadelphia (NAP), Wilmington (SAW), Mobile (SAM), and Portland (NWP). These districts were chosen due to their diverse subaqueous environments and surveying practices. For example, districts exposed to lower energy regimes may see less sediment transport and less shoaling by extension. This could mean that the waterways which are surveyed are less complex and can adequately be surveyed using a simpler single-beam echosounder. The Buffalo district is an example of this, being exposed to no tidal response and lower energy wave spectra. In more dynamic environments like high energy coastal inlets, more shoaling is expected which would lead to more complex bathymetry. A multi-beam echosounder would be best in such cases. The Wilmington district is an example of this as there are numerous high energy, morphologically complex coastal inlets impacted by both ebb and flood shoals. The other districts include various environments spanning low energy, deep, protected ports to high energy, shallow, riparian environments.

| District | District Code | Environment |
|----------|---------------|-------------|
| Buffalo | CELRB | Great Lakes |
| Philadelphia | CENAP | Mid-Atlantic |
| Wilmington | CESAW | Atlantic Coast |
| Mobile | CESAM | Gulf Coast |
| Portland | CENWP | Pacific Northwest |

### Published Data Sources

| Name | Source | Description | Access Method | Data URL | Metadata URL | Details | Data Citation |
|------|--------|-------------|---------------|----------|--------------|---------|---------------|
| eHydro | U.S. Army Corps of Engineers | An enterprise geospatial database used to process and disseminate hydrographic survey data for coastal and inland navigation channels. This data is used to inform navigational charts for waterways across the contiguous U.S., Alaska, Hawaii, and American territories. | ArcGIS REST API or Online | [ArcGIS Web Dashboard](https://www.arcgis.com/apps/dashboards/4b8f2ba307684cf597617bf1b6d2f85d) | [ArcGIS REST Database](https://services7.arcgis.com/n1YM8pTrFmm7L4hs/ArcGIS/rest/services/eHydro_Survey_Data/FeatureServer/0) | Data typically provided in state plane coordinates (US Survey Foot) horizontal projection systems. Vertical projection system typically referenced to local mean lowest low water (MLLW). Point clouds and vectorized survey boundaries provided within an ESRI geodatabase. | USACE. (n.d.). eHydro survey data [Dataset]. ArcGIS REST Services Directory. Retrieved from https://services7.arcgis.com/n1YM8pTrFmm7L4hs/ArcGIS/rest/services/eHydro_Survey_Data/FeatureServer/ |

## Methods

### Data Retrieval
Survey data is downloaded from the USACE eHydro database via the PyGeoHydro Python library, which leverages an ArcGIS REST server to query data by USACE district code and user-defined date search window. SBES surveys are retrieved to test interpolation methods for the historically most common, less dense surveys.

### Survey Characterization
Survey-specific metrics are calculated to characterize each survey prior to interpolation. The Coefficient of Variation (CV) for measured depth, directional anisotropy of depth measurements, and raw data density are computed, as these characteristics are expected to influence interpolation accuracy.

### Interpolation Methods
Six interpolation methods are implemented and compared:

- **TIN (Triangulated Irregular Network):** The current default USACE method. Employs Delaunay triangulation with linear interpolation within triangles. Implemented using SciPy.
- **NN (Natural Neighbor):** Builds on TIN by employing area-weighted interpolation via Voronoi tessellation to reduce artifacts at triangle boundaries. Implemented using SciPy.
- **IDW (Inverse Distance Weighting):** Weights influence of neighboring observations using squared euclidean distance. Implemented using Verde.
- **RBF (Radial Basis Function):** Builds on IDW by implementing a weighting function to define the influence of nearby observations using radial distances. This is applied to each observation and summed to create a single surface. RMSE-based cross-validation is used to determine optimal kernel function and parameters between thin plate splines, cubic polynomial, and quintic polynomials. Implemented using SciPy.
- **OK (Ordinary Kriging):** Models spatial correlation by fitting an empirical variogram to depth data using weighted least squares via GSTools. The optimal variogram shape function (from 9 candidates) is selected by minimizing the Akaike Information Criterion (AIC). Kriging is then executed using PyKrige.
- **RK (Regression Kriging):** Addresses potential nonstationarity by fitting a spline trend surface to depth data using Verde, then applying OK on the spline residuals via PyKrige. Final predictions combine the spline trend and kriged residuals.

The nine variogram shape functions tested include: Spherical, Exponential, Gaussian, Matérn, Stable, Rational, Circular, SuperSpherical, and JBessel (Hole-Effect).

### Output
All interpolation outputs are saved as Cloud Optimized GeoTIFFs (COGs) at 10 ft spatial resolution in the original horizontal and vertical projection system of the input data, bounded by the convex hull survey boundary. For kriging methods, a corresponding kriging variance COG is also saved. The 10 ft resolution is chosen for compatibility with the Corps Shoaling Analysis Tool (CSAT).

### Modeling Framework
- **Cross-validation:** Block 10-fold cross-validation reporting RMSE, MAE, and Mean Error for each method and survey.
- **Residual analysis:** Shapiro-Wilk test, skewness, and kurtosis for normality; Moran's I for spatial independence; Breusch-Pagan test for homoscedasticity.
- **Statistical modeling:** Cross-validated RMSE serves as the response variable in a Linear Mixed Model (LMM) or Generalized LMM (GLMM). Fixed effects include interpolation method, survey characteristics, and their two-way interactions. Survey is treated as a random effect. Post-hoc comparisons using Estimated Marginal Means (EMMs) with Tukey adjustment identify overall and condition-specific best-performing methods. Implemented in R using `lme4` and `emmeans`.

## Repository Structure

```
~/local_data
├── code
│   ├── model_comparisons
│   │   └── # Code and interactive notebooks to calculate and compare interpolation model accuracy using RMSE, MAE, and Mean Error
│   ├── process_data
│   │   └── # Code and interactive notebooks used to process and generate bathymetric surfaces from hydrographic survey point cloud data
│   ├── source_data
│   │   └── # Code and interactive notebooks used to source the data needed, in this case the hydrographic survey point cloud data
│   └── visualize
│       └── # Code to create figures, tables, and maps for the metrics and residuals which are used to identify best interpolation methods
├── figures_tables
│   └── # where the figures and tables presenting the results will be stored
├── model_outputs
│   └── # Where the results of the modeling will be stored. This may include cross-validation metrics in tabular format, bathymetric surfaces in GeoTiff format, etc.
└── src
    └── autodbm
        └── # where the helper functions and classes for the study are located. Functions and classes used to retrieve, process, and assess data for each model
```

### Computational Requirements

This repository is built on [the Pixi package management tool](https://pixi.prefix.dev/latest/). This tool has built-in support for work on multiple platforms (Linux, macOS, Windows, and more), organizes and composes multiple computing environments, manages complex data pipelines via tasks, and has many more features built-in. Additionally, a Dockerfile is provided which will build an x64-based Linux image with Pixi preinstalled. This will allow the user to choose to download Pixi on their machine locally, or spin up a Docker container for better reproducibility.

This repository utilizes the *pixi.toml* manifest to define two Python environments:
- **`process`:** Used for data retrieval, preprocessing, decimation, and interpolation (TIN, NN, IDW, RBF, OK, RK). Key libraries include PyGeoHydro, Verde, GSTools, PyKrige, and SciPy.
- **`eval`:** Used for cross-validation, residual analysis, and statistical modeling. Key libraries include tools for LMM/GLMM modeling and spatial statistics.

- For local installation: x64-based Linux, Windows, and OSX as well as ARM-based OSX are supported. [Pixi local installation steps](https://pixi.sh/latest/installation/)
- For Docker installation: Container will be Linux x64-based, which may lead to slower processing on ARM-based hardware due to additional virtualization layer(s). [Docker installation steps](https://docs.docker.com/engine/install/)

## How to Reproduce

Reproducibility is handled by *Pixi Tasks*, which executes the data pipeline by leveraging the Pixi CLI to run the needed Python files. This is accomplished using the *src/autodbm.reproduce.py* file.

1. Make sure the Pixi manifest and Python environments are established and initialized.
```bash
pixi install -a
```

2. Retrieve eHydro point clouds and generate bathymetric surfaces using the *process* environment. Executing the below from the command line within the repository's parent directory will retrieve the target hydrographic survey data, apply the various interpolation methods, and produce 10 ft spatial resolution bathymetric surface rasters.
```bash
pixi run -e process interpolate
```

3. Evaluate and compare interpolation method accuracy and residuals using the *eval* environment. Executing the below will report the models' validation statistics, assess patterns within model residuals, and output spatial error maps.
```bash
pixi run -e eval analyze
```

4. Run the linear mixed modeling to assess the statistical significance in model accuracy with respect to the survey characteristics. This will be completed using the default environment leveraging the R language.
```bash
pixi run mixedmodel
```

### Data Access

Accessing eHydro data is streamlined using custom-built Python functions. A user can specify a search window (e.g., `01-01-2021` to `01-01-2026`), the USACE civil works district code (e.g., `CESAW` for Wilmington District), or a specific `SurveyId` for a particular hydrographic survey. For this study, search window and USACE district codes were used to collect surveys for testing.

```python
from src.autodbm.processing_help import retrieve_ehydro_data

# Example: retrieve all surveys from Wilmington District between 2021-2026
surveys = retrieve_ehydro_data(
    district_symbol="CESAW",
    start_date="2021-01-01",
    end_date="2026-01-01",
    max_workers = 8
)

# Example: retrieve a particvular surveys from Elizabeth Marine Terminal, Port Newark from 02 MAY, 2023
surveys = retrieve_ehydro_data(
    district_symbol="CENAN",
    surveyId = 'NB_05_PHD_20230502_CS_5289_30',
    max_workers = 8
)
```

## Results

*Results pending. Outputs will include:*
- Bathymetric surface rasters (COGs) for each interpolation method and survey
- Kriging variance surfaces for OK and RK methods
- Block 10-fold cross-validation statistics (RMSE, MAE, Mean Error) by method and survey
- Assessment of model residuals with respect to linearity, spatial independence, normality, and homoscedasticity
- Spatial prediction residual maps
- LMM/GLMM model outputs identifying significant differences between methods
- EMM post-hoc comparison tables
- Decision framework for interpolation method selection based on survey characteristics

## Citation

DOI: ***DOI_PENDING***

## License

*License pending.*

## Contribution Guidelines

*Contribution guidelines pending.*

## References

Andes, L. C., & Cox, A. L. (2017). Rectilinear Inverse Distance Weighting Methodology for Bathymetric Cross-Section Interpolation along the Mississippi River. *Journal of Hydrologic Engineering, 22*(7), 04017014. https://doi.org/10.1061/(ASCE)HE.1943-5584.0001514

Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models using lme4. *Journal of Statistical Software, 67*(1), 1–48. https://doi.org/10.18637/jss.v067.i01

Chegini, T., Li, H.-Y., & Leung, L. R. (2021). HyRiver: Hydroclimate data retriever. *Journal of Open Source Software, 6*(66), 3175. https://doi.org/10.21105/joss.03175

Cressie, N. (1985). Fitting variogram models by weighted least squares. *Journal of the International Association for Mathematical Geology, 17*(5), 563–586. https://doi.org/10.1007/BF01032109

De Andrade, L. C., Silva, A. A. E., Veloso, G. V., Filho, E. I. F., & Ferreira, I. O. (2025). Comparison of deterministic, probabilistic and machine learning-based methods for bathymetric surface modeling. *Modeling Earth Systems and Environment, 11*(1), 6. https://doi.org/10.1007/s40808-024-02189-8

Emery, B. E. (2024). Fair and Reasonable: A Conceptual Insight into USACE Dredge Estimating. *Journal of Waterway, Port, Coastal, and Ocean Engineering, 150*(5), 05024001. https://doi.org/10.1061/JWPED5.WWENG-2104

Henrico, I. (2021). Optimal interpolation method to predict the bathymetry of Saldanha Bay. *Transactions in GIS.* https://doi.org/10.1111/tgis.12759

Lenth, R. V. (2024). emmeans: Estimated marginal means, aka least-squares means (Version 1.10.0) [Computer software]. https://CRAN.R-project.org/package=emmeans

Li, Z., Peng, Z., Zhang, Z., Chu, Y., Xu, C., Yao, S., … Ma, J. (2023). Exploring modern bathymetry: A comprehensive review of data acquisition devices, model accuracy, and interpolation techniques for enhanced underwater mapping. *Frontiers in Marine Science, 10*, 1178845. https://doi.org/10.3389/fmars.2023.1178845

Müller, S., Schüler, L., Zech, A., & Heße, F. (2022). GSTools v1.3: A toolbox for geostatistical modelling in Python. *Geoscientific Model Development, 15*(7), 3161–3182. https://doi.org/10.5194/gmd-15-3161-2022

Murphy, B. S., Band, L. E., & Band, R. B. (2023). PyKrige (Version 1.7.1) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.7887527

Pratomo, D. G., Safira, R. A. D., & Stefani, O. (2023). A comparison of different GIS-based interpolation methods for bathymetric data: Case study of Bawean Island, East Java. *Geodesy and Cartography, 49*(4), 186–194. https://doi.org/10.3846/gac.2023.18250

Šiljeg, A., Lozić, S., & Šiljeg, S. (2015). A comparison of interpolation methods on the basis of data obtained from a bathymetric survey of Lake Vrana, Croatia. *Hydrology and Earth System Sciences, 19*(8), 3653–3666. https://doi.org/10.5194/hess-19-3653-2015

Uieda, L., Soler, S. R., Rampin, R., van Kemenade, H., Turk, M., Shapero, D., Banihirwe, A., & Leeman, J. (2023). Verde (Version 1.8.1) [Computer software]. https://doi.org/10.5281/zenodo.593011

USACE. (2013). *Hydrographic surveying (EM 1110-2-1003).* https://www.publications.usace.army.mil/Portals/76/Publications/EngineerManuals/EM_1110-2-1003.pdf

USACE. (n.d.). eHydro survey data [Dataset]. ArcGIS REST Services Directory. Retrieved January 31, 2026, from https://services7.arcgis.com/n1YM8pTrFmm7L4hs/ArcGIS/rest/services/eHydro_Survey_Data/FeatureServer/

Virtanen, P., Gommers, R., Oliphant, T. E., … SciPy 1.0 Contributors. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nature Methods, 17*(3), 261–272. https://doi.org/10.1038/s41592-019-0686-2
