# Analysis code for "Tracing echoes from the past: a district-level analysis of thyroid cancer incidence and radiation exposure in the post-Chornobyl Ukraine"

Nataliia Levchuk, Jonas Schöley, Laura Ann Cilek, Domantas Jasilionis

*For the possibility of accessing the input data please contact Nataliia Levchuk.*

-   `src`
    -   `00-global_functions.R`: Functions imported by other scripts.
    -   `10-map_templates.R`: Generate geodata and map templates for Ukrainian regional analyses.
    -   `11-model_input.R`: Prepare Ukrainian data on regional thyroid cancer incidence and radiation exposure for modeling.
    -   `20-model_crudesir.R`: Calculate crude regional SIR of thyroid cancer incidence across Ukraine.
    -   `21-model_smoothsir.R`: Calculate smooth regional SIR of thyroid cancer incidence across Ukraine.
    -   `22-model_dosage.R`: Estimate exposure rate ratio of absorbed dosage.
    -   `install_dependencies.R`: Install required packages.
    -   `NL`: Preparation of dosage and population data and additional analyses by Nataliia Levchuk and Laura Ann Cilek.
        -   `Part6_NL.R`: Spatial interpolation of settlement dosage data to get district mean doses and plus interpolation of population weighted doses at the age 15-19.
        -   `Part7_NL.R`: Measuring spatial autocorrelation, Moran I, and Lisa plotting.
        -   `Part9_NL.R`: Calculating district average thyroid doses in 1986 by six population percentiles (0, .10, .25, .50, .75, .90, 1).
        -   `Part10_NL.R`: Step bar plot for dosage data (1986 dosage15+ by six population percentiles).
        -   `Part11_NL.R`: Calculating unsmoothed district SIRs of thyroid cancer 15+ in 2001 by six population percentiles.
        -   `Part12_NL.R`: Step bar plot for unsmoothed SIRs at the age 15+ in 2001 by six population percentiles.
        -   `Part14_NL.R`: Calculating Ukraine’s population weights and standardized incidence rates at the age 15+.

![](README_files/teaser.png)
