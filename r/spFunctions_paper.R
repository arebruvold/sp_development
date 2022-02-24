# README ####
# Provided under a Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) license.
# For now the readme is incomplete, plans for major revision to improve performance and user friendliness.
# Must setup export of raw files using rev3 agilent script such that the format 
# of the .csv files are correct. Output in counts.
# Multiple isotopes of same element is as of now not supported
# Ionic standards format:
#  - 100Ti50Si10Au is read as a mixture of 100 ppb Ti, 50 ppb Si, 10 ppb Au.
#  - Minium two points per calibration, can be 0. (E.g. 0TiSiAu)
#  - NP reference material 
# RM datafile must be specified, only using single RM advised and supported.
# Sample names and datafile names must be unique.
# If density file for isotope is empty, n particles will be empty and bugs may occur.

# Packages ####

library(zoo) #custom rolling windows
library(tidyverse) #backbone
library(furrr) #parallel computing
library(beepr) #alarm when finishing, not essential
library(tictoc) #timer, not essential
plan(multisession, workers = 14)

# Speciation configuration  ####

# Densities, can be specified separately.
# From CRC Handbook of Chemistry and Physics, 97th edition where available
dens_comps <- tibble(
  "isotope" = c("RM", "Au", "Al", "Mn", "Pb", "Fe", "Si", "Ti", "Cr", "Ce", "Zr", "Cu", "Cd", "Ba", "Co", "Ni", "Zn"),
  "density" = c(19.32, 19.32, 2.56, 4.25, 6.29, 4.30, 2.65, 4.17, 5.22, 7.22, 5.68, 6.31, 2, 2, 2, 2, 2),
  "element_fraction" = c(1, 1, 26.982/(39.098+26.982+3*28.085+8*15.999), 1 / 1.6, 1 / 1.465, 1 / 1.59, 1 / 2.139, 1/1.67, 1/1.462, 1/1.23, 91.22/(91.22+16.00*2), 63.55/(63.55+16), 1, 1, 1, 1, 1)
)

mass_fractioner <- function(chemformula, element) {
  # Converts formula into mass fraction for given element.
  # Number after element indicates number of times this element occurs.
  # Case sensitive, number after element, e.g. TiO2, or C6H12O6.
  # Credits to https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee#file-periodic-table-of-elements-csv and
  # https://gist.github.com/robertwb/22aa4dbfb6bcecd94f2176caa912b952 for the periodic table.
  composition <- tibble(symbol = str_split(chemformula, "(?<=[A-Za-z\\d]{1,3})(?=[A-Z])") %>%
    as_vector()) %>%
    mutate(
      n = str_extract(symbol, "[\\d]{1,2}") %>% as.numeric(),
      n = replace_na(n, 1),
      symbol = str_replace_all(symbol, "\\d", "")
    ) %>%
    inner_join(read_csv("periodic_table.csv") %>% janitor::clean_names() %>%
      select(symbol, atomic_mass), by = "symbol") %>%
    mutate(n_x_atomic_mass = atomic_mass * n)

  composition %>%
    filter(symbol == element) %>%
    pull(n_x_atomic_mass) / sum(composition %>% pull(n_x_atomic_mass))
  
}

# Supporting functions ####

# Maximum of density fn estimator - to get the peak of non-parametric, multimodal distribution.
# Could potentially be replaced with rccp equiv for performance gains.
dens_max <- function(INDENS){
  dens_test <- density(INDENS, from = 0, width = 2.5)
  dens_test$x[which.max(dens_test$y)]
}

# Filter for low baselines, in which case density estimates can be imprecise for integer count data.
# Could be improved on, e.g. by use of asymmetric kernel, or by defining the kernel.
low_c_mean <- function(INCOUNT){
 # INCOUNT %>% pmin(2) %>% mean()
  INCOUNT[INCOUNT>7] <- 7 #should be quicker than code above
  INCOUNT %>% mean(trim = 0)
}

# required to calculate threshold from non-integer values
# Credits to stats.stackexchange.com/questions/10926/how-to-calculate-confidence-interval-for-count-data-in-r
PoiCI <- function (num, conf.level = 0.95) {
  a = 1 - conf.level
# lower <- 0.5 * qchisq(a/2, 2*num)
  0.5 * qchisq(1-a/2, 2*num+2)
}

#get list of all count files in directory
folder_rawfinder <- function(folder) {
  # if a .csv is the input, returns a single .csv
  if (folder %>% str_detect("\\.csv")) {
    folder
  } else {
    list.files(
      path = paste(folder), recursive = FALSE,
      pattern = "\\.csv$",
      full.names = TRUE
    ) %>%
      # incase / is added after path.
      str_replace("//", "/") %>%
      as_tibble() %>%
      # only counts / cps files
      filter(str_detect(value, "\\_count|\\_cps")) %>%
      # removes non-unique raw files (generated from Agilent script bug). This could perhaps be more robust.
      # Is no longer needed with new script.
      # distinct(str_remove(value, "(?<=\\/[0-9]{3})[^\\/]*(?=\\_\\d\\.csv$)"), .keep_all = TRUE) %>%
      pull(value) %>%
      as_vector()
  }
}

folder_rawfinder_PE <- function(folder) {
  # if a .csv is the input, returns a single .csv
  # Perkin Elmer files.
  if (folder %>% str_detect("\\.csv")) {
    folder
  } else {
    list.files(
      path = paste(folder), recursive = FALSE,
      pattern = "\\.csv$",
      full.names = TRUE
    ) %>%
      # incase / is added after path.
      str_replace("//", "/") %>%
      as_tibble() %>%
      # only counts / cps files
      # filter(str_detect(value, "\\_count|\\_cps")) %>%
      pull(value) %>%
      as_vector()
  }
}

# Poisson cps generator 
signal_simulator <-
  function(lambda, dwells, isotope = NA_character_) {
    as_tibble(seq(0, dwells - 1, by = 1) / 10000) %>%
      mutate(counts = rpois(dwells, lambda)) %>%
      mutate(
        isotope = isotope,
        sample_name = "Simulated",
        datafile = "Simulated"
      ) %>%
      rename(c(
        "Time" = "value",
        "counts" = "counts",
        "isotope" = "isotope"
      ))
  }

# Simulate from sample baseline
sample_simulator <-
  function(foundNPs) {
    sim_po <- foundNPs %>%
      mutate(
        counts = rpois(n = n(), lambda = baseline),
        sample_name = paste(sample_name, "_sim_po", sep = "")
      )
    sim_nb <- foundNPs %>%
      mutate(
        counts = rnbinom(n = n(), mu = baseline, size = 4.9498),
        sample_name = paste(sample_name, "_sim_nb", sep = "")
      )
    bind_rows(sim_po, sim_nb, foundNPs)
  }


# Single particle analysis ####

# raw count file reader
raw_reshaper <- function(CPS_csv, prop = prop) {
  header_info <- read_lines(CPS_csv, skip = 0, n_max = 2)
  reshaped <- read_csv(
    CPS_csv,
    skip = 4,
    col_types = "dd"
  ) %>%
    rename_with(~ str_replace(.x, ".*197.*", "Au")) %>%
    rename_with(~ str_extract(.x, "[a-zA-Z]+")) %>%
    drop_na() %>%
    slice_head(prop = prop) %>%
    mutate(
      sample_name = header_info[2],
      dataset = header_info[1] %>% 
        str_replace_all("\\\\", "/") %>%
        str_extract("(?<=\\/)[^\\/]*(?=\\/[^\\/]*\\.d$)"),
      datafile = header_info[1] %>%
        str_replace_all("\\\\", "/") %>%
        str_extract("(?<=\\/)\\d{3}[^\\/]*(?=\\.[dD]$)"),
      isotope = colnames(.)[2]
    ) %>%
    rename(counts = 2) %>%
    mutate(counts = round(counts, digits = 0))
}

# raw count file reader for Perkin Elmer
# Change naming to use file name if used for this.
raw_reshaper_PE <- function(CPS_csv, prop = 1) {
  read_csv(CPS_csv,
           skip = 0,
           col_types = "dd") %>%
    select(1) %>%
    drop_na() %>%
    mutate(isotope = colnames(.)[1]) %>% 
    rename(counts = 1) %>%
    slice_head(prop = prop) %>%
    mutate(
      sample_name = CPS_csv %>% str_replace_all("\\\\", "/") %>% 
        str_replace_all(c(" " =  "",
                          "\\\\" = "/",
                        "[\\(\\)]" = "_")) %>% 
        str_extract("(?<=\\/)[^\\/]+(?=\\.csv$)") %>% 
        str_replace("\\.", ""),
      dataset = CPS_csv %>% str_replace_all("\\\\", "/") %>%
        str_extract("(?<=\\/)[^\\/]*(?=\\/[^\\/]+\\.csv$)"),
      datafile = sample_name,
      counts = round(counts, digits = 0)
    ) %>%
    mutate(Time = row_number() / 10000, .before = 1)
}

# baseline and threshold estimation
# To self: Poisson mode: https://www.immagic.com/eLibrary/ARCHIVES/GENERAL/WIKIPEDI/W121109P.pdf
baseline_thr_finder <- function(reshaped) {
  bl_thr <- reshaped %>%
    # robust baseline with kde
    mutate(baseline = rollapply(counts, 301, dens_max, by = 50, fill = NA)) %>%
    # end and and start fill inn due to windowed operation
    mutate(baseline = na.fill(baseline, "extend")) %>%
    # Zero inflation, problem of integers and bounded distribution
    # for kernel density vs bounded distribution is corrected:
    # (Could be removed after kernel density width change?)
     mutate( 
      baseline =
        if_else(
          baseline <= 3,
          rollapply(counts, 301, low_c_mean, fill = NA, by = 50),
          baseline
        )
    ) %>%
    mutate(
      baseline = na.fill(baseline, "extend")
    ) %>%
    mutate(
      mean_baseline = mean(baseline, trim = 0.1),
      # to fix detector dropping making baseline too low
      baseline = if_else(baseline < 0.75 * mean_baseline |
        # and in case of aggregates/ high intensity vs zero/low background
        baseline > 3 * mean_baseline,
        mean_baseline,
        baseline
      )
    ) %>%

    # Poisson threshold for filtering particle events, minimum of 5 counts.
    # (Could consider using qpois + baseline - mean_baseline, yet may
    # not be robust if inhomogenous sample/ surface adsorption.
    # Use mean_baseline for sample comparison, yet never setting 
    # below preexisting baseline)
    mutate(h_thr = rollmax(baseline, 400, by = 50, fill = NA),
           h_thr = na.fill(h_thr, "extend"),
      h_thr = qpois(
        p = (1 - 0.05 / 600000),
        lambda = if_else(h_thr > 0.2, h_thr, 0.2))
    )  %>%   # drop redundant column(s?)
  select(!mean_baseline)
}

# integration of peaks, shape parameters
peak_integration <- function(bl_thr) {
  peaks <- bl_thr %>% 
    mutate(
      c_above_baseline_thr = if_else(
        counts <= baseline,
        0,
        counts - baseline
      ),
      # generates a 1 for each time a value is above the baseline, and the previous point is below.
      peak_n = if_else(
        c_above_baseline_thr != 0 &
          (lag(c_above_baseline_thr) == 0 |
            is.na(lag(
              c_above_baseline_thr
            ))),
        1,
        0
      )
    ) %>%
    # labeling by n event, classify events whether max is above thr or not.
    mutate(
      peak_n = cumsum(peak_n)
    ) %>%
    group_by(peak_n) %>%
    #Peak width error
    mutate(
      peak_max_I = max(counts),
      #peak_width = n(), wrong, counts until next signal..
      peak_width = sum(c_above_baseline_thr > 0),
      above_height_thr = peak_max_I > h_thr 
    ) %>%
    # integration using group by.
    # interpolation when passing baseline could perhaps improve accuracy
    # across size ranges at the expense of computation time.
    # This because the aspect ratio may increase w size. 
    group_by(peak_n) %>%
    mutate(peak_area = sum(c_above_baseline_thr) %>% round(digits = 0)) %>%
    # for signal distribution:
    group_by(peak_area, above_height_thr) %>%
    mutate(
      peak_area_count = length(unique(peak_n))
    ) %>%
    ungroup() %>% 
    select(!c_above_baseline_thr)
  
  peaks
}

# wrapper function
nParticle_finder_blocks <- function(CPS_csv, prop = 1) {
  CPS_csv %>%
    folder_rawfinder() %>%
    future_map_dfr(
      ~ .x %>%
        # read in raw (cps) file
        raw_reshaper(prop = prop) %>%
        # calculate baseline
        baseline_thr_finder() %>%
        # integrate peaks
        peak_integration(),
      .options = furrr_options(seed = TRUE)
    )
}

nParticle_finder_blocks_PE <- function(CPS_csv, prop = 1) {
  CPS_csv %>%
    folder_rawfinder_PE() %>%
    future_map_dfr(
      ~ .x %>%
        # read in raw (cps) file
        raw_reshaper_PE(prop = prop) %>%
        # calculate baseline
        baseline_thr_finder() %>%
        # integrate peaks
        peak_integration(),
      .options = furrr_options(seed = TRUE)
    )
}

# # timing blocks fn
# nParticle_finder_timer <- function(CPS_csv, prop) {
#  { tic("reshaping")
#  reshaped <- CPS_csv %>%
#         # read in raw (cps) file
#         raw_reshaper(prop = prop)
#  toc()}
#   
# { tic("baselines")
#  baselines_thrs <- reshaped %>% 
#         # calculate baseline
#         baseline_thr_finder()
#  toc()}
# {tic("integrations:")
#  peaks_int <- baselines_thrs %>% 
#         # integrate peaks
#         peak_integration()
#  toc()}
# 
# }

# Calibrations & output####
# NOT OPTIMIZED

# getting response factors, standard recognition
signal_response_factor <- function(foundNPs) {

  # Recognize files containing STDs.
  isotope_list <- foundNPs %>% distinct(isotope)

  STDs <- foundNPs %>%
    filter(
      sample_name %>% str_detect(as_vector(isotope_list)),
      sample_name %>% str_detect("\\d{1,3}[A-Z]{1}[a-z]{1}"), #+ not added!
      sample_name %>% str_detect("[A-Z]{1}[a-z]{1}")
    ) %>%
    # calculate mean counts for entire timescan (per dwell)
    group_by(datafile, isotope) %>%
    summarise(
      sample_name = unique(sample_name),
      mean_counts = mean(counts)
    ) %>%
    filter(sample_name %>% str_detect(isotope)) %>%
    # add concentration of each added element
    mutate(
      STD_conc =
        str_extract(
          sample_name,
          paste0("[\\d.]+(?=[\\sa-zA-Z]*", isotope, ")")
        ) %>% as.numeric()
    ) %>% ungroup() %>% 
    # linear model for each isotope to get counts per dwell per ppb aka ug per liter
    group_by(isotope) %>%
    mutate(
      coef_a =
        coef(
          lm(formula = as.numeric(mean_counts) ~ STD_conc)
        )[1],
      coef_b =
        coef(
          lm(formula = as.numeric(mean_counts) ~ STD_conc)
        )[2],
      r2 = summary(
        lm(formula = as.numeric(mean_counts) ~ STD_conc)
      )$r.squared
    ) %>%
    ungroup()
}


# Particle parameters and mass calibration using RM
signal_conc_calibration <- function(foundNPs,
                                    dens_comps,
                                    RM_datafile,
                                    RM_dia = 60,
                                    RM_isotope = "Au",
                                    supp_coef = 1) {
  # First translates to element mass per count (meaning NOT particle mass)
  # This is done by using Au RM to relate its mass to the counts.
  # Then for other isotopes, the ratio between the intensities (b_coef) is used to convert to element mass per count for these...
  #supp_coef: matrix suppression coefficient, 40% reduction in sensitivity corresponds to a factor of 0.6.
  
  # Create table of RM peak areas:
  RM_areas <- foundNPs %>%
    filter(
      above_height_thr == TRUE,
      #sample_name %>% str_detect("RM"),
      datafile %>% str_detect("RM|PerkinElmer"),
      datafile %>% str_detect(RM_datafile),
      isotope == RM_isotope
    ) %>%
    group_by(peak_n) %>%
    summarise(RM_area = mean(peak_area))
  
  # Use density function to retrieve the "mode" of the peak are of the RM. [REVISIT: How does density fn work?] OPTIMIZE
  counts_RM <-
    density(RM_areas$RM_area)$x[which.max(density(RM_areas$RM_area)$y)]
  
  # Get the mass per count by relating the "known" RM mass to counts_RM
  mass_count_kg <- # counts per kg RM element
    (
      dens_comps %>% filter(isotope == "RM") %>% pull(density) * # density RM
        1000 * (4 / 3) * pi * (RM_dia * 10 ^ (-9) / 2) ^ 3 * # RM volume
        (
          dens_comps %>% filter(isotope == "RM") %>% pull(element_fraction)
        )
    ) /
    counts_RM # element fraction divided by RM peak area "mode"
  
  # Response factors
  RFs <- signal_response_factor(foundNPs) %>%
    select(isotope, coef_b) %>%
    distinct(coef_b, .keep_all = TRUE)
  
  # RM response factor
  RF_RM <- RFs %>%
    filter(isotope == RM_isotope) %>%
    distinct(coef_b) %>%
    as.numeric()
  
  # Detector mass flow rate (L/dwell) (equiv to Transport efficiency * sample intake flow rate)
  detector_flow_rate <-
    mass_count_kg * RF_RM * 10 ^ 9 #* 10^9 to account for ug/L to kg/L
  
  # Use signal_response factors to correct the mass per count for each isotope
  np_events <- foundNPs %>%
    group_by(isotope, datafile) %>%
    distinct(peak_n, .keep_all = TRUE) %>%
    # filter(peak_max_I/peak_width >) # Discriminate by aspect ratio?
    inner_join(dens_comps) %>% inner_join(RFs) %>%
    ungroup() %>%
    group_by(isotope) %>%
    # RF_RM / coef_b normalization/ conversion to different isotope
    mutate(peak_mass = mass_count_kg * peak_area * supp_coef * RF_RM  / coef_b) %>%
    group_by(isotope, datafile) %>% 
    mutate(
      ionic_conc = mean(baseline) * supp_coef / coef_b) %>% 
    ungroup() %>%
    
    # Translate to particle size from peak_area using appropriate equation and dens_comp_particles. 
    # Thus PARTICLE mass, including elements not measured on ICP.
    mutate(size_nm = (6 * peak_mass /
                        (pi * density * 1000 * element_fraction)) ^ (1 / 3) * 10 ^ 9) %>%
    # Volumetric flow rate to detector, (liter per 100 us) and
    # Transport efficiency, assuming flow rate of 0.346 mL/min
    mutate(
      TE = detector_flow_rate*60/(100*0.346*10^(-9)),
      detector_flow_rate = detector_flow_rate
    )
  
  # sample pass through. NEED TO BE VALIDATED.
  all_samples <- foundNPs %>% distinct(sample_name, isotope)
  
  np_events_new <- left_join(all_samples, np_events, by = c("sample_name", "isotope"))
  # mutate(TE = na.fill(TE, "extend")) %>%
  # replace(is.na(.), 0)
  #
}

# Particle parameters and mass calibration using RM
signal_conc_calibration_PE <- function(foundNPs,
                                    dens_comps,
                                    RM_datafile,
                                    RM_dia = 60,
                                    RM_isotope = "Au") {
  # First translates to element mass per count (meaning NOT particle mass)
  # This is done by using Au RM to relate its mass to the counts.
  # Then for other isotopes, the ratio between the intensities (b_coef) is used to convert to element mass per count for these...
  
  
  # Create table of RM peak areas:
  RM_areas <- foundNPs %>%
    filter(
      above_height_thr == TRUE,
      datafile %>% str_detect("RM"),
      isotope == RM_isotope
    ) %>%
    group_by(peak_n) %>%
    summarise(RM_area = mean(peak_area))
  
  # Use density function to retrieve the "mode" of the peak are of the RM. [REVISIT: How does density fn work?] OPTIMIZE
  counts_RM <-
    density(RM_areas$RM_area)$x[which.max(density(RM_areas$RM_area)$y)]
  
  # Get the mass per count by relating the "known" RM mass to counts_RM
  mass_count_kg <- # counts per kg RM element
    (
      dens_comps %>% filter(isotope == "RM") %>% pull(density) * # density RM
        1000 * (4 / 3) * pi * (RM_dia * 10 ^ (-9) / 2) ^ 3 * # RM volume
        (
          dens_comps %>% filter(isotope == "RM") %>% pull(element_fraction)
        )
    ) /
    counts_RM # element fraction divided by RM peak area "mode"
  
  # Response factors
  RFs <- signal_response_factor(foundNPs) %>%
    select(isotope, coef_b) %>%
    distinct(coef_b, .keep_all = TRUE)
  
  # RM response factor
  RF_RM <- RFs %>%
    filter(isotope == RM_isotope) %>%
    distinct(coef_b) %>%
    as.numeric()
  
  # Detector mass flow rate (L/dwell) (equiv to Transport efficiency * sample intake flow rate)
  # MAY HAVE BEEN ALTERED, CHECK!
  detector_flow_rate <-
    mass_count_kg * RF_RM * 10 ^ 9 #* 10^9 to account for ug/L to kg/L
  
  # Use signal_response factors to correct the mass per count for each isotope
  np_events <- foundNPs %>%
    group_by(isotope, datafile) %>%
    distinct(peak_n, .keep_all = TRUE) %>%
    # filter(peak_max_I/peak_width >) # Discriminate by aspect ratio?
    inner_join(dens_comps) %>% inner_join(RFs) %>%
    ungroup() %>%
    group_by(isotope) %>%
    # Debug sizes: SWAPPED FROM coef_B / RF_RM 19.08.21
    mutate(peak_mass = mass_count_kg * peak_area * RF_RM / coef_b) %>%
    ungroup() %>%
    
    # Translate to particle size from peak_area using appropriate equation and dens_comp_particles. 
    # Thus PARTICLE mass, including elements not measured on ICP.
    mutate(size_nm = (6 * peak_mass /
                        (pi * density * 1000 * element_fraction)) ^ (1 / 3) * 10 ^ 9) %>%
    # Volumetric flow rate to detector, (liter per 100 us) and
    # Transport efficiency, assuming flow rate of 0.346 mL/min
    mutate(
      TE = detector_flow_rate*60/(100*0.346*10^(-9)),
      detector_flow_rate = detector_flow_rate
    )
  
  # sample pass through. NEED TO BE VALIDATED.
  all_samples <- foundNPs %>% distinct(sample_name, isotope)
  
  np_events_new <- left_join(all_samples, np_events, by = c("sample_name", "isotope"))
  # mutate(TE = na.fill(TE, "extend")) %>%
  # replace(is.na(.), 0)
  #
}

# Generating summary table of samples
nParticle_summariser <- function(foundNPs_sizes, acq_time) {
  above_T <-
    foundNPs_sizes %>%
    filter(above_height_thr == TRUE) %>%
    group_by(isotope, sample_name) %>%
    distinct(peak_n, .keep_all = TRUE) %>%
    summarise(
      TE = mean(TE),
      n_particles = n(),
      mass_conc = sum(peak_mass) * 10^12 /
        (mean(detector_flow_rate) * acq_time * 10000),
      particle_conc = n_particles * 10^0 /
        (mean(detector_flow_rate) * acq_time * 10000),
      density = mean(density),
      mean_size = mean(size_nm),
      median_size = median(size_nm),
      mean_h_thr = mean(h_thr),
      mean_baseline = mean(baseline),
      ionic_conc = mean(ionic_conc)
    )
  all_samples <- foundNPs_sizes %>% distinct(sample_name, isotope)

  left_join(all_samples, above_T, by = c("sample_name", "isotope")) %>%
    mutate(TE = na.fill(TE, "extend")) %>%
    replace(is.na(.), 0) # %>%
  # mutate(across(everything(), ~ .x %>%
  #                 format(
  #                   scientific = TRUE, digits = 2
  #                 )))
}


peakshaper <- function(sample_isotope, peak_width = 20) {
  grp_size <- peak_width * 2 + 1

  sample_isotope_tops <- sample_isotope %>%
    ungroup() %>%
    mutate(peak_top = counts == peak_max_I & above_height_thr == TRUE)


  peaks_and_around <-
    outer(which(sample_isotope_tops$peak_top == TRUE), -peak_width:peak_width, `+`) %>%
    as_tibble(.name_repair = "unique") %>%
    pivot_longer(cols = everything()) %>%
    select(-name) %>%
    filter(value >= 1 &
      value <= nrow(sample_isotope_tops)) %>%
    pull(value) %>%
    as_vector()


  sample_isotope[peaks_and_around, ] %>%
    mutate(peak_grp = as.integer(gl(n(), grp_size, n()))) %>%
    group_by(peak_grp) %>%
    mutate(rel_pos = (Time - median(Time)) * 10000)
}

peakshaper_wrapper <- function(foundNPs, peak_width = 20) {
  foundNPs %>%
    group_by(isotope, sample_name) %>%
    group_modify(~ peakshaper(.x, peak_width)) %>%
    ungroup()
}

peakshape_summary <- function(foundNPs) {
  foundNPs %>%
    group_by(sample_name, isotope) %>%
    filter(above_height_thr == TRUE) %>%
    distinct(peak_n, .keep_all = TRUE) %>%
    summarise(
      mean_peak_area = mean(peak_area),
      mean_peak_width = mean(peak_width),
      mean_peak_height = mean(peak_max_I),
      sd_peak_area = sd(peak_area),
      sd_peak_width = sd(peak_width),
      sd_peak_height = sd(peak_max_I),
      n_peaks = n()
    )
}

# sP plotters ####

# Extract b coefficients by distinct, use these for calibration
# NOT OPTIMIZED

peakshape_plotter <- function(foundNPs) {
  foundNPs %>%
    select(sample_name, counts, rel_pos, peak_grp) %>%
    ungroup() %>%
    mutate(rel_pos = round(rel_pos, digits = 0)) %>%
    dplyr::group_by(sample_name, rel_pos) %>%
    dplyr::summarise(
      mean_counts = mean(counts),
      sd_counts = sd(counts)
    ) %>%
    ggplot(aes(x = rel_pos, y = mean_counts)) +
    geom_line(alpha = 0.3, size = 1) +
    geom_ribbon(aes(
      ymin = mean_counts - sd_counts,
      ymax = mean_counts + sd_counts, fill = sample_name
    ), alpha = 0.2) +
    facet_wrap(~sample_name)
}


signal_rf_plotter <- function(signal_response_factors) {
  signal_response_factors %>%
    ggplot(aes(STD_conc, mean_counts)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x) +
    labs(
      x = paste0("Concentration [ug/L]"),
      y = "Intensity [counts/dwell]"
    ) +
    geom_text(
      x = -Inf, y = Inf, hjust = 0, vjust = 1,
      aes(label = paste("slope:", round(coef_b, digits = 3)))
    ) +
    geom_text(
      x = -Inf, y = Inf, hjust = 0, vjust = 2.5,
      aes(label = paste("intercept:", round(coef_a, digits = 3)))
    ) +
    geom_text(
      x = -Inf, y = Inf, hjust = 0, vjust = 4,
      aes(label = paste("r2:", round(r2, digits = 5)))
    ) +
    facet_wrap(~isotope, scales = "free")
}

# NOTE HERE IS COUNTED PEAK AREAS FROM FULL DATASET IF YOU SLICE!!!

signal_dist_plotter <- function(foundParticles, peak_area_thr = 0) {
  sdist_labels <- foundParticles %>%
    group_by(isotope, sample_name) %>%
    distinct(peak_n, .keep_all = TRUE) %>%
    summarise(
      mean_baseline = mean(baseline),
      nparticles = sum(above_height_thr == TRUE)
    )

  foundParticles %>%
    distinct(isotope, sample_name, peak_n, .keep_all = TRUE) %>%
    filter(peak_area >=
      peak_area_thr & peak_area > 0) %>%
    ggplot(aes(x = peak_area, fill = above_height_thr, alpha = above_height_thr)) +
    geom_bar(width = 0.005, position = position_dodge()) +
    scale_y_continuous(trans = "pseudo_log") +
    scale_x_log10() +
    labs(title = "Intensity-discriminated signal distributions") +
    scale_alpha_discrete(range = c(0.3, 1)) +
    geom_text(
      data = sdist_labels,
      x = -Inf,
      y = Inf,
      hjust = 0,
      vjust = 1,
      inherit.aes = FALSE,
      aes(label = paste("particles:", nparticles))
    ) +
    geom_text(
      data = sdist_labels,
      x = -Inf,
      y = Inf,
      hjust = 0,
      vjust = 2.5,
      inherit.aes = FALSE,
      aes(label = paste("BL:", mean_baseline %>% round(2)))
    ) +
    # Labeller? https://newbedev.com/how-to-change-facet-labels
    facet_grid(isotope ~ sample_name)
}


## In the case of no particles above height thr, will stop, this due to empty object sent to aesthetics
signal_found_plotter <- function(foundParticles,
                                 prop = 1,
                                 t_0 = 0, t_e = 600,
                                 i_min = 0,
                                 i_max = Inf,
                                 area_below = Inf) {
  if (ncol(foundParticles) <= 16) {
    plot_Found <- foundParticles %>%
      group_by(isotope, datafile) %>%
      slice_head(prop = prop) %>%
      ungroup() %>%
      filter(
        Time > t_0,
        Time < t_e,
        peak_max_I >= i_min,
        peak_max_I <= i_max
      ) %>%
      ggplot(aes(x = Time, y = counts)) +
      # geom_line(size = 0.2) +
      theme_minimal() +
      geom_line() +
      geom_line(aes(x = Time, y = baseline), color = "blue", size = 2) +
      geom_line(aes(x = Time, y = h_thr), color = "red") +
      geom_point(
        data = foundParticles %>% group_by(isotope, datafile) %>%
          slice_head(prop = prop) %>% distinct(peak_n, .keep_all = TRUE) %>%
          mutate(Time_event = if_else(above_height_thr == TRUE &
            peak_area < area_below, Time, NA_real_)),
        aes(
          x = Time_event,
          y = 0
        ),
        #   label = round(Time_event, digits = 2),
        #   angle = 90,
        #   hjust = 1
        # ),
        size = 3,
        color = "green"
      ) +
      # geom_text(
      #   data = foundParticles %>% group_by(isotope, datafile) %>%
      #     # Use FILTER instead of the mutate
      #     slice_head(prop = prop) %>% distinct(peak_n, .keep_all = TRUE) %>%
      #     filter(above_height_thr == TRUE,
      #            peak_area < area_below),
      #   aes(
      #     x = Time,
      #     y = 0,
      #     label = round(peak_area, digits = 0),
      #     angle = 90
      #   ),
      #   size = 3,
      #   color = "green"
      # ) +
      labs(title = "Signal processing") +
      # # labels/text, OLD, NEED HELP WITH THIS:
      # labs(subtitle = paste(
      #   foundParticles %>% filter(above_height_thr == "Yes") %>%
      #     distinct(peak_n) %>% nrow(),
      #   " particle events",
      #   ", min intensity threshold: ", min(foundParticles$h_thr),
      #   "\n", foundParticles$isotope[1], ":  ",
      #   foundParticles$datafile, ":  ",
      #   foundParticles$dataset,
      #   sep = ""
      # )) +
      facet_wrap(isotope ~ sample_name, scales = "free")
  }

  if (ncol(foundParticles) > 16) {
    plot_Found <- foundParticles %>%
      group_by(isotope, datafile) %>%
      slice_head(prop = prop) %>%
      ungroup() %>%
      filter(
        Time > t_0,
        Time < t_e,
        peak_max_I >= i_min,
        peak_max_I <= i_max
      ) %>%
      ggplot(aes(x = Time, y = counts)) +
      # geom_line(size = 0.2) +
      theme_minimal() +
      geom_line() +
      geom_line(aes(x = Time, y = baseline), color = "blue", size = 2) +
      geom_line(aes(x = Time, y = h_thr), color = "red") +
      # geom_point(
      #   data = foundParticles %>% group_by(isotope, datafile) %>%
      #     slice_head(prop = prop) %>%  distinct(peak_n, .keep_all = TRUE) %>%
      #     mutate(Time_event = if_else(above_height_thr == "Yes", Time, NA_real_)),
      #   aes(
      #     x = Time_event,
      #     y = 2)
      #   ,
      #   #   label = round(Time_event, digits = 2),
      #   #   angle = 90,
      #   #   hjust = 1
      #   # ),
      #   size = 2,
      #   color = "green") +
      geom_text(
        data = foundParticles %>% group_by(isotope, datafile) %>%
          # Use FILTER instead of the mutate
          slice_head(prop = prop) %>% distinct(peak_n, .keep_all = TRUE) %>%
          filter(above_height_thr == TRUE),
        aes(
          x = Time,
          y = -5,
          label = round(size_nm, digits = 0),
          angle = 90
        ),
        size = 3,
        color = "green"
      ) +
      labs(title = "Signal processing") +
      # # labels/text, OLD, NEED HELP WITH THIS:
      # labs(subtitle = paste(
      #   foundParticles %>% filter(above_height_thr == "Yes") %>%
      #     distinct(peak_n) %>% nrow(),
      #   " particle events",
      #   ", min intensity threshold: ", min(foundParticles$h_thr),
      #   "\n", foundParticles$isotope[1], ":  ",
      #   foundParticles$datafile, ":  ",
      #   foundParticles$dataset,
      #   sep = ""
      # )) +
      facet_wrap(isotope ~ sample_name, scales = "free")
  }
  plot_Found
}


