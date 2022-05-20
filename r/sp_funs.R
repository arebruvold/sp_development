
{
  library(zoo) # custom rolling windows
  library(tidyverse) # backbone
  library(RcppRoll)
  library(furrr) # parallel computing
  library(assertr)
  plan(multisession, workers = 14)
  library(tictoc)
}


# supporting functions ####

# Maximum of density fn estimator - to get the peak of non-parametric, multimodal distribution.
# Could potentially be replaced with rccp equiv for performance gains.
dens_max <- function(INDENS){
  dens_test <- density(INDENS, from = 0, width = 2.5)
  dens_test$x[which.max(dens_test$y)]
}

# Filter for low baselines, in which case density estimates can be imprecise for integer count data.
# Could be improved on, e.g. by use of asymmetric kernel, or by defining the kernel.
low_c_mean <- function(INCOUNT){
  INCOUNT[INCOUNT>5] <- 5
  INCOUNT %>% mean(trim = 0.05)
}

# required to calculate threshold from non-integer values
# Credits to stats.stackexchange.com/questions/10926/how-to-calculate-confidence-interval-for-count-data-in-r
PoiCI <- function (num, conf.level = 0.95) {
  a = 1 - conf.level
# lower <- 0.5 * qchisq(a/2, 2*num)
  0.5 * qchisq(1-a/2, 2*num+2)
}

# sp functions ####

## file handling ####

# finds count files in given folder only, not recursive.
# - in: folder
# - out: n x 1 tibble with file paths, n = # count files.
sp_rawfinder <- function(folder){
  tryCatch(
    expr = {
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
          enframe(name = NULL, value = "filepath") %>% 
          # only counts / cps files
          filter(str_detect(filepath, "\\_count|\\_cps"))
        # pull(value)
      }
    },
    error = function(e){
      print(
        sprintf("An error occurred in sp_rawfinder at %s : %s",
                Sys.time(),
                e)
      )
    })
}

# raw count file reader
# - in: file path to .csv
# - out: m x 1 tibble with signal intensity count data, m = # dwells/ data points.
sp_reader <- function(CPS_csv) {
  tryCatch(
    expr = {
      reshaped <- read_csv(
        CPS_csv,
        skip = 4,
        col_select = 2,
        col_types = "d"
      ) %>% 
        rename(counts = 1) %>% drop_na() %>% 
        mutate(counts = round(counts, digits = 0))
      
      return(reshaped)
    },
    error = function(e) {
      print(
        sprintf(
          "An error occurred in sp_reader at %s : %s",
          Sys.time(),
          e
        )
      )
    }
  )
}

# classifies .csv files into SAMPLE, RM or STD, adds dataset, sample nae, datafile, isotope.
# - in: n x 1 tibble with file paths
# - out: tibble with metadata on the files.
# - todo: faster implementation using list.files into vector pipe into a map with a function that read lines, pipe into new map make function to take each lines and extract info with regex mutate.
sp_classifier <- function(folder, RM_string = "RM"){
  tryCatch(
    expr = {
      folder %>%
        sp_rawfinder() %>% 
        ungroup() %>%
        rowwise() %>% # consider rewrite map for performance
        mutate(
          sample_name = read_lines(filepath, skip = 0, n_max = 2)[2],
          dataset = read_lines(filepath, skip = 0, n_max = 2)[1] %>%
            str_replace_all("\\\\", "/") %>%
            str_extract("(?<=\\/)[^\\/]*(?=\\/[^\\/]*\\.d$)"),
          sample_name = read_lines(filepath, skip = 0, n_max = 2)[2],
          dataset = read_lines(filepath, skip = 0, n_max = 2)[1] %>%
            str_replace_all("\\\\", "/") %>%
            str_extract("(?<=\\/)[^\\/]*(?=\\/[^\\/]*\\.d$)"),
          datafile = read_lines(filepath, skip = 0, n_max = 2)[1] %>%
            str_replace_all("\\\\", "/") %>%
            str_extract("(?<=\\/)\\d{3}[^\\/]*(?=\\.[dD]$)"),
          isotope = read_csv(filepath,
                             skip = 4,
                             n_max = 0,
                             col_select = 2,
          ) %>% names() %>% str_replace(".*197.*", "Au") %>%
            str_extract("[A-Z]{1}[a-z]{1}"),
          type =
            case_when(
              str_detect(sample_name, RM_string) ~ "RM",
              str_detect(sample_name, "\\d{1,3}[A-Z]{1}[a-z]{1}") ~ "STD",
              TRUE ~ "SAMPLE"
            )
        ) %>%
        ungroup()
    },
    error = function(e){
      print(
        sprintf("An error occurred in sp_classifier at %s : %s",
                Sys.time(),
                e)
      )
    })
}



classified <- sp_classifier("~/sp-data/Fordefjorden/", RM_string = "AuRM500ng/LM")


## sp peak discrimination ####

# calculates baseline and poisson height threshold. Used within sp_peak_discriminator.
# - in: m x 1 tibble with count data
# - out: 3 x 1 tibble with count data, baseline and height threshold.
# - todo(?): reimplement check for fluctuating baseline.
sp_baseline_thr <- function(reshaped) {
  # median_counts <- median(reshaped$counts, na.rm = TRUE)
  bl_thr <- reshaped %>%
    drop_na() %>%
    # robust baseline with rolling median or kde
    mutate( # baseline = rollapply(counts, 301, dens_max, by = 50, fill = NA),
      baseline = roll_median(counts, n = 301, by = 50, fill = NA),
      # end and and start fill inn due to windowed operation
      baseline = na.fill(baseline, "extend")
    ) %>%
    # signal drop or clogging warnings:
    # verify(baseline > median_counts * 0.75 - 2) %>%
    # verify(baseline < median_counts * 1.25 + 2) %>%
    # Poisson threshold for filtering particle events, minimum of 5 counts.
    # (Could consider using qpois + baseline - mean_baseline, yet may
    mutate(
      h_thr = rollmax(baseline, 400, by = 50, fill = NA),
      h_thr = na.fill(h_thr, "extend"),
      h_thr = qpois(
        p = (1 - 0.05 / 600000),
        lambda = if_else(h_thr > 0.2, h_thr, 0.2)
      )
    )
}


# Discriminates and calculates peak parameters, numbers.
# - in: m x 1 tibble with count data.
# - out: p x 5 tibble with peak_n, peak_time, peak_max, peak_width and peak_area
sp_peak_discriminator <- function(sp_read) {
  peaks <- sp_baseline_thr(sp_read) %>%
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
      ),
      # labeling by n event, classify events whether max is above thr or not.
      peak_n = cumsum(peak_n),
      peak_time = row_number()
    ) %>%
    group_by(peak_n) %>%
    # Peak width error(?)
    mutate(
      peak_max = max(counts),
      peak_width = sum(c_above_baseline_thr > 0),
      above_height_thr = peak_max > h_thr,
      # integration using group by.
      # interpolation when passing baseline could perhaps improve accuracy
      # across size ranges at the expense of computation time.
      # This because the aspect ratio may increase w size.
      peak_area = sum(c_above_baseline_thr) %>% round(digits = 0)
      )
    %>%
    # # for signal distribution:
    # group_by(peak_area, above_height_thr) %>%
    # mutate(
    #   peak_area_count = length(unique(peak_n))
    # ) %>%
    ungroup() %>%
    select(!c_above_baseline_thr) %>%
    distinct(peak_n, .keep_all = TRUE) %>%
    filter(above_height_thr == TRUE) %>%
    select(peak_n, peak_time, peak_max, peak_width, peak_area)

  peaks
}


# adds peak data to datafiles.
# in: n x 5, sample_name, dataset, datafile, isotope, type
# out: tibble with nested peak data (n x 7, nested: m x 5)
sp_particles <- function(classified) {
  # ideas for optimization(?): https://www.brodrigues.co/blog/2021-03-19-no_loops_tidyeval/
  classified %>%
    mutate(peaks = future_map(filepath, ~ sp_reader(.) %>% sp_peak_discriminator()))
}


# calibration ####

# calculates intercept, response from calibration curves. Used within sp_calibration.
# - in: classified tibble
# - out: tibble showing STDs used and response, intercept, r2.
# - todo: error checking: minimum 2 points, non-negative values, all isotopes present.
sp_response <- function(classified) {
  tryCatch(
    expr = {
      classified %>%
        filter(
          type == "STD",
          str_detect(sample_name, isotope)
        ) %>%
        rowwise() %>%
        mutate(
          mean_counts = read_csv(filepath,
            skip = 4,
            col_types = "d", col_select = 2
          ) %>%
            drop_na() %>%
            rename(counts = 1) %>%
            as_vector() %>%
            mean(),
          conc_ppb = str_extract(sample_name, str_extract(
            sample_name,
            paste0("[\\d.]+(?=[\\sa-zA-Z]*", isotope, ")")
          )) %>% as.numeric()
        ) %>%
        group_by(isotope) %>%
        mutate(
          intercept =
            coef(
              lm(formula = as.numeric(mean_counts) ~ conc_ppb)
            )[1],
          response =
            coef(
              lm(formula = as.numeric(mean_counts) ~ conc_ppb)
            )[2],
          r2 = summary(
            lm(formula = as.numeric(mean_counts) ~ conc_ppb)
          )$r.squared
        ) %>%
        ungroup()
    },
    error = function(e) {
      print(
        sprintf(
          "An error occurred in foo at %s : %s",
          Sys.time(),
          e
        )
      )
    }
  )
}

# calculates signal per mass for RM. Used in sp_calibration, needs sp_peak_discriminator and sp_reader.
# - in: classified tibble
# - out: single value mass per count.
# - todo:
#     - error checking: certain range.
sp_mass_cal <- function(classified,
                           RM_dia = 60,
                           RM_density = 19.32,
                           RM_isotope = "Au",
                           element_fraction = 1) {
  tryCatch(
    expr = {
      
        RM_areas <- classified %>%
          filter(
            type == "RM",
            str_detect(sample_name, isotope)
          ) %>%
          pull(filepath) %>% sp_reader() %>% 
          sp_peak_discriminator()
        
        RM_area_mode <- density(
          RM_areas$peak_area
        )$x[which.max(density(RM_areas$peak_area)$y)]
        
        # Get the mass per count by relating the "known" RM mass to counts_RM
        mass_count <- # counts per kg RM element
          (
            RM_density * # density RM [g/cm3]
              1000 * (4 / 3) * pi * (RM_dia * 10^(-9) / 2)^3 * # RM volume
              (
                element_fraction
              )
          ) /
          RM_area_mode # element fraction divided by RM peak area "mode"
        
      return(mass_count)
    },
    error = function(e) {
      print(
        sprintf(
          "An error occurred in sp_mass_cal at %s : %s",
          Sys.time(),
          e
        )
      )
    }
  )
}

sp_conc_cal <- function(classified,
                        RM_dia = 60,
                        RM_density = 19.32,
                        RM_isotope = "Au",
                        element_fraction = 1,
                        RM_conc = 0) {
  tryCatch(
    expr = {
      
      RM_areas <- classified %>%
        filter(
          type == "RM",
          str_detect(sample_name, isotope)
        ) %>%
        pull(filepath) %>% sp_reader() %>% 
        sp_peak_discriminator()
      
      
      # Get the mass per count by ... 
      mass_count <- # counts per kg RM element
 
      
      return(mass_count)
    },
    error = function(e) {
      print(
        sprintf(
          "An error occurred in sp_mass_cal at %s : %s",
          Sys.time(),
          e
        )
      )
    }
  )
}



# outputs calibration data: detector flow rate, mass per signal, response.
# - in: classified
# - out: tibble w detector flow rate, mass per signal, response for each isotope
sp_calibration <- function(classified,
                     RM_dia = 60,
                     RM_density = 19.32,
                     RM_isotope = "Au",
                     element_fraction = 1) {
  tryCatch(
    expr = {
      mass_signal_RM <- sp_mass_cal(classified,
        RM_dia = 60,
        RM_density = 19.32,
        RM_isotope = "Au",
        element_fraction = 1
      ) # for RM

      responses <- sp_response(classified)

      response_RM <- responses %>%
        filter(isotope == RM_isotope) %>%
        distinct(response) %>%
        as.numeric()

      calib <- responses %>%
        distinct(isotope, response, intercept, r2) %>%
        mutate(
          mass_signal = mass_signal_RM * response_RM / response,
          detector_flow_rate = mass_signal_RM * response_RM * 10^9 # L per dwell
        )

      return(calib)
    },
    error = function(e) {
      print(
        sprintf(
          "An error occurred in sp_calib at %s : %s",
          Sys.time(),
          e
        )
      )
    }
  )
}

tic("sp_calibration")
calibrated <- sp_calibration(classified)
toc("sp_calibration")




# adds mass to peaks data
# in: tibble w nested peak data (n x 7, nested: m x 5)
# out: tibble w nested peak data (n x 6, new nested: m x 6)
sp_particle_mass <- function(sp_particled, calibrated) {
  mass <- sp_particled %>%
    mutate(particles = future_map2(
      peaks, isotope,
      ~ .x %>% mutate(particle_mass = peak_area *
        calibrated %>%
          filter(isotope == .y) %>%
          pull(mass_signal) %>%
          as.numeric())
    ))
}


#bruke mutate = map(peaks) og en tibble av fns for å generere all output?
  
sp_summary() <- function(sp_particle_massed, ){}

# Output fn ####
# use sp_particles and sp_calibration output to get sizes. Output nParticle summariser type dataframe, with mass distributions.
# also take in dens_comps type dataframe to get sizes.
# keep also filepath for possibility of validating the spectrum.

# Data:
# - n particles
# - particle number concentration
# - particle mass concentration
# - ionic concentration (sum counts - sum peak areas) / n_dwells
# - Median mass, kde mass
# - Median size, kde size
# - TE
# - Particle dataframe:
#     - Masses, sizes, times, width, max, 



# TODO ####
# add calibrations to get mass and sizes
# Option for concentration calibration.

# rmd for browsing?


# - sum counts
# - ionic counts (sum counts - sum peak areas )/ n_dwells


  
