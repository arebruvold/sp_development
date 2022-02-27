
library(zoo) #custom rolling windows
library(tidyverse) #backbone
library(RcppRoll)
library(furrr) #parallel computing
library(assertr)
plan(multisession, workers = 14)
library(tictoc)


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
  INCOUNT[INCOUNT>5] <- 5 #should be quicker than code above
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

## sp calibration####

# finds count files only
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
        sprintf("An error occurred in foo at %s : %s",
                Sys.time(),
                e)
      )
    })
}


# better implementation: list.files into vector pipe into a map with a function that read lines, pipe into new map make function to take each lines and get the info I want with regex mutate.
# classifies all found cps count files into RM, STD or SAMPLE
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



sp_response(classified) %>% View()


sp_mass_signal <- function(classified,
                           RM_dia = 60,
                           RM_density = 19.32,
                           RM_isotope = "Au",
                           element_fraction = 1) {
  tryCatch(
    expr = {
      #
      RM_areas <- classified %>%
        filter(
          type == "RM",
          str_detect(sample_name, isotope)
        ) %>%
        pull(filepath) %>% sp_reader() %>% 
        sp_peaks()

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
          "An error occurred in sp_mass_signal at %s : %s",
          Sys.time(),
          e
        )
      )
    }
  )
}


sp_calib <- function(classified,
                     RM_dia = 60,
                     RM_density = 19.32,
                     RM_isotope = "Au",
                     element_fraction = 1) {
  tryCatch(
    expr = {
      mass_signal_RM <- sp_mass_signal(classified,
        RM_dia = 60,
        RM_density = 19.32,
        RM_isotope = "Au",
        element_fraction = 1
      ) # for RM

      responses <- sp_response(classified) #

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

## sp peak discrimination ####

hei <- sp_calib(classified) %>% View()

CPS_csv <- "~/sp-data/dev-sp/008_RM-197-197_count_1.csv"

# raw count file reader
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


hoi <- sp_reader("~/sp-data/dev-sp/008_RM-197-197_count_1.csv")

reshaped <- hoi

baseline_thr_finder <- function(reshaped) {
  median_counts <- median(reshaped$counts, na.rm = TRUE)

  bl_thr <- reshaped %>%
    drop_na() %>%
    # robust baseline with rolling median or kde
    mutate( # baseline = rollapply(counts, 301, dens_max, by = 50, fill = NA),
      baseline = roll_median(counts, n = 301, by = 50, fill = NA),
      # end and and start fill inn due to windowed operation
      baseline = na.fill(baseline, "extend")
    ) %>%
    #  signal drop or clogging warning
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


# - peak location (time start or max peak time)

sp_peaks <- function(sp_read) {
  peaks <- baseline_thr_finder(sp_read) %>% 
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
      peak_max = max(counts),
      #peak_width = n(), wrong, counts until next signal..
      peak_width = sum(c_above_baseline_thr > 0),
      above_height_thr = peak_max > h_thr 
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
    select(!c_above_baseline_thr) %>% 
    distinct(peak_n, .keep_all = TRUE) %>% 
    filter(above_height_thr == TRUE) %>% 
    select(peak_n, peak_max, peak_width, peak_area)
  
  peaks
}


# see: https://www.brodrigues.co/blog/2021-03-19-no_loops_tidyeval/ quo?
sp_particles <- function(classified) {
  classified %>%
    mutate(peaks = future_map(filepath, ~ sp_reader(.) %>% sp_peaks()))
}

# TODO ####
# add calibrations to get mass and sizes
# rmd for browsing?


# - sum counts
# - ionic counts (sum counts - sum peak areas )/ n_dwells





# using mass per count AND sp response
  
# keep for now, needed for sp_mass signal: ####
  
# nParticle_finder_blocks <- function(CPS_csv, prop = 1) {
#   CPS_csv %>%
#     folder_rawfinder() %>%
#     future_map_dfr(
#       ~ .x %>%
#         # read in raw (cps) file
#         raw_reshaper(prop = prop) %>%
#         # calculate baseline
#         baseline_thr_finder() %>%
#         # integrate peaks
#         peak_integration(),
#       .options = furrr_options(seed = TRUE)
#     )
# }
  
