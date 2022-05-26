
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

#used for size calculations, to be replaced/ combined with mass_fractioner.
dens_comps <- tibble(
  "isotope" = c("RM", "Au", "Al", "Mn", "Pb", "Fe", "Si", "Ti", "Cr", "Ce", "Zr", "Cu", "Cd", "Ba", "Co", "Ni", "Zn"),
  "density" = c(19.32, 19.32, 2.56, 4.25, 6.29, 4.30, 2.65, 4.17, 5.22, 7.22, 5.68, 6.31, 2, 2, 2, 2, 2),
  "element_fraction" = c(1, 1, 26.982/(39.098+26.982+3*28.085+8*15.999), 1 / 1.6, 1 / 1.465, 1 / 1.59, 1 / 2.139, 1/1.67, 1/1.462, 1/1.23, 91.22/(91.22+16.00*2), 63.55/(63.55+16), 1, 1, 1, 1, 1)
)

#not implemented yet.
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



# classified <- sp_classifier("~/sp-data/Fordefjorden/", RM_string = "AuRM500ng/LM")


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
    ) %>%
    # # for signal distribution (REMOVED):
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
# out: n x 6, adds nested p x 5: peak_n, peak_time. peak_max, peak_width, peak_area  
sp_peaker <- function(classified) {
  # ideas for optimization(?): https://www.brodrigues.co/blog/2021-03-19-no_loops_tidyeval/
  classified %>%
    mutate(
      peaks = future_map(filepath, ~ sp_reader(.) %>% sp_peak_discriminator()),
      mean_counts = future_map(filepath, ~ sp_reader(.) %>% pull(counts) %>% mean())
    )
}


## calibration ####

# calculates intercept, response from calibration curves. Used within sp_calibration.
# - in: n x 6, filepath, sample_name, dataset, datafile, isotope, type
# - out: s x 11, STDs used and conc_ppb, mean counts, response, intercept, r2.
# - todo: error checking: minimum 2 points, non-negative values, all isotopes present, minimum R.
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
          "An error occurred in sp_response at %s : %s",
          Sys.time(),
          e
        )
      )
    }
  )
}

# # calculates signal per mass for RM. Used in sp_calibration, needs sp_peak_discriminator and sp_reader.
# # - in: n x 6 ("classified"), filepath, sample_name, dataset, datafile, isotope, type
# # - out: single value mass per count.
# # - todo: error checking, e.g. unique RM, NPs present in RM, certain range.
# sp_mass_signal_rmer <- function(classified,
#                            RM_dia,
#                            RM_density,
#                            RM_isotope,
#                            element_fraction) {
#   tryCatch(
#     expr = {
#       
#         RM_areas <- classified %>%
#           filter(
#             type == "RM",
#             str_detect(sample_name, isotope)
#           ) %>%
#           pull(filepath) %>% sp_reader() %>% 
#           sp_peak_discriminator()
#         
#         RM_area_mode <- density(
#           RM_areas$peak_area
#         )$x[which.max(density(RM_areas$peak_area)$y)]
#         
#         # Get the mass per count by relating the "known" RM mass to counts_RM
#         mass_count <- # counts per kg RM element
#           (
#             RM_density * # density RM [g/cm3]
#               1000 * (4 / 3) * pi * (RM_dia * 10^(-9) / 2)^3 * # RM volume
#               (
#                 element_fraction
#               )
#           ) /
#           RM_area_mode # element fraction divided by RM peak area "mode"
#         
#       return(mass_count)
#     },
#     error = function(e) {
#       print(
#         sprintf(
#           "An error occurred in sp_mass_signal_rmer at %s : %s",
#           Sys.time(),
#           e
#         )
#       )
#     }
#   )
# }

# outputs calibration data: detector flow rate, mass per signal, response.
# - in: n x 6 ("classified"), filepath, sample_name, dataset, datafile, isotope, type
# - out: tibble w detector flow rate, mass per signal, response for each isotope
sp_calibrator <- function(classified,
                     RM_dia,
                     RM_density,
                     RM_isotope,
                     element_fraction) {
  tryCatch(
    expr = {
      #calculates mass per count using known RM mass and relating to KDE max peak area.
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
      
      mass_signal_RM <- # counts per kg RM element
        (
          RM_density * # density RM [g/cm3]
            1000 * (4 / 3) * pi * (RM_dia * 10^(-9) / 2)^3 * # RM volume
            (
              element_fraction
            )
        ) /
        RM_area_mode

      responses <- sp_response(classified)

      response_RM <- responses %>% 
        filter(isotope == RM_isotope) %>%
        distinct(response) %>%
        as.numeric()

      calib <- responses %>%
        distinct(isotope, response, intercept, r2) %>%
        mutate(
          mass_signal = mass_signal_RM * response_RM / response, # counts per kg element
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

## output ####

sp_outputer <- function(peaked, 
                        calibration_data,
                        acq_time,
                        dens_comps,
                        RM_isotope,
                        sample_intake_rate) {
  peaked %>%
    ungroup() %>%
    mutate(
      peaks = future_map2(
        peaks, isotope,
        ~ .x %>% mutate(
          particle_mass = peak_area *
              (filter(calibration_data, isotope == .y) %>%
              pull(mass_signal) %>%
              as.numeric()),
          particle_size = (6 * particle_mass /
            (pi * 
              (filter(dens_comps, isotope == .y) %>%
              pull(density, element_fraction) %>%
              Reduce(`*`, .) %>%
              as.numeric()) * 1000))^(1 / 3) * 10^9
        )
      ),
      summary_data = future_map2(
        peaks, isotope,
        ~ .x %>%
          summarise(
            n_particles = n(),
            mean_size = mean(particle_size),
            median_size = median(particle_size),
            mass_conc = sum(particle_mass) * 10^12 / (calibration_data %>%
              filter(isotope == RM_isotope) %>%
              pull(detector_flow_rate) %>%
              as.numeric() * acq_time * 10000), # ng/L
            particle_conc = n_particles / (calibration_data %>%
              filter(isotope == RM_isotope) %>%
              pull(detector_flow_rate) %>%
              as.numeric() * acq_time * 10000), # #/L
            detector_flow_rate = calibration_data %>%
              filter(isotope == RM_isotope) %>%
              pull(detector_flow_rate) %>%
              as.numeric(),
            transport_efficiency = detector_flow_rate * 1000 * 60 * 10000 / sample_intake_rate # L/dwell to ml/min conversion
          )
      ), # Richard: why does this not work?
      # total_conc = calibration_data %>%
      #   filter(isotope == isotope) %>%
      #   pull(mass_signal, detector_flow_rate) %>%
      #   Reduce(`*`, .) %>%
      #   as.numeric() * mean_counts %>% as.numeric()
      total_conc = future_map2(
        mean_counts, isotope,
        ~ .x / (calibration_data %>%
          filter(isotope == .y) %>%
          pull(response) %>%
          Reduce(`*`, .))
      )
    ) %>% unnest(summary_data)
}

sp_wrapper <- function(csv_folder,
                       acq_time = 60,
                       sample_intake_rate = 0.346,
                       RM_string = "RM",
                       RM_dia = 60,
                       RM_isotope = "Au",
                       RM_density = 19.32,
                       element_fraction = 1,
                       dens_comps = dens_comps) {
  
  classified <- csv_folder %>%
    sp_classifier(RM_string)
  
  peaked <- sp_peaker(classified)
  
  calibrated <- sp_calibrator(classified,
                              RM_dia,
                              RM_density,
                              RM_isotope,
                              element_fraction
  )
  
  sp_output <- sp_outputer(peaked, 
                           calibrated,
                           acq_time,
                           dens_comps,
                           RM_isotope,
                           sample_intake_rate
  )
  
  return(sp_output)
}

## validation ####
test %>%
  filter(reduce(map(time_intervals, near, x= row_number(), tol = 2),`|`))

# Takes in sample file and the particle output, selects the smallest particles in terms of area and max intensity and reveals them in the spectrum with context.
# Todo:
#   - change to smallest 10 particles
#   - show context of particles, 
sp_spectrum_validation <- function(filepath, sp_output) {
  peaks_data <- sp_output %>%
    filter(filepath == filepath) %>%
    pluck("peaks") %>% map_df(.f = as_tibble) 
  # %>% 
  #   filter(quantile(peak_area, n_prop) > peak_area)
  
  time_intervals <- peaks_data %>% pull(peak_time)
  
raw_data <- sp_reader(filepath) %>%
  filter(reduce(map(time_intervals, near, x = row_number(), tol = 50), `|`))


  ggplot(data = raw_data, aes(as.factor(row.names(raw_data)), counts)) +
    geom_line()
  
  
  # +
  #   geom_point(data = sp_output %>%
  #     filter(filepath == filepath) %>% pluck("peaks"), aes(x = peak_time, y = 0), size = 4, color = "green")
}
  


# signal_found_plotter <- function(foundParticles,
#                                  prop = 1,
#                                  t_0 = 0, t_e = 600,
#                                  i_min = 0,
#                                  i_max = Inf,
#                                  area_below = Inf) {
#   if (ncol(foundParticles) <= 16) {
#     plot_Found <- foundParticles %>%
#       group_by(isotope, datafile) %>%
#       slice_head(prop = prop) %>%
#       ungroup() %>%
#       filter(
#         Time > t_0,
#         Time < t_e,
#         peak_max_I >= i_min,
#         peak_max_I <= i_max
#       ) %>%
#       ggplot(aes(x = Time, y = counts)) +
#       # geom_line(size = 0.2) +
#       theme_minimal() +
#       geom_line() +
#       geom_line(aes(x = Time, y = baseline), color = "blue", size = 2) +
#       geom_line(aes(x = Time, y = h_thr), color = "red") +
#       geom_point(
#         data = foundParticles %>% group_by(isotope, datafile) %>%
#           slice_head(prop = prop) %>% distinct(peak_n, .keep_all = TRUE) %>%
#           mutate(Time_event = if_else(above_height_thr == TRUE &
#                                         peak_area < area_below, Time, NA_real_)),
#         aes(
#           x = Time_event,
#           y = 0
#         ),
#         #   label = round(Time_event, digits = 2),
#         #   angle = 90,
#         #   hjust = 1
#         # ),
#         size = 3,
#         color = "green"
#       ) +
#       # geom_text(
#       #   data = foundParticles %>% group_by(isotope, datafile) %>%
#       #     # Use FILTER instead of the mutate
#       #     slice_head(prop = prop) %>% distinct(peak_n, .keep_all = TRUE) %>%
#       #     filter(above_height_thr == TRUE,
#       #            peak_area < area_below),
#       #   aes(
#       #     x = Time,
#       #     y = 0,
#       #     label = round(peak_area, digits = 0),
#       #     angle = 90
#     #   ),
#     #   size = 3,
#     #   color = "green"
#     # ) +
#     labs(title = "Signal processing") +
#       # # labels/text, OLD, NEED HELP WITH THIS:
#       # labs(subtitle = paste(
#       #   foundParticles %>% filter(above_height_thr == "Yes") %>%
#       #     distinct(peak_n) %>% nrow(),
#       #   " particle events",
#       #   ", min intensity threshold: ", min(foundParticles$h_thr),
#       #   "\n", foundParticles$isotope[1], ":  ",
#       #   foundParticles$datafile, ":  ",
#       #   foundParticles$dataset,
#       #   sep = ""
#       # )) +
#     facet_wrap(isotope ~ sample_name, scales = "free")
#   }
#   
#   if (ncol(foundParticles) > 16) {
#     plot_Found <- foundParticles %>%
#       group_by(isotope, datafile) %>%
#       slice_head(prop = prop) %>%
#       ungroup() %>%
#       filter(
#         Time > t_0,
#         Time < t_e,
#         peak_max_I >= i_min,
#         peak_max_I <= i_max
#       ) %>%
#       ggplot(aes(x = Time, y = counts)) +
#       # geom_line(size = 0.2) +
#       theme_minimal() +
#       geom_line() +
#       geom_line(aes(x = Time, y = baseline), color = "blue", size = 2) +
#       geom_line(aes(x = Time, y = h_thr), color = "red") +
#       # geom_point(
#       #   data = foundParticles %>% group_by(isotope, datafile) %>%
#       #     slice_head(prop = prop) %>%  distinct(peak_n, .keep_all = TRUE) %>%
#       #     mutate(Time_event = if_else(above_height_thr == "Yes", Time, NA_real_)),
#       #   aes(
#       #     x = Time_event,
#       #     y = 2)
#       #   ,
#       #   #   label = round(Time_event, digits = 2),
#       #   #   angle = 90,
#       #   #   hjust = 1
#     #   # ),
#     #   size = 2,
#     #   color = "green") +
#     geom_text(
#       data = foundParticles %>% group_by(isotope, datafile) %>%
#         # Use FILTER instead of the mutate
#         slice_head(prop = prop) %>% distinct(peak_n, .keep_all = TRUE) %>%
#         filter(above_height_thr == TRUE),
#       aes(
#         x = Time,
#         y = -5,
#         label = round(size_nm, digits = 0),
#         angle = 90
#       ),
#       size = 3,
#       color = "green"
#     ) +
#       labs(title = "Signal processing") +
#       # # labels/text, OLD, NEED HELP WITH THIS:
#       # labs(subtitle = paste(
#       #   foundParticles %>% filter(above_height_thr == "Yes") %>%
#       #     distinct(peak_n) %>% nrow(),
#       #   " particle events",
#       #   ", min intensity threshold: ", min(foundParticles$h_thr),
#       #   "\n", foundParticles$isotope[1], ":  ",
#       #   foundParticles$datafile, ":  ",
#       #   foundParticles$dataset,
#       #   sep = ""
#       # )) +
#     facet_wrap(isotope ~ sample_name, scales = "free")
#   }
#   plot_Found
# }

# TODO ####
# - add comments
# - make output function work yet produce NAs if calibration and/or dens_comps are empty.
# - validate
# - increase speed, classifier especially.
# - error testing using assertr package.
# - add diagnostic spectrum plotter, combine with identifying particles by time signature.
#     - extract vector from peaks data and feed to tis function, selectively plotting only the ones chosen.
# - Make convenient rmd or shiny app/ setup server for automated data processing.
# - Readme, note only same acq time per isotope.
# error if multiple RMs.




  
