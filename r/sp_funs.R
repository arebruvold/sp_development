

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

# sp functions ####

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



classified <- sp_classifier("~/sp-data/dev-sp/", RM_string = "AuRM500ng/LM")

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
        pull(filepath) %>%
        nParticle_finder_blocks() %>%
        filter(
          above_height_thr == TRUE
        ) %>%
        group_by(peak_n) %>%
        summarise(RM_area = mean(peak_area))

      RM_area_mode <- density(
        RM_areas$RM_area
      )$x[which.max(density(RM_areas$RM_area)$y)]

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

sp_mass_signal(classified)

# calculate mass flow rate and thus TE ####
# using mass per count AND sp response
  
# Keep for now: ####
  
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
  
