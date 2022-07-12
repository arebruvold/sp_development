source("r/sp_funs.R")

#Check if a continuous poisson with noise and a poiosson with higher rate  occuring for some subsequent dwells results in similar plots as shown by Schardt et al.

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

test_sign <- signal_simulator(0.001, 100000) %>%
  mutate(pstart = rpois(100000, 0.01)) %>%
  mutate(row_n = row_number())

# now make the particle simulator, generating only signals on a particle event times, though then for consecutive dwells.


# conditional mutate extract row numbers that have a poisson process ion cloud event
# extract vector of particle events
# # use near to check if current row is enar any of the values in the vector

p_rows <- test_sign %>%
  mutate(pstart = rpois(100000, 0.01)) %>%
  filter(pstart == 1) %>% 
  pull(row_n)

test_sign %>% filter(near(row_n, )
)




