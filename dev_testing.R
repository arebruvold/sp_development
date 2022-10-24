# Development testing



test <- sp_outputer(peaked, calibrated, acq_time = 60, dens_comps = dens_comps,RM_isotope = "Au", sample_intake_rate = 0.311)

test_nas <- sp_wrapper("~/sp-data/22_102021_aersol_gas_dilution_optimization/22_1021_Aerosol_gas_tuning_Au_d000_n105.b/", RM_string = "Au RM 60 nm 200ng/L", dens_comps = dens_comps,RM_isotope = "Au", sample_intake_rate = 0.311)


test_working <- sp_wrapper("~/sp-data/22_102021_aersol_gas_dilution_optimization/22_1020_Aerosol_gas_tuning_04_060.b/", RM_string = "Au RM 60 nm 200ng/L", dens_comps = dens_comps,RM_isotope = "Au", sample_intake_rate = 0.311)

class_working <- sp_classifier("~/sp-data/22_102021_aersol_gas_dilution_optimization/22_1020_Aerosol_gas_tuning_04_060.b/", RM_string = "Au RM 60 nm 200ng/L", dens_comps = dens_comps,RM_isotope = "Au", sample_intake_rate = 0.311)
