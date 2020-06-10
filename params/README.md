# Params

Contains text files with information concering the earthquake chronology data that is used to read in the data in `../data/`

Format is:

`filename = ` Relative path to data file

`chron_type = ` Can be `Age2Sigma` for simple text file format or `OxCal` is OxCal output is used.

`sigma_level = ` Can be `1` or `2` depending on whether the uncertainties provided in by data in the simple text format represent 1 or 2 sigma.

`event_order = ` Described the order in which event ages are read in from the input data. Can be 'Forwards` or `Backward` in the case of simple text data format. If OxCal format is used, should be a list of event names relating to OxCal output defining chronological order, e.g. `['D', 'C', 'B', 'A'] `

`tectonic_region = ` One of `Intraplate_noncratonic`, `Intraplate_cratonic`, `Near_plate_boundary`, `Plate_boundary_network`, `Plate_boundary_master`.

`faulting_style = ` One of `Reverse`, `Normal` or `Strike_slip`.

`slip_rate_mean_lower_upper = ` List of mean, lower and upper bound estimates of slip-rates in mm/yr. Bounds are treated as 95% intervals. e.g. `[1.0, 0.5, 1.5]`
