#Data

Contains file with raw earthquake event timing data. The code is designed to read data in mutliple formats, as different papers in the literature will present their results differently. .txt files with filenames ending in 'simple' imply only basic information about the uncertainty distribution of the event dates is given (mean, standard deviation, 95% confidence bounds etc). .csv file contain full OxCal output describing the posterior pdfs of event ages. Data should be parsed in conjunction with the information in the associated parameter file in `../params/`.

Formats for 'simple' files can be:
`Date1	Date2`
which gives upper and lower uncertainty bounds as calendar dates (positive dates are AD, negative BC);
or
`Date	Uncertainty`
which gives a mean date and some measure of uncertainty (1 or 2 sigma);
or
`Age1	Age2`
which gives upper and lower uncertainty bounds as ages before present (i.e. before 1950);
or
`Age	Uncertainty`
which gives a mean age and some measure of uncertainty (1 or 2 sigma).

An optional third column `Certain` can be added. This contains a 1 or 0 for each row depending on whether the event is certain (1) or uncertain (1).

Whether the uncertainty ranges represent 1 or two sigma is defined by the parameter `sigma_level` in the parameter file. 
