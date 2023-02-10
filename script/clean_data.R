# clean_data.R
# February 9, 2023
# Cleaning data from CA Open Data portal

library(data.table)

# Case counts
cases.dat <- fread('data/cases.csv')

# Vaccination counts
vax.dat <- fread('data/vax.csv')

# Merge data
covid19.dat <- merge(cases.dat[area_type == 'County',.(date, county = area, cases, deaths)],
			vax.dat[,.(date = administered_date, county, fully_vaccinated)],
			by = c('date', 'county'))


san_francisco.dat <- covid19.dat[county == "San Francisco",]
setorder(san_francisco.dat, date)
write.csv(san_francisco.dat, "data/san_francisco.csv", row.names = F)
