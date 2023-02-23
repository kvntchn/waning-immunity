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

san_francisco.dat <- fread("data/san_francisco.csv")
san_francisco.dat <- san_francisco.dat[date >= as.Date("2021-06-15"),]
san_francisco.dat[,time := 1:.N]
san_francisco.dat[, cases := zoo::rollmean(cases, 7, na.pad = T, align = 'left')]
san_francisco.dat <- san_francisco.dat[!is.na(cases),]
write.csv(san_francisco.dat, "data/san_francisco_clean.csv", row.names = F)
