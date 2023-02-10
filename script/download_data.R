# download_data.R
# February 9, 2023
# Downloading data from CA Open Data portal

library(data.table)

# Case counts
cases.url <- paste0(
	"https://data.chhs.ca.gov/dataset/",
	"f333528b-4d38-4814-bebb-12db1f10f535/",
	"resource/",
	"046cdd2b-31e5-4d34-9ed3-b48cdbc4be7a/",
	"download/covid19cases_test.csv"
)
cases.dat <- fread(cases.url)
write.csv(cases.dat, 'data/cases.csv', row.names = F)

# Vaccination counts
vax.url <- paste0(
	"https://data.chhs.ca.gov/dataset/",
	"e283ee5a-cf18-4f20-a92c-ee94a2866ccd/",
	"resource/",
	"130d7ba2-b6eb-438d-a412-741bde207e1c/",
	"download/covid19vaccinesbycounty.csv"
)
vax.dat <- fread(vax.url)
write.csv(vax.dat, 'data/vax.csv', row.names = F)
