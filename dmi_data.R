source("libs_and_funcs.R")

#Get data from DMI metObs REST API

#Parameters
stat <- "06174"
limit <- "30000"

date_range1 <- "2019-08-01T00:00:00Z/2019-09-30T00:00:00Z"
date_range2 <- "2020-03-01T00:00:00Z/2020-06-01T00:00:00Z"

#Base URL
base <- paste0("https://dmigw.govcloud.dk/v2/metObs/collections/observation/items?limit=", limit)

#First request
request1 <- paste0(base, "&parameterId=wind_speed", "&stationId=", stat, "&datetime=", date_range1, "&api-key=", dmi_key)
response1 <- GET(request1)
json1 <- content(response1, as="text") 

#Second request
request2 <- paste0(base, "&parameterId=wind_speed", "&stationId=", stat, "&datetime=", date_range2, "&api-key=", dmi_key)
response2 <- GET(request2)
json2 <- content(response2, as="text") 

#Time zone is UTC+02
df1 <- fromJSON(json1)$features$properties %>% 
  mutate(datetime = ymd_hms(observed) - 2*60*60) %>% 
  select(datetime, wnd = value) %>% 
  arrange(datetime)

df2 <- fromJSON(json2)$features$properties %>% 
  mutate(datetime = ymd_hms(observed) - 2*60*60) %>% 
  select(datetime, wnd = value) %>% 
  arrange(datetime)

write_csv(df1, "data/wind_1.csv")
write_csv(df2, "data/wind_2.csv")
