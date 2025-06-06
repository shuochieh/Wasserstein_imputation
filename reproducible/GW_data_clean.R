library(rmatio)
library(lubridate)
library(gridExtra)
library(maps)
library(ggplot2)

dta = read.mat("./real_data/gw_temp/groundwater/gw_2020_67.mat")
dta = dta[[1]]

data_span = matrix(0, nrow = length(dta), ncol = 6)
for (i in 1:length(dta)) {
  GMTtime = date_decimal(dta[[i]][[5]], "GMT")
  cat("Series", i, ": from", year(head(GMTtime, 1)), "to", 
      year(tail(GMTtime, 1)), "\n")
  data_span[i,] = c(year(head(GMTtime, 1)), month(head(GMTtime, 1)), day(head(GMTtime, 1)),
                    year(tail(GMTtime, 1)), month(tail(GMTtime, 1)), day(tail(GMTtime, 1)))
}

check_sufficient = function (dta, yr, m, threshold = 480) {
  n = length(dta[[5]])
  date_info = matrix(0, nrow = n, ncol = 3)
  dates = date_decimal(dta[[5]], "GMT")
  
  idx = which(year(dates) == yr)
  idx = intersect(which(month(dates) == m), idx)

  suff = length(idx) > threshold
  
  if (suff) {
    return (list("suff" = suff, "idx" = sort(idx), "mean" = mean(dta[[6]][idx])))
  }
  
  return (list("suff" = suff, "idx" = sort(idx)))
}

check_and_fill = function (dta, yr, m, threshold = 480) {
  temp = check_sufficient (dta, yr, m, threshold)
  
  if (temp$suff) {
    return (temp$mean)
  } else {
    return (NA)
  }
}

## Get monthly data
apply(data_span, FUN = min, MARGIN = 2) 
apply(data_span, FUN = max, MARGIN = 2) 

monthly_data = matrix(NA, nrow = length(dta), ncol = ((2020 - 1992 + 1) * 12))
for (i in 1:length(dta)) {
  counter = 1
  for (yr in 1992:2020) {
    # res = rep(NA, 12)
    for (m in 1:12) {
      # res[m] = check_and_fill(dta[[i]], yr, m)
      monthly_data[i, counter] = check_and_fill(dta[[i]], yr, m)
      cat("   ", yr, "/", m, "\n")
      counter = counter + 1
    }
  }
  cat("Series", i, "done\n")
}

colSums(is.na(monthly_data))

# discard first 9 months and the last 4 months (no series has data)
# Essentially, monthly data from 1992/10 to 2020/8
gw_dta = monthly_data[,-c(1:9, 345:348)]
n_missing = rowSums(is.na(gw_dta))

select_stations = which(n_missing <= 127)

# Plot the stations geographic locations
long = rep(NA, length(select_stations))
lat = rep(NA, length(select_stations))
for (ii in 1:length(select_stations)) {
  i = select_stations[ii]
  long[ii] = dta[[i]][[3]]
  lat[ii] = dta[[i]][[4]]
}

loc <- data.frame(
  long = long,
  lat = lat
)
taiwan_map <- map_data("world", region = "Taiwan")
ggplot() +
  geom_polygon(data = taiwan_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "black") +
  xlim(120,122) +
  geom_point(data = loc, aes(x = long, y = lat), color = "brown3", size = 1.5) +
  # geom_vline(xintercept = c(121), linetype = "dashed", color = "blue") +
  # geom_hline(yintercept = c(24.9, 22.9), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "", x = "Longitude", y = "Latitude")

# Define regional stations
north_stations = intersect(which(long > 121), which(lat > 24.9))
east_stations = intersect(which(long > 121), which(lat < 24.9))
west_stations = intersect(which(long < 121), which(lat > 22.9))
south_stations = which(lat < 22.9)

nloc <- data.frame(
  long = long[north_stations],
  lat = lat[north_stations]
)
eloc <- data.frame(
  long = long[east_stations],
  lat = lat[east_stations]
)
wloc <- data.frame(
  long = long[west_stations],
  lat = lat[west_stations]
)
sloc <- data.frame(
  long = long[south_stations],
  lat = lat[south_stations]
)
taiwan_map <- map_data("world", region = "Taiwan")
ggplot() +
  geom_polygon(data = taiwan_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "black") +
  xlim(120,122) + 
  geom_point(data = nloc, aes(x = long, y = lat, color = "North"), size = 1.5) +  
  geom_point(data = wloc, aes(x = long, y = lat, color = "West"), size = 1.5) +   
  geom_point(data = sloc, aes(x = long, y = lat, color = "South"), size = 1.5) +  
  geom_point(data = eloc, aes(x = long, y = lat, color = "East"), size = 1.5) +   
  scale_color_manual(values = c("North" = "blue3", "West" = "brown3", "South" = "magenta", "East" = "seagreen4")) +
  theme_minimal() +
  labs(title = "", x = "Longitude", y = "Latitude", color = "")

# Prelim line plots
par(mfrow = c(2, 2))

temp = select_stations[north_stations]
temp = gw_dta[temp,]
temp = temp - rowMeans(temp, na.rm = TRUE)
temp = temp / apply(temp, 1, sd, na.rm = TRUE)
for (ii in 1:nrow(temp)) {
  if (ii == 1) {
    plot(temp[ii,], type = "l", lwd = 0.5, xlab = "", ylab = "",
         main = "north",
         ylim = range(temp, na.rm = TRUE),
         xaxt = "n")
  } else {
    lines(temp[ii,], lwd = 0.5)
  }
}

temp = select_stations[west_stations]
temp = gw_dta[temp,]
temp = temp - rowMeans(temp, na.rm = TRUE)
temp = temp / apply(temp, 1, sd, na.rm = TRUE)
for (ii in 1:nrow(temp)) {
  if (ii == 1) {
    plot(temp[ii,], type = "l", lwd = 0.5, xlab = "", ylab = "",
         main = "west",
         ylim = range(temp, na.rm = TRUE),
         xaxt = "n")
  } else {
    lines(temp[ii,], lwd = 0.5)
  }
}

temp = select_stations[south_stations]
temp = gw_dta[temp,]
temp = temp - rowMeans(temp, na.rm = TRUE)
temp = temp / apply(temp, 1, sd, na.rm = TRUE)
for (ii in 1:nrow(temp)) {
  if (ii == 1) {
    plot(temp[ii,], type = "l", lwd = 0.5, xlab = "", ylab = "",
         main = "south",
         ylim = range(temp, na.rm = TRUE),
         xaxt = "n")
  } else {
    lines(temp[ii,], lwd = 0.5)
  }
}

temp = select_stations[east_stations]
temp = gw_dta[temp,]
temp = temp - rowMeans(temp, na.rm = TRUE)
temp = temp / apply(temp, 1, sd, na.rm = TRUE)
for (ii in 1:nrow(temp)) {
  if (ii == 1) {
    plot(temp[ii,], type = "l", lwd = 0.5, xlab = "", ylab = "",
         main = "east",
         ylim = range(temp, na.rm = TRUE),
         xaxt = "n")
  } else {
    lines(temp[ii,], lwd = 0.5)
  }
}

# All station plots

temp = select_stations
temp = gw_dta[temp,]
temp = temp - rowMeans(temp, na.rm = TRUE)
temp = temp / apply(temp, 1, sd, na.rm = TRUE)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
for (ii in 1:nrow(temp)) {
  if (ii == 1) {
    plot(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
         y = temp[ii,], type = "l", lwd = 0.5, xlab = "", ylab = "",
         main = "Groundwater level (normalized)",
         ylim = range(temp, na.rm = TRUE))
  } else {
    lines(x = seq(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
          y = temp[ii,], lwd = 0.5)
  }
}

# Missing patterns
GW = gw_dta[select_stations,]

par(mfrow = c(1, 1))
hist(rowMeans(is.na(GW)), breaks = 15, xlab = "Missing percentages", main = "")
summary(rowMeans(is.na(GW)))


missing_percentage <- colMeans(is.na(GW))
data <- data.frame(
  dates = seq.Date(from = as.Date("1992-10-01"), to = as.Date("2020-08-01"), by = "month"),
  missing_percentage = missing_percentage
)
ggplot(data, aes(x = dates, y = missing_percentage)) +
  geom_line(color = "blue") +
  scale_x_date(
    date_breaks = "2 years",  
    date_labels = "%Y"      
  ) +
  labs(x = "Year", y = "Missing Percentage", title = "") +
  theme(
    axis.title.x = element_text(size = 14),  # Increase x-axis label size
    axis.title.y = element_text(size = 14),  # Increase y-axis label size
    axis.text.x = element_text(size = 12),   # Increase x-axis tick text size
    axis.text.y = element_text(size = 12)    # Increase y-axis tick text size
  )



# Save data
write.table(GW, "./real_data/GW_select.csv", row.names = F, col.names = F, sep = ",")













