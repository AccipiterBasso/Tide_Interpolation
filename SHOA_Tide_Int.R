# --------------------------------------------------------------------------
# Author: Enzo Basso
# Affiliation: Bird Ecology Lab
# Email: ebassoq@gmail.com
# Date: 2025-12-08
# Description: Code to interpolate the tide based on the formula 
# provided by the Servicio Hidrográfico y Oceanográfico de la Armada (SHOA)
# --------------------------------------------------------------------------
# --- Versions used ---
# R version 4.5.0 (2025-04-11)
# zoo 1.8.14

# --- Data example ---
# Description: Data taken from the 2018 tide tables provided by SHOA
# Each row represents a tide event (high or low tide)
# Notes:
# 2 days of records (03-01-2018 and 04-01-2018)
# 4 tides per day (approx. 2 high and 2 low tides), consistent with a semidiurnal tide pattern
data <- data.frame(
  Date = c('03-01-2018', '03-01-2018', '03-01-2018', '03-01-2018', 
           '04-01-2018', '04-01-2018', '04-01-2018', '04-01-2018'),
  Time = c('02:14', '08:55', '14:09', '20:42', 
           '03:02', '09:43', '15:47', '21:32'),
  Tide = c(2.66, 0.06, 1.98, 0.44, 2.65, 0.08, 1.96, 0.48)
)
head(data)

# --- Library loading ---
library('zoo')
# --- Step 1: Calculate differences and decimal hours ---
diff_a = function(data) {
  data$Datetime = as.POSIXct(paste(data$Date, data$Time), format = '%d-%m-%Y %H:%M')
  diff_min = c(diff(as.numeric(data$Datetime)) / 60, NA)
  data$Diff_a = ifelse(is.na(diff_min), NA, sprintf('%02d:%02d', diff_min %/% 60, round(diff_min %% 60)))
  data$a = ifelse(is.na(diff_min), NA, diff_min / 60)
  data$c = c(abs(diff(data$Tide)), NA)
  return(data)
}

# --- Step 2: Interpolate per minute (linear interpolation) ---
int_pol = function(data) {
  sec = seq(min(data$Datetime), max(data$Datetime), by = '1 min')
  dat_int = data.frame(
    Datetime = sec,
    Date = format(sec, '%d-%m-%Y'),
    Time = format(sec, '%H:%M'),
    Tide = NA,
    Diff_a = NA,
    a = NA,
    c = NA
    )
  idx = match(data$Datetime, dat_int$Datetime)
  dat_int$Tide[idx] = data$Tide
  dat_int$Diff_a[idx] = data$Diff_a
  dat_int$a[idx] = data$a
  dat_int$c[idx] = data$c
  #.
  dat_int$a = na.locf(dat_int$a, na.rm = FALSE)
  dat_int$a = na.locf(dat_int$a, fromLast = TRUE)
  dat_int$c = na.locf(dat_int$c, na.rm = FALSE)
  dat_int$c = na.locf(dat_int$c, fromLast = TRUE)
  
  return(dat_int)
}

# --- Step 3: Half-step interpolation (halving between two known tide points) ---
int_hal = function(data) {
  Tide = data$Tide
  indna = which(!is.na(Tide))
  Tide_int = Tide
  for(i in 1:(length(indna) - 1)) {
    val_i = Tide[indna[i]]
    val_e = Tide[indna[i + 1]]
    indbet = (indna[i] + 1):(indna[i + 1] - 1)
    if(length(indbet) > 0) {
      mid = floor(length(indbet)/2)
      if(mid > 0) {
        Tide_int[indbet[1:mid]] = val_i
        Tide_int[indbet[(mid + 1):length(indbet)]] = val_e
      }
    }
  }
  data$Tide = Tide_int
  return(data)
}

# --- Step 4: Calculate time difference to nearest tide peak ---
diff_b = function(data_int, data) {
  data_int$Diff_b = NA
  data_int$b = NA
  n_datetime = as.POSIXct(paste(data$Date, data$Time), format = '%d-%m-%Y %H:%M')
  for(i in 1:nrow(data_int)) {
    time_int = data_int$Datetime[i]
    diff_sec = abs(as.numeric(n_datetime) - as.numeric(time_int))
    idx_n = which.min(diff_sec)
    diff_min = diff_sec[idx_n] / 60
    h = diff_min %/% 60
    min = round(diff_min %% 60)
    data_int$Diff_b[i] = sprintf('%02d:%02d', h, min)
    data_int$b[i] = diff_min / 60
  }
  return(data_int)
}

# --- Step 5: Calculate tidal correction using SHOA formula ---
corr = function(data) {
  data$cor = NA
  for(i in 1:nrow(data)) {
    if(!is.na(data$a[i]) && !is.na(data$c[i])) {
      a_c = data$a[i]
      c_c = data$c[i]
    } else if (i > 1 && !is.na(data$a[i-1]) && !is.na(data$c[i-1])) {
      a_c = data$a[i-1]
      c_c = data$c[i-1]
    } else {
      a_c = data$a[!is.na(data$a)][1]
      c_c = data$c[!is.na(data$c)][1]
    }
      b_c = data$b[i]
      if(!is.na(b_c)) {
        #. SHOA formula
        data$cor[i] = c_c / 2 * (1 - cos(180 * b_c / a_c * pi / 180))
      }
  }
  return(data)
}

# --- Step 6:  Interpolate tide with SHOA correction ---

int_tides = function(data) {
  data$Diff_tide =  NA
  uni_val = unique(data$Tide)
  less_than = uni_val[1] < uni_val[2]
  Tide = data$Tide
  Cor = data$cor
  c_oper = ifelse(less_than, 'sum', 'subtract')
  for(i in 1:nrow(data)) {
    if(i > 1 && Tide[i] != Tide[i-1]) {
      c_oper = ifelse(c_oper == 'sum', 'subtract', 'sum')
    }
    if(c_oper == 'sum') {
      data$Diff_tide[i] = Tide[i] + Cor[i]
    } else {
      data$Diff_tide[i] = Tide[i] - Cor[i]
    }
  }
  return(data)
}

# --- Step 7: Complete pipeline ---
com_pip = function(data) {
  data = diff_a(data)
  int_dat1 = int_pol(data)
  int_dat2 = int_hal(int_dat1)
  int_dat3 = diff_b(int_dat2, data)
  int_dat4 = corr(int_dat3)
  int_dat5 = int_tides(int_dat4)
  return(int_dat5)
}

# --- Run complete pipeline ---
output = com_pip(data)
head(output)

# --- Plot output ---
plot(output$Datetime, output$Diff_tide, type = 'b', 
     xlab = 'Time (HH:MM)', ylab = 'Height (m)', 
     main = '', col = '#CCCCCC', pch = 16, cex = 0.4, xaxt = 'n')
#. Points: High tide and low tide values obtained from tide tables
t_val = sapply(data$Tide, function(v) which.min(abs(output$Diff_tide - v)))
points(output$Datetime[t_val], data$Tide, col = c('#619cff', '#fb8072'), pch = 19)
#. Vertical line indicating the change of day (midnight)
hm = format(output$Datetime, '%H:%M')
mn = which(hm == '00:00')
abline(v = output$Datetime[mn], lty = 3, lwd = 2, col = '#669900')
#. High tide and low tide times obtained from tide tables
axis(1, at = output$Datetime[t_val], labels = data$Time)
#. Add formula annotation on top-left corner with lowercase delta 
formula_text <- expression(delta == frac(c, 2) * '[' * 1 - cos(180 * frac(b, a)) * ']')
text(x = min(output$Datetime + 3600), y = max(output$Diff_tide -0.2, na.rm=TRUE), 
     labels = formula_text, pos = 4, cex = 0.8, col = 'black')

# ---  Example: Using one month of tide data ---
#. NOTE: This function is designed to interpolate tide data by month
#. Using multiple months together may give incorrect estimates at month boundaries due to missing values
#. The data correspond to the tide tables for the port of Ancud in January 2018, obtained from SHOA
dt = read.csv('Ancud_Jan_Tides.csv', sep = ',', header = TRUE)

#. Run complete pipeline 
out = com_pip(dt)

#. Plot output
plot(out$Datetime, out$Diff_tide, type = 'b', 
     xlab = 'Time (HH:MM)', ylab = 'Height (m)', 
     main = '', col = '#CCCCCC', pch = 16, cex = 0.4, xaxt = 'n')
#. High tide and low tide times obtained from tide tables
t_val = sapply(dt$Tide, function(v) which.min(abs(out$Diff_tide - v)))
#. Vertical line indicating the change of day (midnight)
abline(v = out$Datetime[mn], lty = 3, lwd = 1, col = '#669900')
hm = format(out$Datetime, '%H:%M')
mn = which(hm == '00:00')
axis(1, at = out$Datetime[t_val], labels = dt$Time)









