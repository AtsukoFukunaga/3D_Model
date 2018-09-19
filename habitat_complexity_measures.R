library(raster)
library(rgeos)
library(ggplot2)


########## terrain fucntion returns slope, aspect, curvature #########
terrain_fun <- function(data, cell_size) {
  d0 <- as.matrix(data)
  if (ncol(d0) == 1 | nrow(d0) == 1) {
    terrain_list <- list(NA, NA, NA, NA, NA, NA)
    names(terrain_list) <- c("mean_slope", "mean_aspect", "circular_mean_aspect", 
                             "mean_curvature", "mean_profile_curvature", "mean_plan_curvature")
  } else {
    da <-  cbind(matrix(NA, nrow = nrow(d0), ncol = 1), 
                 rbind(matrix(NA, nrow = 1, ncol = ncol(d0) - 1), 
                       matrix(d0[1 : (nrow(d0) - 1), 1 : (ncol(d0) - 1)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1)))
    db <- rbind(matrix(NA, nrow = 1, ncol = ncol(d0)), 
                matrix(d0[1 : (nrow(d0) - 1), ], nrow = nrow(d0) - 1, ncol = ncol(d0)))
    dc <- cbind(rbind(matrix(NA, nrow = 1, ncol = ncol(d0) - 1), 
                      matrix(d0[1 : (nrow(d0) - 1), 2 : ncol(d0)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1)),
                matrix(NA, nrow = nrow(d0), ncol = 1))
    dd <- cbind(matrix(NA, nrow = nrow(d0), ncol = 1), 
                matrix(d0[, 1 : (ncol(d0) - 1)], nrow = nrow(d0), ncol = ncol(d0) -1))
    df <- cbind(matrix(d0[, 2 : ncol(d0)], nrow = nrow(d0), ncol = ncol(d0) - 1), 
                matrix(NA, nrow = nrow(d0), ncol = 1))
    dg <- cbind(matrix(NA, nrow = nrow(d0), ncol = 1), 
                rbind(matrix(d0[2 : nrow(d0), 1 : (ncol(d0) - 1)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1), 
                      matrix(NA, nrow = 1, ncol = ncol(d0) - 1)))
    dh <- rbind(matrix(d0[2 : nrow(d0), ], nrow = nrow(d0) - 1, ncol = ncol(d0)), 
                matrix(NA, nrow = 1, ncol = ncol(d0)))
    di <- cbind(rbind(matrix(d0[2 : nrow(d0), 2 : ncol(d0)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1), 
                      matrix(NA, nrow = 1, ncol = (ncol(d0) - 1))), 
                matrix(NA, nrow = nrow(d0), ncol = 1))
    
    ## slope
    
    x_rate <- ((dc + (2 * df) + di) - (da + (2 * dd) + dg)) / (8 * cell_size)
    y_rate <- ((dg + (2 * dh) + di) - (da + (2 * db) + dc)) / (8 * cell_size)
    
    slope_degrees <- atan(sqrt(x_rate ^ 2 + y_rate ^ 2)) * 57.29578
    mean_slope <- mean(slope_degrees, na.rm = TRUE)
    
    ## aspect
    
    x_rate_aspect <- ((dc + (2 * df) + di) - (da + (2 * dd) + dg)) / 8
    y_rate_aspect <- ((dg + (2 * dh) + di) - (da + (2 * db) + dc)) / 8
    
    aspect <- 57.29578 * atan2(y_rate_aspect, (-1 * x_rate_aspect))
    
    aspect_conversion <- function(x) {
      if (is.na(x)) {
        cell <- NA
      } else if (x < 90.0) {
        cell <- 90.0 - x
      } else {
        cell <- 450 - x
      }
      return(cell)
    }
    
    aspect_conv <- apply(aspect, c(1, 2), aspect_conversion)
    mean_aspect <- mean(aspect_conv, na.rm = TRUE)
    
    # convert aspect_conv to radians
    aspect_rad <- aspect_conv * (pi / 180)
    mean_sin <- mean(sin(aspect_rad), na.rm = TRUE)
    mean_cos <- mean(cos(aspect_rad), na.rm = TRUE)
    mean_aspect_degree <- 57.29578 * atan2(mean_sin, mean_cos)
    mean_aspect_2 <- mean_aspect_degree %% 360
    
    ## curvature
    
    coefD <- (((dd + df) / 2) - d0) / (cell_size ^ 2)
    coefE <- (((db + dh) / 2) - d0) / (cell_size ^ 2)
    coefF <- (dc + dg -da - di) / (4 * (cell_size ^ 2))
    coefG = (df - dd) / 2 * cell_size
    coefH = (db - dh) / 2 * cell_size
    
    curvature <- -2 * (coefD + coefE)
    mean_curvature <- mean(curvature, na.rm = TRUE)
    
    prof_curv <- 2 * ((coefD * coefG ^ 2 + coefE * coefH ^ 2 + coefF * coefG * coefH) / 
                         (coefG ^ 2 + coefH ^ 2))
    mean_prof_curv <- mean(prof_curv, na.rm = TRUE)
    
    plan_curv <- -2 * ((coefD * coefH ^ 2 + coefE * coefG ^ 2 - coefF * coefG * coefH) / 
                        (coefG ^ 2 + coefH ^ 2))
    mean_plan_curv <- mean(plan_curv, na.rm = TRUE)

    terrain_list <- list(mean_slope, mean_aspect, mean_aspect_2, mean_curvature, 
                         mean_prof_curv, mean_plan_curv)
    names(terrain_list) <- c("mean_slope", "mean_aspect", "circular_mean_aspect", 
                             "mean_curvature", "mean_profile_curvature", "mean_plan_curvature")
  }
  
  return(terrain_list)
  
}


########## process files #########

file_path <- "~/Habitat_complexity_test/DEMs/"
file_names <- c("LIS_4254")
file_suf <- "_1_1.tif"

files <- paste(file_path, file_names, file_suf, sep = "")
resolution <- 0.01  # DEM resolution = 1cm
return_resolution_factor <- 1  # resolution factor for return values, use 1, 2, 4, 8, 16, 32, 64. 1 for 1 cm, 2 for 2 cm etc.

hab_data <- data.frame(file_name = character(0), fd64 = numeric(0), fd128 = numeric(0),
                       planerS64 = numeric(0), planerS128 = numeric(0),
                       rugosity = numeric(0), max_h = numeric(0), min_h = numeric(0),
                       diff_h = numeric(0), sd_h = numeric(0),
                       mean_slope = numeric(0), mean_aspect = numeric(0), 
                       circular_mean_aspect = numeric(0), 
                       mean_curvature = numeric(0), 
                       mean_profile_curvature = numeric(0), 
                       mean_plan_curvature = numeric(0))

for (k in 1:length(files)) {
  ras <- raster(files[k])
  
  print(file_names[k])
  
  extent64 <- aggregate(ras, fac = 64, fun = mean, expand = FALSE, na.rm = FALSE)
  extent64_uniform <- extent64 > -Inf
  extent64_polygon <- rasterToPolygons(extent64_uniform, dissolve = TRUE)
  
  extent128 <- aggregate(ras, fac = 128, fun = mean, expand = FALSE, na.rm = FALSE)
  extent128_uniform <- extent128 > -Inf
  extent128_polygon <- rasterToPolygons(extent128_uniform, dissolve = TRUE)
  
  ras_clip_64 <- mask(crop(ras, extent(extent64_polygon)), extent64_polygon)
  ras_clip_128 <- mask(crop(ras, extent(extent128_polygon)), extent128_polygon)
  
  dat64 <- data.frame(fac = c(1, 2, 4, 8, 16, 32, 64), s_area = NA, area = NA, 
                      max_height = NA, min_height = NA, height_std = NA,
                      mean_slope = NA, mean_aspect = NA, circular_mean_aspect = NA, 
                      mean_curvature = NA, mean_profile_curvature = NA, mean_plan_curvature = NA)
  dat128 <- data.frame(fac = c(1, 2, 4, 8, 16, 32, 64, 128), s_area = NA, area = NA, 
                       max_height = NA, min_height = NA, height_std = NA,
                       mean_slope = NA, mean_aspect = NA, circular_mean_aspect = NA, 
                       mean_curvature = NA, mean_profile_curvature = NA, mean_plan_curvature = NA)
  
  dat <- list(dat64, dat128)
  
  
  for (j in 1:length(dat)) {
    
    if (j == 1) {
      for (i in 1:nrow(dat[[j]])) {
        
        print(dat[[j]]$fac[i])
        
        if (dat[[j]]$fac[i] == 1) {
          reef <- ras_clip_64
        } else {
          reef <- aggregate(ras_clip_64, fac = dat[[j]]$fac[i], fun = mean, expand = FALSE, na.rm = FALSE)
        }
        reef_g <- as(reef, "SpatialGridDataFrame")
        terrain_res <- terrain_fun(reef, resolution)
        dat[[j]]$s_area[i] <- surfaceArea(reef_g)
        dat[[j]]$area[i] <- sum(!is.na(reef_g@data)) * (dat[[j]]$fac[i] * 0.01) ^ 2
        dat[[j]]$max_height[i] <- max(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$min_height[i] <- min(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$height_std[i] <- sd(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$mean_slope[i] <- terrain_res$mean_slope
        dat[[j]]$mean_aspect[i] <- terrain_res$mean_aspect
        dat[[j]]$circular_mean_aspect[i] <- terrain_res$circular_mean_aspect
        dat[[j]]$mean_curvature[i] <- terrain_res$mean_curvature
        dat[[j]]$mean_profile_curvature[i] <- terrain_res$mean_profile_curvature
        dat[[j]]$mean_plan_curvature[i] <- terrain_res$mean_plan_curvature
      }
    } else if (j == 2) {
      for (i in 1:nrow(dat[[j]])) {
        
        print(dat[[j]]$fac[i])
        
        if (dat[[j]]$fac[i] == 1) {
          reef <- ras_clip_128
        } else {
          reef <- aggregate(ras_clip_128, fac = dat[[j]]$fac[i], fun = mean, expand = FALSE, na.rm = FALSE)
        }
        reef_g <- as(reef, "SpatialGridDataFrame")
        terrain_res <- terrain_fun(reef, resolution)
        dat[[j]]$s_area[i] <- surfaceArea(reef_g)
        dat[[j]]$area[i] <- sum(!is.na(reef_g@data)) * (dat[[j]]$fac[i] * 0.01) ^ 2
        dat[[j]]$max_height[i] <- max(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$min_height[i] <- min(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$height_std[i] <- sd(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$mean_slope[i] <- terrain_res$mean_slope
        dat[[j]]$mean_aspect[i] <- terrain_res$mean_aspect
        dat[[j]]$circular_mean_aspect[i] <- terrain_res$circular_mean_aspect
        dat[[j]]$mean_curvature[i] <- terrain_res$mean_curvature
        dat[[j]]$mean_profile_curvature[i] <- terrain_res$mean_profile_curvature
        dat[[j]]$mean_plan_curvature[i] <- terrain_res$mean_plan_curvature      }
    } else {
      break
    }
  }
  
  p64 <- ggplot(dat[[1]], aes(x = log(fac * 0.01), y = log(s_area))) +
    geom_point() +
    geom_smooth(method = "lm", se = F) +
    xlim(c(-5, 0)) + ylim(c(4.4, 6.5)) + 
    labs(x = expression(paste("log(", delta, ")")),
         y = expression(paste("logS(", delta, ")")))
  
  p128 <- ggplot(dat[[2]], aes(x = log(fac * 0.01), y = log(s_area))) +
    geom_point() +
    geom_smooth(method = "lm", se = F) + 
    xlim(c(-5, 1)) + ylim(c(4.4, 6.5)) + 
    labs(x = expression(paste("log(", delta, ")")),
         y = expression(paste("logS(", delta, ")"))) +
    theme_bw()
  
  # png("FD_plot", height = 1.2, width = 1.2, units = 'in', 
  #     type = "windows", res = 300)
  # p128
  # dev.off()
  
  # plot(p64)
  # plot(p128)
  
  d64 <- lm(log(s_area/area) ~ log(fac * 0.01), data = dat[[1]])
  slope64 <- coef(d64)[[2]]
  fd64 <- 2 - slope64
  
  d128 <- lm(log(s_area/area) ~ log(fac * 0.01), data = dat[[2]])
  slope128 <- coef(d128)[[2]]
  fd128 <- 2 - slope128
  
  fac <- which(dat[[1]]$fac == return_resolution_factor)
  
  rugosity64 <- dat[[1]]$s_area[fac]/dat[[1]]$area[fac]
  
  planerS64 <- dat[[1]]$area[fac]
  planerS128 <- dat[[2]]$area[fac]
  
  max_h <- dat[[1]]$max_height[fac]
  min_h <- dat[[1]]$min_height[fac]
  diff_h <- max_h - min_h
  sd_h <- dat[[1]]$height_std[fac]
  
  mean_slope <- dat[[1]]$mean_slope[fac]
  mean_aspect <- dat[[1]]$mean_aspect[fac]
  circular_mean_aspect <- dat[[1]]$circular_mean_aspect[fac]
  mean_curvature <- dat[[1]]$mean_curvature[fac]
  mean_profile_curvature <- dat[[1]]$mean_profile_curvature[fac]
  mean_plan_curvature <- dat[[1]]$mean_plan_curvature[fac]

  temp <- data.frame(file_name = file_names[k], fd64 = fd64, fd128 = fd128,
                     planerS64 = planerS64, planerS128 = planerS128,
                     rugosity = rugosity64, max_h = max_h, min_h = min_h,
                     diff_h = diff_h, sd_h = sd_h, 
                     mean_slope = mean_slope, mean_aspect = mean_aspect, 
                     circular_mean_aspect = circular_mean_aspect, 
                     mean_curvature = mean_curvature, 
                     mean_profile_curvature = mean_profile_curvature, 
                     mean_plan_curvature = mean_plan_curvature)
  
  hab_data <- rbind(hab_data, temp)
  
}

write.csv(hab_data, "habitat_complexity.csv", row.names = FALSE)
