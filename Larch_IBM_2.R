library(ggplot2)
library(patchwork)

# --- Simulate DEM: single peak from a bivariate normal ---
n  <- 500          # grid size
cs <- 30           # cell size (m)
H  <- 1000         # peak height (m)
sx <- 2000; sy <- 3000; rho <- 0.5   # std devs (m) and correlation

xs <- (seq_len(n) - n / 2) * cs
ys <- (seq_len(n) - n / 2) * cs
zx <- xs / sx
zy <- ys / sy

# rows = y (row 1 = north), cols = x
dem <- H * exp(-1 / (2 * (1 - rho^2)) *
                 (outer(rep(1, n), zx^2) - 2 * rho * outer(zy, zx) + outer(zy^2, rep(1, n))))

# --- Hillshade (Horn's method, as in ESRI/GDAL) ---
hillshade <- function(dem, cs, altitude = 45, azimuth = 180) {
  n <- nrow(dem); m <- ncol(dem)
  p <- rbind(dem[1, ], dem, dem[n, ])   # pad edges by replication
  p <- cbind(p[, 1], p, p[, m])
  r <- function(i, j) p[i:(i + n - 1), j:(j + m - 1)]
  a <- r(1,1); b <- r(1,2); c <- r(1,3)
  d <- r(2,1);              f <- r(2,3)
  g <- r(3,1); h <- r(3,2); i <- r(3,3)
  
  dzdx <- ((c + 2*f + i) - (a + 2*d + g)) / (8 * cs)
  dzdy <- ((g + 2*h + i) - (a + 2*b + c)) / (8 * cs)
  
  slope  <- atan(sqrt(dzdx^2 + dzdy^2))
  aspect <- atan2(dzdy, -dzdx)
  zenith <- (90 - altitude) * pi / 180
  az     <- ((360 - azimuth + 90) %% 360) * pi / 180
  
  hs <- 255 * (cos(zenith) * cos(slope) +
                 sin(zenith) * sin(slope) * cos(az - aspect))
  pmax(hs, 0)
}

hs <- hillshade(dem, cs)   # 500 x 500 matrix, values 0-255
dim(hs)




# --- Plot ---
image(t(hs[n:1, ]), col = gray.colors(256, start = 0, end = 1),
      useRaster = TRUE, axes = FALSE, asp = 1)

# Convert a matrix to long format (row 1 plotted at the top, as in the earlier plot)
mat_to_df <- function(m, value_name) {
  df <- data.frame(
    x = rep(seq_len(ncol(m)), each = nrow(m)),
    y = rep(rev(seq_len(nrow(m))), times = ncol(m)),
    v = as.vector(m)
  )
  names(df)[3] <- value_name
  df
}

dem_df <- mat_to_df(dem, "elevation")
hs_df  <- mat_to_df(hs,  "hillshade")

theme_map <- theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"))

p_dem <- ggplot(dem_df, aes(x, y, fill = elevation)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", name = "Elevation (m)") +
  coord_equal(expand = FALSE) +
  labs(title = "Simulated DEM") +
  theme_map

p_hs <- ggplot(hs_df, aes(x, y, fill = hillshade)) +
  geom_raster() +
  scale_fill_gradient(low = "black", high = "white", name = "Hillshade") +
  coord_equal(expand = FALSE) +
  labs(title = "Hillshade") +
  theme_map

p_dem + p_hs



# --- Sprinkle Larch seeds and see if they germinate ---
Pointnum=100

#set maximum hillshade for good habitat for larch germination
max_hs=150

#define outputs
output.X=matrix(nrow = Pointnum, ncol = 1)
output.Y=matrix(nrow = Pointnum, ncol = 1)

#set a counter
outcounter=0

#main for loop
for (i in 1:100000){
  #generate a random candidate point
  x.candidate=runif(1, min=1, max=ncol(hs))
  y.candidate=runif(1, min=1, max=nrow(hs))
  
  #calculate the probability of retaining the candidate point
  if(hs[y.candidate,x.candidate]>max_hs){
    prob=0.001
  }
  else
    prob=1-(hs[y.candidate,x.candidate]/max_hs)
  
  #roll the dice to see if you keep the candidate point
  testval=runif(1,min=0,max=1)  
  if (testval<prob){
    outcounter=outcounter+1
    keep=1
    output.X[outcounter]=x.candidate
    output.Y[outcounter]=y.candidate
  }
  
  #if you've reached your target number of points, break from the for loop
  if(outcounter==Pointnum){
    break
  }
}

# Points data frame 
pts <- data.frame(x = output.X, y = nrow(hs) + 1 - output.Y)

p_hs2 <- p_hs +
  geom_point(data = pts, aes(x = x, y = y),
             shape = 1, color = "red", size = 2, stroke = 0.8,
             inherit.aes = FALSE)

p_dem + p_hs2




#--- Next Generation ------

#Each larch gets to produce one surviving seed to seed the next generation
#Their dispersal width is set by disp_width
disp_width=10

#main for loop
new.X <- numeric(Pointnum)
new.Y <- numeric(Pointnum)

for (j in 1:Pointnum){
  for (i in 1:100000){
    #generate a random candidate point near parent j
    x.candidate=runif(1, min=output.X[j]-disp_width, max=output.X[j]+disp_width)
    y.candidate=runif(1, min=output.Y[j]-disp_width, max=output.Y[j]+disp_width)
    
    #keep the candidate inside the grid
    x.candidate=min(max(x.candidate, 1), ncol(hs))
    y.candidate=min(max(y.candidate, 1), nrow(hs))
    
    #calculate the probability of retaining the candidate point
    if(hs[y.candidate,x.candidate]>max_hs){
      prob=0.001
    } else
      prob=1-(hs[y.candidate,x.candidate]/max_hs)
      #prob=0.001
    
    #roll the dice to see if you keep the candidate point
    if (runif(1,min=0,max=1)<prob){
      new.X[j]=x.candidate
      new.Y[j]=y.candidate
      break
    }
  }
}

# Plot the new generation (yellow) on top of the original points (red)
pts_new <- data.frame(x = new.X, y = nrow(hs) + 1 - new.Y)

p_hs3 <- p_hs2 +   # p_hs2 from line 134 already includes the red points
  geom_point(data = pts_new, aes(x = x, y = y),
             shape = 1, color = "yellow", size = 2, stroke = 0.8,
             inherit.aes = FALSE)

p_dem + p_hs3
