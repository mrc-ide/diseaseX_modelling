# -------------------------------------------------------------------------
# Generate spatial coordinates for infection offspring
# -------------------------------------------------------------------------
# parent_x, parent_y : numeric(1)  – coordinates of the parent
# n_offspring        : integer(1)  – how many offspring to place
# spatial_kernel     : function(n) – returns n random *distances*  (1-D PDF)
# kernel_is_radial   : logical(1)  – set TRUE if spatial_kernel is *already*
#                                    the correct 2-D radial sampler (e.g. rgamma(n, 2, σ))
# oversample_factor  : integer(1)  – only used when kernel_is_radial = FALSE.
#                                    Draw m = n * oversample_factor candidates and
#                                    resample with weights ∝ r.  A value of 5–10
#                                    is usually plenty.
# -------------------------------------------------------------------------
spatial_calc <- function(parent_x,
                         parent_y,
                         n_offspring,
                         spatial_kernel,
                         kernel_is_radial = FALSE,
                         oversample_factor = 10L) {
  
  stopifnot(n_offspring  > 0,
            oversample_factor >= 1,
            is.function(spatial_kernel))
  
  ## ----------------------------------------------------------------------
  ## 1. Draw radial distances --------------------------------------------
  ## ----------------------------------------------------------------------
  if (kernel_is_radial) {
    # user supplies a generator that already samples the 2-D radial law
    distance <- spatial_kernel(n_offspring)
    
  } else {
    # user supplies a *1-D* generator  f(r)  → need to up-weight large r
    m          <- n_offspring * oversample_factor
    candidates <- spatial_kernel(m)
    
    if (any(candidates < 0))
      stop("spatial_kernel produced negative distances; check the function.")
    
    # importance resampling: weights ∝ r  → target PDF ∝ r * f(r)
    w          <- candidates
    distance   <- sample(candidates, n_offspring,
                         replace = TRUE,
                         prob     = w)
  }
  
  ## ----------------------------------------------------------------------
  ## 2. Draw directions and compute coordinates (vectorised) --------------
  ## ----------------------------------------------------------------------
  theta <- runif(n_offspring, 0, 2 * pi)
  x     <- parent_x + distance * cos(theta)
  y     <- parent_y + distance * sin(theta)
  
  ## ----------------------------------------------------------------------
  ## 3. Return tidy data frame -------------------------------------------
  ## ----------------------------------------------------------------------
  data.frame(distance         = distance,
             x_coordinate     = x,
             y_coordinate     = y,
             overall_distance = sqrt(x^2 + y^2))  # distance from (0,0); tweak if needed
}

## 1. Uniform *line* density in [0, R]  (needs correction)
R                <- 1
exp_surface_rad  <- function(n) rexp(n, R)

coords1 <- spatial_calc(0, 0, 1e4, exp_surface_rad)
hist(coords1$overall_distance, breaks = 20)

## 2. Exponential surface kernel with scale σ  (already radial)
sigma            <- 1
exp_surface_rad  <- function(n) rgamma(n, 2, sigma)   # shape 2  ⇒ radial law

coords2 <- spatial_calc(0, 0, 5e4, exp_surface_rad,
                        kernel_is_radial = TRUE)

hist(coords1$overall_distance, breaks = 100)
hist(coords2$overall_distance, breaks = 100)


