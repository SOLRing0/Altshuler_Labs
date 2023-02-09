
################################### Packages ###################################
## Specify the packages of interest
packages <-
  c("ggplot2",
    "tidyverse",
    "tidyr",
    "dplyr",
    "tibble",
    "readr",
    "readxl",
    "ggthemes",
    "magick",
    "grid",
    "cowplot",
    "metR",
    "viridis",
    "R.matlab",
    "gaussplotR")

## Load or install & load all packages
package.check <- lapply(
  packages,
  FUN = function(x)
  {
    if (!require(x, character.only = TRUE))
    {
      install.packages(x, dependencies = TRUE,
                       repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)

#### import truncated data ####
lm_steady_dir <-
  "./data/spatemp/LM_ER_steadystate_200bins_bin100_to_bin200/"

lm_steady_files <-
  list.files(lm_steady_dir, recursive = TRUE)

lm_steady_filepaths <-
  paste0(
    lm_steady_dir,
    lm_steady_files
  )

lm_steady_basenames <-
  basename(lm_steady_files) %>%
  str_remove("_steadystate.csv")

## extract the relevant stuff
lm_steady_dat <-
  map(lm_steady_filepaths, read_csv, col_names = TRUE) 
names(lm_steady_dat) <- lm_steady_basenames


## cell1_norms
lm_steady_cell1_norms <-
  lm_steady_dat %>%
  map(select, logSF, logTF, ends_with("cell1_norm")) %>%
  map(drop_na)
vec <- NULL
for (i in 1:length(lm_steady_cell1_norms)) {
  vec[i] <- ncol(lm_steady_cell1_norms[[i]])
}
lm_steady_cell1_norms <-lm_steady_cell1_norms[vec > 2]
names(lm_steady_cell1_norms) <-
  paste0(names(lm_steady_cell1_norms),
         "-cell1")

## cell2_norms
lm_steady_cell2_norms <-
  lm_steady_dat %>%
  map(select, logSF, logTF, ends_with("cell2_norm")) %>%
  map(drop_na)
vec <- NULL
for (i in 1:length(lm_steady_cell2_norms)) {
  vec[i] <- ncol(lm_steady_cell2_norms[[i]])
}
lm_steady_cell2_norms <-lm_steady_cell2_norms[vec > 2]
names(lm_steady_cell2_norms) <-
  paste0(names(lm_steady_cell2_norms),
         "-cell2")

## cell3_norms
lm_steady_cell3_norms <-
  lm_steady_dat %>%
  map(select, logSF, logTF, ends_with("cell3_norm")) 
vec <- NULL
for (i in 1:length(lm_steady_cell3_norms)) {
  vec[i] <- ncol(lm_steady_cell3_norms[[i]])
}
lm_steady_cell3_norms <-lm_steady_cell3_norms[vec > 2]
names(lm_steady_cell3_norms) <-
  paste0(names(lm_steady_cell3_norms),
         "-cell3")

## cell4_norms
lm_steady_cell4_norms <-
  lm_steady_dat %>%
  map(select, logSF, logTF, ends_with("cell4_norm"))
vec <- NULL
for (i in 1:length(lm_steady_cell4_norms)) {
  vec[i] <- ncol(lm_steady_cell4_norms[[i]])
}
lm_steady_cell4_norms <-lm_steady_cell4_norms[vec > 2]
names(lm_steady_cell4_norms) <-
  paste0(names(lm_steady_cell4_norms),
         "-cell4")

#### get all the data ####
lm_steady_all_cell_norms <-
  c(lm_steady_cell1_norms,
    lm_steady_cell2_norms,
    lm_steady_cell3_norms,
    lm_steady_cell4_norms)  %>%
  .[order(names(.))]

#### specify cell of interest ####
names(lm_steady_all_cell_norms)
#i = 2

## extract data
dat <- 
  lm_steady_all_cell_norms[[i]] %>%
  drop_na() %>%
  as.data.frame()
colnames(dat) <-
  c("X_values", "Y_values", "response")

## contour plot
dat %>%
  ggplot(aes(x = X_values,
             y = Y_values,
             z = response)) +
  metR::geom_contour_fill(aes(fill = ..level..),
                          color = "black",
                          size = 0.03,
                          bins = 15) +
  xlab("log2 of SF (cpd)") +
  ylab("log2 of TF (Hz)") +
  ggtitle(
    paste0(
      names(lm_steady_all_cell_norms)[i],
      " contour")
  ) +
  coord_fixed(expand = FALSE) +
  theme_classic() +
  theme(plot.title = element_text(size = 8))


datj <- 
  dat %>%
  filter(X_values >=  -5 - 0.1) %>%
  filter(X_values <=  0 + 0.1) %>%
  filter(Y_values >=  -1 - 0.1) %>%
  filter(Y_values <=   4 + 0.1)

## truncated contour plot
datj %>%
  ggplot(aes(x = X_values,
             y = Y_values,
             z = response)) +
  metR::geom_contour_fill(aes(fill = ..level..),
                          color = "black",
                          size = 0.03,
                          bins = 15) +
  xlab("log2 of SF (cpd)") +
  ylab("log2 of TF (Hz)") +
  ggtitle(
    paste0(
      names(lm_steady_all_cell_norms)[i],
      " contour")
  ) +
  coord_fixed(expand = FALSE) +
  theme_classic() +
  theme(plot.title = element_text(size = 8))

uncon_mod <-
  gaussplotR::fit_gaussian_2D(
    data = datj,
    method = "elliptical_log",
    constrain_amplitude = TRUE,
    constrain_orientation = "unconstrained", 
    maxiter = 10000,
    minFactor = 0.000000488281
  )

ind_mod <-
  gaussplotR::fit_gaussian_2D(
    data = datj,
    method = "elliptical_log",
    constrain_amplitude = TRUE,
    constrain_orientation = -1, 
    maxiter = 10000,
    minFactor = 0.000000488281
  )

spd_mod <-
  gaussplotR::fit_gaussian_2D(
    data = datj,
    method = "elliptical_log",
    constrain_amplitude = TRUE,
    constrain_orientation = 0, 
    maxiter = 10000,
    minFactor = 0.000000488281
  )


