#### data to superimpose all PSTH ####
## both species

superimp_dir <- "./data/spatemp/LM_ER_fullsweep_200bins_bin001_to_bin200/"
superimp_bins <- 200

all_superimp_files <-
  list.files(superimp_dir, recursive = TRUE)

allcell_psth_superimp <-
  paste0(
    superimp_dir,
    all_superimp_files[str_which(all_superimp_files,"-ALL_CELLS-SUBSET_PSTH_MEAN.csv")]
  )

allcell_spont_superimp <-
  paste0(
    superimp_dir,
    all_superimp_files[str_which(all_superimp_files,"-spontaneous_rate.csv")]
  )

allcell_ER_superimp <-
  paste0(
    superimp_dir,
    all_superimp_files[str_which(all_superimp_files,"-STER_data.csv")]
  )

allcell_basenames_superimp <-
  basename(allcell_psth_superimp) %>%
  str_remove("-ST-ALL_CELLS-SUBSET_PSTH_MEAN.csv")

## extract the relevant stuff
psth_superimp_dat <-
  map(allcell_psth_superimp, read_csv, col_names = FALSE) 
names(psth_superimp_dat) <- allcell_basenames_superimp

spont_superimp_dat <-
  map(allcell_spont_superimp, read_csv, col_names = FALSE) 
names(spont_superimp_dat) <- allcell_basenames_superimp

ER_superimp_dat <-
  map(allcell_ER_superimp, read_csv, col_names = FALSE) 
names(ER_superimp_dat) <- allcell_basenames_superimp


### Zebra finches only, LM ####

all_zf_files <- 
  all_superimp_files[grep("^STC", all_superimp_files)]


allcell_psth_superimp <-
  paste0(
    superimp_dir,
    all_zf_files[str_which(all_zf_files,"-ALL_CELLS-SUBSET_PSTH_MEAN.csv")]
  )

allcell_spont_superimp <-
  paste0(
    superimp_dir,
    all_zf_files[str_which(all_zf_files,"-spontaneous_rate.csv")]
  )

allcell_ER_superimp <-
  paste0(
    superimp_dir,
    all_zf_files[str_which(all_zf_files,"-STER_data.csv")]
  )

allcell_basenames_superimp <-
  basename(allcell_psth_superimp) %>%
  str_remove("-ST-ALL_CELLS-SUBSET_PSTH_MEAN.csv")

## extract the relevant stuff
psth_superimp_dat <-
  map(allcell_psth_superimp, read_csv, col_names = FALSE) 
names(psth_superimp_dat) <- allcell_basenames_superimp

spont_superimp_dat <-
  map(allcell_spont_superimp, read_csv, col_names = FALSE) 
names(spont_superimp_dat) <- allcell_basenames_superimp

ER_superimp_dat <-
  map(allcell_ER_superimp, read_csv, col_names = FALSE) 
names(ER_superimp_dat) <- allcell_basenames_superimp


## figure out cell counts
superimp_cell_counts <- NULL
for (i in 1:length(psth_superimp_dat)) {
  superimp_cell_counts[[i]] <-
    psth_superimp_dat[[i]] %>%
    select(-c(X1, X2, X3)) %>%
    ncol()
}
names(superimp_cell_counts) <- allcell_basenames_superimp
# make a df
superimp_cell_countdf <- 
  do.call(rbind, superimp_cell_counts) %>%
  as.data.frame() %>%
  rownames_to_column(var = "basename") %>%
  transmute(
    basename = basename,
    cell_count = V1/200
  )

## set up name vectors
superimp_cell1_names <-
  paste0(allcell_basenames_superimp, "-cell1")

superimp_cell2_names <-
  superimp_cell_countdf %>%
  filter(cell_count > 1) %>%
  select(basename) %>%
  as_vector() %>%
  paste0("-cell2")
superimp_cell2_basenames <-
  superimp_cell_countdf %>%
  filter(cell_count > 1) %>%
  select(basename) %>%
  as_vector()
names(superimp_cell2_basenames) <- NULL

## restrict additional data sets
cell2_select <- names(psth_superimp_dat) %in% superimp_cell2_basenames
psth_superimp_cell2_file <- psth_superimp_dat[cell2_select]
spont_superimp_cell2_file <- spont_superimp_dat[cell2_select]
# ER_superimp_cell2_file <- ER_superimp_dat[cell2_select]

## FIRST CELL
spont_superimp_cell1_dat <- NULL
for (i in 1:length(spont_superimp_dat)) {
  spont_superimp_cell1_dat[[i]] <-
    spont_superimp_dat[[i]][1, 1] %>%
    rename(spont = X1) %>%
    as_vector()
  colnames(spont_superimp_cell1_dat[[i]]) <- NULL
}
names(spont_superimp_cell1_dat) <- superimp_cell1_names

psth_superimp_cell1_dat <- NULL
for (i in 1:length(psth_superimp_dat)) {
  psth_superimp_cell1_dat[[i]] <-
    psth_superimp_dat[[i]] %>%
    filter(!X3 == 0.000668) %>%
    #select(4:103) %>%
    select(4:203) %>%
    as.matrix()
  ## subtract spont rate
  psth_superimp_cell1_dat[[i]] <-
    (psth_superimp_cell1_dat[[i]] - spont_superimp_cell1_dat[[i]])
  ## normalize
  psth_superimp_cell1_dat[[i]] <-
    psth_superimp_cell1_dat[[i]]/max(psth_superimp_cell1_dat[[i]])
  colnames(psth_superimp_cell1_dat[[i]]) <- NULL
}
names(psth_superimp_cell1_dat) <- superimp_cell1_names


## addition
psth_superimp_cell1_sumdat <-
  Reduce('+', psth_superimp_cell1_dat)
## averaging
psth_superimp_cell1_array <- simplify2array(psth_superimp_cell1_dat)
psth_superimp_cell1_avgdat <- apply(psth_superimp_cell1_array, c(1,2), mean)

## SECOND CELL
spont_superimp_cell2_dat <- NULL
for (i in 1:length(spont_superimp_cell2_file)) {
  spont_superimp_cell2_dat[[i]] <-
    spont_superimp_cell2_file[[i]][2, 1] %>%
    rename(spont = X1) %>%
    as_vector()
  colnames(spont_superimp_cell2_dat[[i]]) <- NULL
}
names(spont_superimp_cell2_dat) <- superimp_cell2_names
psth_superimp_cell2_dat <- NULL
for (i in 1:length(psth_superimp_cell2_file)) {
  psth_superimp_cell2_dat[[i]] <-
    psth_superimp_cell2_file[[i]] %>%
    filter(!X3 == 0.000668) %>%
    #select(104:203) %>%
    select(204:403) %>%
    as.matrix()
  ## subtract spont rate
  psth_superimp_cell2_dat[[i]] <-
    (psth_superimp_cell2_dat[[i]] - spont_superimp_cell2_dat[[i]])
  ## normalize
  psth_superimp_cell2_dat[[i]] <-
    psth_superimp_cell2_dat[[i]]/max(psth_superimp_cell2_dat[[i]])
  colnames(psth_superimp_cell2_dat[[i]]) <- NULL
}
names(psth_superimp_cell2_dat) <- superimp_cell2_names


psth_superimp_cell2_sumdat <-
  Reduce('+', psth_superimp_cell2_dat)
## averaging
psth_superimp_cell2_array <- simplify2array(psth_superimp_cell2_dat)
psth_superimp_cell2_avgdat <- apply(psth_superimp_cell2_array, c(1,2), mean)

## all together now
allsupercells_sumdat_lists <-
  c(psth_superimp_cell1_dat,
    psth_superimp_cell2_dat)
## averaging
allsupercells_array <- simplify2array(allsupercells_sumdat_lists)
allsupercells_array_avgdat_zf <- apply(allsupercells_array, c(1,2), mean)

psth_superimp_allsupercells_sumdat <-
  Reduce('+', allsupercells_sumdat_lists)


### ZF ALL TOGETHER plot ####
allsupercells_psth_letters <- seq(1:superimp_bins)
allsupercells_max_zf <- max(allsupercells_array_avgdat_zf)
allsupercells_plotz <- NULL
allsupercells_datlist <- NULL
for (i in 1:nrow(allsupercells_array_avgdat_zf)) {
  allsupercells_datlist[[i]] <-
    tibble(lets = allsupercells_psth_letters,
           psth_dat = allsupercells_array_avgdat_zf[i,])
  allsupercells_plotz[[i]] <-
    ggplot(allsupercells_datlist[[i]], aes(x = lets, y = psth_dat)) +
    geom_col(fill = col_zb,
             color = col_zb,
             size = 0#0.07
    ) +
    scale_y_continuous(breaks = seq(0, 0.30#allsupercells_max_zf
                                    )) +
    coord_cartesian(ylim = c(0, 0.30#allsupercells_max_zf
                             ),
                    expand = FALSE) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(colour = "black", size = 0.1),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.ticks = element_line(colour = "black", size = 0.1),
      axis.ticks.x = element_blank(),
      axis.ticks.length = unit(0.05, "cm")
    )
}

allzfcells_psthplot <-
  plot_grid(plotlist = rev(allsupercells_plotz),
            nrow = 6)

### Hummers only, LM ####

all_hb_files <- 
  all_superimp_files[grep("^CALSTC", all_superimp_files)]


allcell_psth_superimp <-
  paste0(
    superimp_dir,
    all_hb_files[str_which(all_hb_files,"-ALL_CELLS-SUBSET_PSTH_MEAN.csv")]
  )

allcell_spont_superimp <-
  paste0(
    superimp_dir,
    all_hb_files[str_which(all_hb_files,"-spontaneous_rate.csv")]
  )

allcell_ER_superimp <-
  paste0(
    superimp_dir,
    all_hb_files[str_which(all_hb_files,"-STER_data.csv")]
  )

allcell_basenames_superimp <-
  basename(allcell_psth_superimp) %>%
  str_remove("-ST-ALL_CELLS-SUBSET_PSTH_MEAN.csv")

## extract the relevant stuff
psth_superimp_dat <-
  map(allcell_psth_superimp, read_csv, col_names = FALSE) 
names(psth_superimp_dat) <- allcell_basenames_superimp

spont_superimp_dat <-
  map(allcell_spont_superimp, read_csv, col_names = FALSE) 
names(spont_superimp_dat) <- allcell_basenames_superimp

ER_superimp_dat <-
  map(allcell_ER_superimp, read_csv, col_names = FALSE) 
names(ER_superimp_dat) <- allcell_basenames_superimp


## figure out cell counts
superimp_cell_counts <- NULL
for (i in 1:length(psth_superimp_dat)) {
  superimp_cell_counts[[i]] <-
    psth_superimp_dat[[i]] %>%
    select(-c(X1, X2, X3)) %>%
    ncol()
}
names(superimp_cell_counts) <- allcell_basenames_superimp
# make a df
superimp_cell_countdf <- 
  do.call(rbind, superimp_cell_counts) %>%
  as.data.frame() %>%
  rownames_to_column(var = "basename") %>%
  transmute(
    basename = basename,
    cell_count = V1/200
  )

## set up name vectors
superimp_cell1_names <-
  paste0(allcell_basenames_superimp, "-cell1")

superimp_cell2_names <-
  superimp_cell_countdf %>%
  filter(cell_count > 1) %>%
  select(basename) %>%
  as_vector() %>%
  paste0("-cell2")
superimp_cell2_basenames <-
  superimp_cell_countdf %>%
  filter(cell_count > 1) %>%
  select(basename) %>%
  as_vector()
names(superimp_cell2_basenames) <- NULL

superimp_cell3_names <-
  superimp_cell_countdf %>%
  filter(cell_count > 2) %>%
  select(basename) %>%
  as_vector() %>%
  paste0("-cell3")
superimp_cell3_basenames <-
  superimp_cell_countdf %>%
  filter(cell_count > 2) %>%
  select(basename) %>%
  as_vector()

superimp_cell4_names <-
  superimp_cell_countdf %>%
  filter(cell_count > 3) %>%
  select(basename) %>%
  as_vector() %>%
  paste0("-cell4")
superimp_cell4_basenames <-
  superimp_cell_countdf %>%
  filter(cell_count > 3) %>%
  select(basename) %>%
  as_vector()

## restrict additional data sets
cell2_select <- names(psth_superimp_dat) %in% superimp_cell2_basenames
psth_superimp_cell2_file <- psth_superimp_dat[cell2_select]
spont_superimp_cell2_file <- spont_superimp_dat[cell2_select]
# ER_superimp_cell2_file <- ER_superimp_dat[cell2_select]

cell3_select <- names(psth_superimp_dat) %in% superimp_cell3_basenames
psth_superimp_cell3_file <- psth_superimp_dat[cell3_select]
spont_superimp_cell3_file <- spont_superimp_dat[cell3_select]
# ER_superimp_cell3_file <- ER_superimp_dat[cell3_select]

cell4_select <- names(psth_superimp_dat) %in% superimp_cell4_basenames
psth_superimp_cell4_file <- psth_superimp_dat[cell4_select]
spont_superimp_cell4_file <- spont_superimp_dat[cell4_select]
# ER_superimp_cell4_file <- ER_superimp_dat[cell4_select]

## FIRST CELL
spont_superimp_cell1_dat <- NULL
for (i in 1:length(spont_superimp_dat)) {
  spont_superimp_cell1_dat[[i]] <-
    spont_superimp_dat[[i]][1, 1] %>%
    rename(spont = X1) %>%
    as_vector()
  colnames(spont_superimp_cell1_dat[[i]]) <- NULL
}
names(spont_superimp_cell1_dat) <- superimp_cell1_names

psth_superimp_cell1_dat <- NULL
for (i in 1:length(psth_superimp_dat)) {
  psth_superimp_cell1_dat[[i]] <-
    psth_superimp_dat[[i]] %>%
    filter(!X3 == 0.000668) %>%
    #select(4:103) %>%
    select(4:203) %>%
    as.matrix()
  ## subtract spont rate
  psth_superimp_cell1_dat[[i]] <-
    (psth_superimp_cell1_dat[[i]] - spont_superimp_cell1_dat[[i]])
  psth_superimp_cell1_dat[[i]] <-
    psth_superimp_cell1_dat[[i]]/max(psth_superimp_cell1_dat[[i]])
  colnames(psth_superimp_cell1_dat[[i]]) <- NULL
}
names(psth_superimp_cell1_dat) <- superimp_cell1_names

## addition
psth_superimp_cell1_sumdat <-
  Reduce('+', psth_superimp_cell1_dat)
## averaging
psth_superimp_cell1_array <- simplify2array(psth_superimp_cell1_dat)
psth_superimp_cell1_avgdat <- apply(psth_superimp_cell1_array, c(1,2), mean)

## SECOND CELL
spont_superimp_cell2_dat <- NULL
for (i in 1:length(spont_superimp_cell2_file)) {
  spont_superimp_cell2_dat[[i]] <-
    spont_superimp_cell2_file[[i]][2, 1] %>%
    rename(spont = X1) %>%
    as_vector()
  colnames(spont_superimp_cell2_dat[[i]]) <- NULL
}
names(spont_superimp_cell2_dat) <- superimp_cell2_names

psth_superimp_cell2_dat <- NULL
for (i in 1:length(psth_superimp_cell2_file)) {
  psth_superimp_cell2_dat[[i]] <-
    psth_superimp_cell2_file[[i]] %>%
    filter(!X3 == 0.000668) %>%
    #select(104:203) %>%
    select(204:403) %>%
    as.matrix()
  ## subtract spont rate
  psth_superimp_cell2_dat[[i]] <-
    (psth_superimp_cell2_dat[[i]] - spont_superimp_cell2_dat[[i]])
  psth_superimp_cell2_dat[[i]] <-
    psth_superimp_cell2_dat[[i]]/max(psth_superimp_cell2_dat[[i]])
  colnames(psth_superimp_cell2_dat[[i]]) <- NULL
}
names(psth_superimp_cell2_dat) <- superimp_cell2_names

psth_superimp_cell2_sumdat <-
  Reduce('+', psth_superimp_cell2_dat)
## averaging
psth_superimp_cell2_array <- simplify2array(psth_superimp_cell2_dat)
psth_superimp_cell2_avgdat <- apply(psth_superimp_cell2_array, c(1,2), mean)


## THIRD CELL
spont_superimp_cell3_dat <- NULL
for (i in 1:length(spont_superimp_cell3_file)) {
  spont_superimp_cell3_dat[[i]] <-
    spont_superimp_cell3_file[[i]][3, 1] %>%
    rename(spont = X1) %>%
    as_vector()
  colnames(spont_superimp_cell3_dat[[i]]) <- NULL
}
names(spont_superimp_cell3_dat) <- superimp_cell3_names

psth_superimp_cell3_dat <- NULL
for (i in 1:length(psth_superimp_cell3_file)) {
  psth_superimp_cell3_dat[[i]] <-
    psth_superimp_cell3_file[[i]] %>%
    filter(!X3 == 0.000668) %>%
    #select(204:303) %>%
    select(404:603) %>%
    as.matrix()
  ## subtract spont rate
  psth_superimp_cell3_dat[[i]] <-
    (psth_superimp_cell3_dat[[i]] - spont_superimp_cell3_dat[[i]])
  psth_superimp_cell3_dat[[i]] <-
    psth_superimp_cell3_dat[[i]]/max(psth_superimp_cell3_dat[[i]])
  colnames(psth_superimp_cell3_dat[[i]]) <- NULL
}
names(psth_superimp_cell3_dat) <- superimp_cell3_names

psth_superimp_cell3_sumdat <-
  Reduce('+', psth_superimp_cell3_dat)
## averaging
psth_superimp_cell3_array <- simplify2array(psth_superimp_cell3_dat)
psth_superimp_cell3_avgdat <- apply(psth_superimp_cell3_array, c(1,2), mean)

## FOURTH CELL
spont_superimp_cell4_dat <- NULL
for (i in 1:length(spont_superimp_cell4_file)) {
  spont_superimp_cell4_dat[[i]] <-
    spont_superimp_cell4_file[[i]][4, 1] %>%
    rename(spont = X1) %>%
    as_vector()
  colnames(spont_superimp_cell4_dat[[i]]) <- NULL
}
names(spont_superimp_cell4_dat) <- superimp_cell4_names
psth_superimp_cell4_dat <- NULL
for (i in 1:length(psth_superimp_cell4_file)) {
  psth_superimp_cell4_dat[[i]] <-
    psth_superimp_cell4_file[[i]] %>%
    filter(!X3 == 0.000668) %>%
    #select(304:403) %>%
    select(604:803) %>%
    as.matrix()
  ## subtract spont rate
  psth_superimp_cell4_dat[[i]] <-
    (psth_superimp_cell4_dat[[i]] - spont_superimp_cell4_dat[[i]])
  psth_superimp_cell4_dat[[i]] <-
    psth_superimp_cell4_dat[[i]]/max(psth_superimp_cell4_dat[[i]])
  colnames(psth_superimp_cell4_dat[[i]]) <- NULL
}
names(psth_superimp_cell4_dat) <- superimp_cell4_names

psth_superimp_cell4_sumdat <-
  Reduce('+', psth_superimp_cell4_dat)
## averaging
psth_superimp_cell4_array <- simplify2array(psth_superimp_cell4_dat)
psth_superimp_cell4_avgdat <- apply(psth_superimp_cell4_array, c(1,2), mean)

allsupercells_sumdat_lists <-
  c(psth_superimp_cell1_dat,
    psth_superimp_cell2_dat, 
    psth_superimp_cell3_dat,
    psth_superimp_cell4_dat)
## averaging
allsupercells_array <- simplify2array(allsupercells_sumdat_lists)
allsupercells_array_avgdat_hb <- apply(allsupercells_array, c(1,2), mean)

psth_superimp_allsupercells_sumdat <-
  Reduce('+', allsupercells_sumdat_lists)


#### HB ALL TOGETHER plot ####
allsupercells_psth_letters <- seq(1:superimp_bins)
allsupercells_max_hb <- max(allsupercells_array_avgdat_hb)
allsupercells_plotz <- NULL
allsupercells_datlist <- NULL
for (i in 1:nrow(allsupercells_array_avgdat_hb)) {
  allsupercells_datlist[[i]] <-
    tibble(lets = allsupercells_psth_letters,
           psth_dat = allsupercells_array_avgdat_hb[i,])
  allsupercells_plotz[[i]] <-
    ggplot(allsupercells_datlist[[i]], aes(x = lets, y = psth_dat)) +
    geom_col(fill = col_hb,
             color = col_hb,
             size = 0#0.07
    ) +
    scale_y_continuous(breaks = seq(0, 0.30#allsupercells_max_hb
    )) +
    coord_cartesian(ylim = c(0, 0.30#allsupercells_max_hb
    ),
    expand = FALSE) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(colour = "black", size = 0.1),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.ticks = element_line(colour = "black", size = 0.1),
      axis.ticks.x = element_blank(),
      axis.ticks.length = unit(0.05, "cm")
    )
}


allhbcells_psthplot <-
  plot_grid(plotlist = rev(allsupercells_plotz),
            nrow = 6)


#### plots ####
panel_4B <- allzfcells_psthplot
panel_4A <- allhbcells_psthplot
