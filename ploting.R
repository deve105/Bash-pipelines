pacman::p_load(tidyverse, tidyplots, zoo, signal)


## To add HTLV-1 coverage to our dataset in google drive
coverage_htlv = read_tsv("data/2504_mosdepth_global.tsv")
htlvcov = coverage_htlv |>
    dplyr::filter(chrom == "total") |>
    dplyr::rename(sample = `22_SS01`) |>
    dplyr::mutate(sample = stringr::str_remove(sample, "^[0-9]+_")) |>
    dplyr::filter(!stringr::str_detect(sample, "ATL|1339")) |>
    dplyr::select(mean, sample)

## Reading the database to save it later
htlv_db = googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/14eZCI8lRsfIebXu-tpAxbuPEOhE6moAnTDKKKUQyNOI/edit?gid=1412798579#gid=1412798579", sheet = "Sheet3")

## Left joining
newdb = htlv_db |>
    left_join(htlvcov, by=c("ID2"="sample"))
## Writing Sheet in Google Drive
googlesheets4::write_sheet(newdb, ss="https://docs.google.com/spreadsheets/d/14eZCI8lRsfIebXu-tpAxbuPEOhE6moAnTDKKKUQyNOI/edit#gid=1412798579", sheet = "Sheet3")


## Re-reading the dataset
htlv_db = googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/14eZCI8lRsfIebXu-tpAxbuPEOhE6moAnTDKKKUQyNOI/edit?gid=1412798579#gid=1412798579", sheet = "Sheet3")

htlv_db
## Calculating the correlation between PVl and meanCOV
htlv_db = htlv_db |>
    mutate(meanCOV=as.numeric(meanCOV),
           PVL=as.numeric(PVL2)) 

htlv_db |>
tidyplots:: tidyplot(x=meanCOV, y=PVL2) |>
    tidyplots::add_data_points_jitter() |>
    tidyplots::adjust_x_axis_title("Mean HTLV-1 Coverage", face="bold") |>
    tidyplots::add_test_pvalue(test="spearman") |>
    tidyplots::adjust_y_axis_title("Proviral Load (% PBMC)", face="bold") # |>
    tidyplots::adjust_size(width = 8, height = 7, unit = "cm") #|>
    tidyplots::save_plot(
        filename = "2504_PVL_vs_MeanCOV.tiff",
        background = "transparent",
        dpi = 1200
    )


labelx = paste("Spearman rho =", 
                          round(cor(htlv_db$meanCOV, htlv_db$PVL, 
                                    method="spearman", use="complete.obs"), 2),
                          "\np = <2.2e-16") #, 
                          #round(cor.test(htlv_db$meanCOV, htlv_db$PVL, 
                                        #method="spearman")$p.value, 7))
labelx
cor.test(htlv_db$meanCOV, htlv_db$PVL, 
                                        method="spearman")
htlv_db |>
    tidyplots::tidyplot(x=meanCOV, y=PVL, color=disease_2) |>
    tidyplots::add_data_points_beeswarm(size=1.5, preserve="total", alpha=.8) |>
    tidyplots::add(geom_smooth(
      method = "lm", 
      se = FALSE, 
      color = "black", 
      linewidth=.5,
      inherit.aes=FALSE, aes(x=meanCOV, y=PVL),
      alpha = 0.4
      )) |>
    tidyplots::add_annotation_text(
      labelx,
      x = min(htlv_db$meanCOV, na.rm=TRUE) + 180,
      y = max(htlv_db$PVL, na.rm=TRUE) - 1,
    ) |>
     tidyplots::adjust_x_axis_title("Mean HTLV-1 Coverage", face="bold") |>
    tidyplots::adjust_y_axis_title("Proviral Load (% PBMC)", face="bold") |>
    tidyplots::adjust_legend_title("Last Disease \nStatus") |>
    tidyplots::adjust_colors(new_colors=c("#377EB8", "#E41A1C", "#a6611a")) |>
    tidyplots::adjust_font(face="bold") |>
    tidyplots::save_plot("output/2508_PVL_vs_Coverage.tiff", bg="transparent", dpi = 600)
### PVL vs Gini 
labely = paste("Spearman rho =", 
                          round(cor(htlv_db$gini, htlv_db$PVL, 
                                    method="spearman", use="complete.obs"), 2),
                          "\np =",
                          round(cor.test(htlv_db$gini, htlv_db$PVL, 
                                        method="spearman")$p.value, 3))
labely
cor.test(htlv_db$meanCOV, htlv_db$PVL, 
                                        method="spearman")

htlv_db |>
    tidyplots::tidyplot(x=gini, y=PVL, color=disease_2) |>
    tidyplots::add_data_points_beeswarm(size=1.5, preserve="total", alpha=.8) |>
    tidyplots::add(geom_smooth(
      method = "lm", 
      se = FALSE, 
      color = "black", 
      linewidth=.5,
      inherit.aes=FALSE, aes(x=gini, y=PVL),
      alpha = 0.4
      )) |>
    tidyplots::add_annotation_text(
      labely,
      x = min(htlv_db$gini, na.rm=TRUE) +0.3,
      y = max(htlv_db$PVL, na.rm=TRUE) - 1,
    ) |>
     tidyplots::adjust_x_axis_title("Oligoclonality Index (OCI)", face="bold") |>
    tidyplots::adjust_y_axis_title("Proviral Load (% PBMC)", face="bold") |>
    tidyplots::adjust_legend_title("Last Disease \nStatus") |>
    tidyplots::adjust_colors(new_colors=c("#377EB8", "#E41A1C", "#a6611a")) |>
    tidyplots::adjust_font(face="bold") |>
    tidyplots::save_plot("output/2508_PVL_vs_gini.tiff", bg="transparent", dpi = 600)




global_cov = read_tsv("data/2504_mosdepth_global_perbase.tsv", col_names=FALSE) |>
    dplyr::rename(position=2, coverage=4, sample=5) |>
    dplyr::group_by(sample) |>
    #dplyr::filter(position>776, position<8301) |>
    dplyr::mutate(smooth=zoo::rollmedian(coverage, k = 101, fill = NA)) |>
    dplyr::ungroup(sample) |>
    dplyr::mutate(cohort=ifelse(str_detect(sample, "IRID|SS"), "Peru", "Miyazaki"))

global_cov
new_colors=c("#E69F00","#C02D45")
global_cov |>
    tidyplot(x=position, y=smooth, color=cohort) |>
    add_median_area(alpha=.5) |>
    adjust_y_axis(transform = "log10", limits=c(10,1000), labels = scales::trans_format("log10", scales::math_format(10^.x))) |>
    adjust_x_axis(breaks=c(0, 820, 2780, 5200, 6670, 8380, 9050), limits=c(0,9050)) |>
    tidyplots::adjust_y_axis_title("Median HTLV-1 Coverage", face="bold") |>
    tidyplots::adjust_x_axis_title("HTLV-1 Genome Position", face="bold") |>
    adjust_colors(new_colors=c("#E69F00","#C02D45"), labels=c("Peru", "Miyazaki-Japan")) |> 
    remove_legend() |>
    tidyplots::adjust_size(width = 12, height = 7, unit = "cm") |>
    tidyplots::save_plot(
        filename = "2504_supplementaryC_wholecov_peru_japan.png",
        background = "transparent",
        dpi = 300
    )
defective = c("IRID006", "IRID086", "IRID067")

global_cov |>
    dplyr::filter(sample%in%defective) |>
    tidyplot(x=position, y=smooth) |>
    add_median_area(alpha=.5) |>
    adjust_y_axis(transform = "log10", limits=c(10,1500), labels = scales::trans_format("log10", scales::math_format(10^.x))) |>
    adjust_x_axis(breaks=c(0, 2780, 5200, 9050), limits=c(0,9100)) |>
    tidyplots::adjust_y_axis_title("HTLV-1 Coverage", face="bold") |>
    tidyplots::remove_x_axis_title() |>
    tidyplots::adjust_title(face="bold") |>
    tidyplots::adjust_colors("#C02D45") |>
    tidyplots::adjust_size(width = 6, height = 15, unit = "cm") |>
    tidyplots::split_plot(by="sample", guides="auto", ncol=1, nrow=3) |>
    tidyplots::save_plot(
        filename = "2504_deffectiveviruses.png",
        background = "transparent",
        dpi = 300
    )


samples <- unique(global_cov$sample)
samples
batch_size <- 28
batches <- split(samples, ceiling(seq_along(samples) / batch_size))

batches
for (i in seq_along(batches)) {
    batch <- batches[[i]]
    batch_data <- global_cov |>
        dplyr::filter(sample %in% batch) |>
        mutate(smooth=smooth+1)
    tidyplot(batch_data, x=position, y=smooth) |>
    add_median_area(alpha=.5) |>
    adjust_y_axis(transform = "log10", limits=c(1,1500), labels = scales::trans_format("log10", scales::math_format(10^.x))) |>
    adjust_x_axis(breaks=c(0, 2780, 5200, 9050), limits=c(0,9100)) |>
    tidyplots::adjust_y_axis_title("HTLV-1 Coverage", face="bold") |>
    tidyplots::remove_x_axis_title() |>
    tidyplots::adjust_title(face="bold") |>
    tidyplots::adjust_colors("#C02D45") |>
    tidyplots::split_plot(by=sample, guides="auto", ncol=4, nrow=7) |>
    save_plot(
        filename = paste0("2504_combined_", i, ".png"),
        background = "transparent",
        dpi = 300)
}


## Plotting by regions

regions_cov <- read_tsv("2504new_combined_output.tsv", col_names = FALSE) |>
    dplyr::rename(position = 2, region = 4, coverage = 5, sample = 6) |>
    mutate(country=case_when(
        grepl("IRID", sample) ~ "Peru",
        grepl("ATLSkin", sample) ~ "Japan",
        .default = NA
    )) 
unique(regions_cov$country)

regions_cov |>
    filter(!is.na(country)) |>
    tidyplot(x=region, y=coverage, color=country) |>
    add_median_dot() |>
    add_ci95_errorbar() |>
    adjust_x_axis(rotate_labels  = 90) |>
    save_plot(
        filename = "2504_combined_Japanese_output_regions.png",
        background = "transparent",
        dpi = 300
    )
