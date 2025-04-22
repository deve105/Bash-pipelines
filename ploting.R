pacman::p_load(tidyverse, tidyplots, zoo, signal)

global_cov = read_tsv("data/2504_mosdepth_global_perbase.tsv", col_names=FALSE) |>
    dplyr::rename(position=2, coverage=4, sample=5) |>
    dplyr::group_by(sample) |>
    #dplyr::filter(position>776, position<8301) |>
    dplyr::mutate(smooth=zoo::rollmedian(coverage, k = 101, fill = NA)) |>
    dplyr::ungroup(sample) |>
    dplyr::mutate(cohort=ifelse(str_detect(sample, "IRID|SS"), "Peru", "Miyazaki"))

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
