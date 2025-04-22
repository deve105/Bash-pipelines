pacman::p_load(tidyverse, tidyplots, zoo, signal)

read_tsv("data/2504_combined_output.tsv", col_names=FALSE) |>
    dplyr::rename(position=2, coverage=4, sample=5) |>
    
    tidyplot(x=position, y=coverage) |>
    add_median_line() #|>
    tidyplots::save_plot(
        filename = "2504_combined_output.png",
        background = "transparent",
        dpi = 300
    )


datatset = read_tsv("2504_combined_output.tsv", col_names=FALSE) |>
    dplyr::rename(position=2, coverage=4, sample=5) 


tidyplot(dataset, x=position, y=coverage) |>
    add_median_line() |>
    tidyplots::save_plot(
        filename = "2504_combined_output.png",
        background = "transparent",
        dpi = 300
    )


samples <- unique(datatset$sample)
samples
batch_size <- 10
batches <- split(samples, ceiling(seq_along(samples) / batch_size))
batches
for (i in seq_along(batches)) {
    batch <- batches[[i]]
    batch_data <- datatset %>% filter(sample %in% batch)
    
    tidyplot(batch_data, x=position, y=coverage) |>
        add_median_area() |>
        add_title(paste0("Batch", i)) |>
        split_plot(by=sample,
            ncol = 2,
            nrow = 5
            ) |> 
        save_plot(
        filename = paste0("2504_combined_output_batch_", i, ".png"),
        background = "transparent",
        dpi = 300
    )
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
