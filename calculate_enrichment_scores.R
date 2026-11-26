library(massiveGST)
library(progress)


omics_names <- c("phospho", "pro", 'rna')
n_markers_values <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000)

data_dir <- "/lustre/home/acct-medsjia/medsjia-wuzr/yy/prcc_multiomics/results/raw/"
signatures_dir <- "/lustre/home/acct-medsjia/medsjia-wuzr/yy/prcc_multiomics/results/signatures/"
save_dir <- "/lustre/home/acct-medsjia/medsjia-wuzr/yy/prcc_multiomics/results/enrichment_scores/"

for (omics_name in omics_names) {
    # load data
    data_file <- paste0(data_dir, omics_name, ".csv")
    signatures_files_basename <- paste0(signatures_dir, omics_name, "_nmf")
    data <- read.csv(data_file, sep = ",", header = TRUE)
    rownames(data) <- data[[1]]
    data <- data[-1]
    markers <- list()
    for (i in 1:3) {
        markers[[i]] <- read.csv(paste0(signatures_files_basename, i, ".csv"), sep = ",", header = TRUE)
        rownames(markers[[i]]) <- markers[[i]][[1]]
        markers[[i]] <- markers[[i]][-1]
        markers[[i]] = markers[[i]][order(markers[[i]][["FDR"]]),]
    }


    pbar <- progress_bar$new(total = length(n_markers_values), format = paste0("Processing ", toupper(omics_name), " [:bar] :percent eta: :eta"), clear = FALSE)
    for (n_markers in n_markers_values) {
        pbar$tick()
        # get marker genes
        marker_genes <- list()
        for (i in 1:3) {
            marker_genes[[i]] <- rownames(markers[[i]])[1:n_markers]
        }

        # calculate enrichment scores
        enrichment_scores <- list(data.frame(), data.frame(), data.frame())
        for (sample in colnames(data)) {
            gene_exp <- data[, sample]
            names(gene_exp) <- rownames(data)
            res <- massiveGST(gene_exp, marker_genes)
            enrichment_scores[[1]][sample, colnames(res)] <- res["1",]
            enrichment_scores[[2]][sample, colnames(res)] <- res["2",]
            enrichment_scores[[3]][sample, colnames(res)] <- res["3",]
        }

        output_dir <- paste0(save_dir, "/n_markers_", n_markers, "/", omics_name, "/")
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        for (i in 1:3) {
            write.csv(enrichment_scores[[i]], paste0(output_dir, "/nmf", i, ".csv"), row.names = TRUE)
        }
    }
}

print("Done!")
