#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input PAF file [required]", metavar="FILE"),
  
  make_option(c("-o", "--output"), type="character", default="synteny_plot.png",
              help="Output plot file [default: synteny_plot.png]", metavar="FILE"),
  
  make_option(c("-q", "--min_quality"), type="integer", default=20,
              help="Minimum mapping quality [default: 20]", metavar="INT"),
  
  make_option(c("-l", "--min_length"), type="integer", default=1000,
              help="Minimum alignment length for plotting [default: 1000]", metavar="INT"),
  
  make_option(c("-f", "--min_filter_length"), type="integer", default=100,
              help="Minimum alignment length for initial filtering [default: 100]", metavar="INT"),
  
  make_option(c("-n", "--num_chromosomes"), type="integer", default=NULL,
              help="Number of top chromosomes to include (by size) [default: all]", metavar="INT"),
  
  make_option(c("-w", "--width"), type="double", default=12,
              help="Plot width in inches [default: 12]", metavar="FLOAT"),
  
  make_option(c("--height"), type="double", default=10,
              help="Plot height in inches [default: 10]", metavar="FLOAT"),
  
  make_option(c("--min_point_size"), type="double", default=0.1,
              help="Minimum point size [default: 0.1]", metavar="FLOAT"),
  
  make_option(c("--max_point_size"), type="double", default=3.0,
              help="Maximum point size [default: 3.0]", metavar="FLOAT"),
  
  make_option(c("--title"), type="character", default="",
              help="Plot title [default: auto-generated]", metavar="STRING"),
  
  make_option(c("--format"), type="character", default="png",
              help="Output format: pdf, png, svg [default: png]", metavar="STRING"),
  
  make_option(c("--dpi"), type="integer", default=300,
              help="DPI for raster formats [default: 300]", metavar="INT"),
  
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="Verbose output"),
  
  make_option(c("--x_label"), type="character", default="Target Genome",
              help="X-axis label [default: Target Genome]", metavar="STRING"),
  
  make_option(c("--y_label"), type="character", default="Query Genome",
              help="Y-axis label [default: Query Genome]", metavar="STRING"),
  
  make_option(c("--cdna_paf"), type="character", default=NULL,
              help="PAF file with cDNA mappings to target genome [optional]", metavar="FILE"),
  
  make_option(c("--cdna_paf_query"), type="character", default=NULL,
              help="PAF file with cDNA mappings to query genome [optional]", metavar="FILE"),
  
  make_option(c("--cdna_id"), type="character", default=NULL,
              help="File with cDNA ID mappings (snp, chr, pos, phenotype columns) [optional]", metavar="FILE"),
  
  make_option(c("--cdna_size"), type="double", default=2,
              help="Size of cDNA markers [default: 2]", metavar="FLOAT"),
  
  make_option(c("--cdna_alpha"), type="double", default=0.8,
              help="Transparency of cDNA markers [default: 0.8]", metavar="FLOAT")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Create D-GENIES style synteny plots from PAF alignment files with dots and optional cDNA mapping", 
                          add_help_option=FALSE)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input PAF file is required (--input)", call.=FALSE)
}

if (!file.exists(opt$input)) {
  stop(paste("Input file does not exist:", opt$input), call.=FALSE)
}

# Check cDNA files
if (!is.null(opt$cdna_paf) && !file.exists(opt$cdna_paf)) {
  stop(paste("cDNA PAF file (target) does not exist:", opt$cdna_paf), call.=FALSE)
}

if (!is.null(opt$cdna_paf_query) && !file.exists(opt$cdna_paf_query)) {
  stop(paste("cDNA PAF file (query) does not exist:", opt$cdna_paf_query), call.=FALSE)
}

if (!is.null(opt$cdna_id) && !file.exists(opt$cdna_id)) {
  stop(paste("cDNA ID file does not exist:", opt$cdna_id), call.=FALSE)
}

vlog <- function(msg) {
  if (opt$verbose) {
    cat(paste("[", Sys.time(), "] ", msg, "\n", sep=""))
  }
}

##########################################################################################
## Functions
##########################################################################################

read_paf_flexible <- function(file) {
  vlog(paste("Reading PAF file:", file))
  
  first_line <- readLines(file, n = 1)
  num_cols <- length(strsplit(first_line, "\t")[[1]])
  vlog(paste("Detected", num_cols, "columns in PAF file"))
  
  col_names <- c("query_name", "query_length", "query_start", "query_end", "strand",
                "target_name", "target_length", "target_start", "target_end",
                "matches", "alignment_length", "mapping_quality")
  
  if (num_cols > 12) {
    extra_cols <- paste0("extra_", 1:(num_cols - 12))
    col_names <- c(col_names, extra_cols)
  }
  
  paf_data <- read.table(file, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                        col.names=col_names, fill=TRUE, comment.char="")
  
  numeric_cols <- c("query_length", "query_start", "query_end", "target_length", 
                   "target_start", "target_end", "matches", "alignment_length", "mapping_quality")
  paf_data[numeric_cols] <- lapply(paf_data[numeric_cols], as.numeric)
  
  vlog(paste("Read", nrow(paf_data), "alignments"))
  return(paf_data)
}

read_cdna_id_file <- function(file) {
  vlog(paste("Reading cDNA ID file:", file))
  
  # Read the file and check if it has headers
  first_line <- readLines(file, n = 1)
  if (grepl("snp", first_line, ignore.case = TRUE)) {
    # File has headers
    cdna_data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    # File doesn't have headers, assume the format from your example
    cdna_data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                           col.names = c("snp", "chr", "pos", "phenotype"))
  }
  
  vlog(paste("Read", nrow(cdna_data), "cDNA entries"))
  return(cdna_data)
}

process_cdna_mappings <- function(cdna_paf_target, cdna_paf_query, cdna_id_file, target_chrom_info, query_chrom_info) {
  if (is.null(cdna_id_file)) {
    return(NULL)
  }
  
  vlog("Processing cDNA mappings...")
  
  # Read cDNA ID file
  cdna_ids <- read_cdna_id_file(cdna_id_file)
  
  cdna_target_data <- NULL
  cdna_query_data <- NULL
  
  # Process target genome cDNA mappings
  if (!is.null(cdna_paf_target)) {
    vlog("Processing target genome cDNA mappings...")
    cdna_paf_t <- read_paf_flexible(cdna_paf_target)
    
    cdna_target_merged <- cdna_paf_t %>%
      left_join(cdna_ids, by = c("query_name" = "snp"), relationship = "many-to-many") %>%
      filter(!is.na(phenotype))
    
    if (nrow(cdna_target_merged) > 0) {
      cdna_target_data <- cdna_target_merged %>%
        left_join(target_chrom_info, by = "target_name") %>%
        filter(!is.na(cumulative_start)) %>%
        mutate(
          adjusted_target_start = target_start + cumulative_start,
          adjusted_target_end = target_end + cumulative_start,
          adjusted_target_mid = (adjusted_target_start + adjusted_target_end) / 2,
          genome = "target"
        )
      vlog(paste("Processed", nrow(cdna_target_data), "target cDNA mappings"))
    }
  }
  
  # Process query genome cDNA mappings
  if (!is.null(cdna_paf_query)) {
    vlog("Processing query genome cDNA mappings...")
    cdna_paf_q <- read_paf_flexible(cdna_paf_query)
    
    cdna_query_merged <- cdna_paf_q %>%
      left_join(cdna_ids, by = c("query_name" = "snp"), relationship = "many-to-many") %>%
      filter(!is.na(phenotype))
    
    if (nrow(cdna_query_merged) > 0) {
      cdna_query_data <- cdna_query_merged %>%
        left_join(query_chrom_info, by = c("target_name" = "query_name")) %>%
        filter(!is.na(cumulative_start)) %>%
        mutate(
          adjusted_query_start = target_start + cumulative_start,
          adjusted_query_end = target_end + cumulative_start,
          adjusted_query_mid = (adjusted_query_start + adjusted_query_end) / 2,
          genome = "query"
        )
      vlog(paste("Processed", nrow(cdna_query_data), "query cDNA mappings"))
    }
  }
  
  return(list(target = cdna_target_data, query = cdna_query_data))
}

plot_synteny_dgenies_style <- function(paf_data, min_length = 1000, min_point_size = 0.1, max_point_size = 3.0,
                                      x_label = "Target Genome", y_label = "Query Genome",
                                      cdna_data = NULL, cdna_size = 2, cdna_alpha = 0.8) {
    
  paf_clean <- paf_data %>%
    filter(alignment_length >= min_length) %>%
    mutate(
      identity = matches / alignment_length,
      strand_type = ifelse(strand == "+", "Forward", "Reverse")
    )
  
  vlog(paste("After length filtering:", nrow(paf_clean), "alignments"))
  
  # Calculate scaling factors for size and alpha
  mapped_bases_range <- range(paf_clean$matches, na.rm = TRUE)
  quality_range <- range(paf_clean$mapping_quality, na.rm = TRUE)
  
  vlog(paste("Mapped bases range:", mapped_bases_range[1], "-", mapped_bases_range[2]))
  vlog(paste("Mapping quality range:", quality_range[1], "-", quality_range[2]))
  
  orientation_stats <- paf_clean %>%
    group_by(target_name, query_name, strand) %>%
    summarise(
      total_length = sum(alignment_length),
      count = n(),
      .groups = 'drop'
    ) %>%
    group_by(target_name, query_name) %>%
    slice_max(total_length, n = 1) %>%
    select(target_name, query_name, predominant_strand = strand)

  # order y numerically 
  target_chrom_info <- paf_clean %>%
    group_by(target_name) %>%
    summarise(
      max_target_length = max(target_length),
      .groups = 'drop'
    ) %>%
    mutate(
      scaffold_num = as.numeric(gsub("\\D", "", target_name))
    ) %>%
    arrange(scaffold_num) %>%  # Simple numerical order: 1, 2, 3, 4, 5...
    mutate(
      cumulative_start = lag(cumsum(max_target_length), default = 0),
      cumulative_end = cumsum(max_target_length),
      midpoint = (cumulative_start + cumulative_end) / 2
    )
  
  vlog(paste("Target chromosomes (y-axis, numerical order):", paste(head(target_chrom_info$target_name, 20), collapse=", "), "..."))
  
  # order x by synteny
  query_chrom_info <- paf_clean %>%
    left_join(target_chrom_info, by = "target_name") %>%
    mutate(
      linearized_target_pos = target_start + cumulative_start
    ) %>%
    group_by(query_name) %>%
    summarise(
      max_query_length = max(query_length),
      median_target_position = median(linearized_target_pos, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(median_target_position) %>% 
    mutate(
      cumulative_start = lag(cumsum(max_query_length), default = 0),
      cumulative_end = cumsum(max_query_length),
      midpoint = (cumulative_start + cumulative_end) / 2
    )
  
  vlog(paste("Query chromosomes (x-axis, synteny order):", paste(head(query_chrom_info$query_name, 20), collapse=", "), "..."))
  
  # Process cDNA data if provided
  cdna_processed <- process_cdna_mappings(opt$cdna_paf, opt$cdna_paf_query, opt$cdna_id, target_chrom_info, query_chrom_info)
  
  paf_oriented <- paf_clean %>%
    left_join(orientation_stats, by = c("target_name", "query_name")) %>%
    left_join(target_chrom_info, by = "target_name", suffix = c("", "_target")) %>%
    left_join(query_chrom_info, by = "query_name", suffix = c("", "_query")) %>%
    mutate(
      corrected_query_start = ifelse(predominant_strand == "-", 
                                     query_length - query_end, 
                                     query_start),
      corrected_query_end = ifelse(predominant_strand == "-", 
                                   query_length - query_start, 
                                   query_end),
      adjusted_target_start = target_start + cumulative_start,
      adjusted_target_end = target_end + cumulative_start,
      adjusted_query_start = corrected_query_start + cumulative_start_query,
      adjusted_query_end = corrected_query_end + cumulative_start_query,
      adjusted_target_mid = (adjusted_target_start + adjusted_target_end) / 2,
      adjusted_query_mid = (adjusted_query_start + adjusted_query_end) / 2,
      target_name = factor(target_name, levels = target_chrom_info$target_name),
      query_name = factor(query_name, levels = query_chrom_info$query_name)
    )
      
  p <- ggplot(paf_oriented, aes(x = adjusted_target_mid, y = adjusted_query_mid, 
                               size = matches, alpha = mapping_quality)) +
    geom_point(color = "black") +
    
    geom_vline(xintercept = target_chrom_info$cumulative_end, color = "lightgrey", alpha = 0.6, linewidth = 0.5, linetype = "dashed") +
    geom_hline(yintercept = query_chrom_info$cumulative_end, color = "lightgrey", alpha = 0.6, linewidth = 0.5, linetype = "dashed") +

    scale_size_continuous(
      name = "Mapped Bases",
      range = c(min_point_size, max_point_size),
      labels = function(x) {
        ifelse(x >= 1e6, paste0(round(x/1e6, 1), "M"),
               ifelse(x >= 1e3, paste0(round(x/1e3, 1), "K"), 
                      as.character(x)))
      },
      guide = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(alpha = 0.7),
        order = 1
      )
    ) +
    
    scale_alpha_continuous(
      name = "Mapping Quality",
      range = c(0.1, 0.9),
      guide = guide_legend(
        title.position = "top", 
        title.hjust = 0.5,
        override.aes = list(size = 2),
        order = 2
      )
    ) +
    
    scale_x_continuous(
      labels = function(x) paste0(round(x/1e6, 1), "M"),
      n.breaks = 8,
      guide = guide_axis(check.overlap = TRUE),
      sec.axis = sec_axis(~ ., 
                         breaks = target_chrom_info$midpoint,
                         labels = target_chrom_info$target_name,
                         name = NULL,
                         guide = guide_axis(check.overlap = TRUE))
    ) +
    
    scale_y_continuous(
      labels = function(y) paste0(round(y/1e6, 1), "M"),
      n.breaks = 8,
      guide = guide_axis(check.overlap = TRUE),
      sec.axis = sec_axis(~ ., 
                         breaks = query_chrom_info$midpoint,
                         labels = query_chrom_info$query_name,
                         name = NULL,
                         guide = guide_axis(check.overlap = TRUE))
    ) +
    
    labs(
      x = x_label,
      y = y_label
    ) +
    
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      #panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      #panel.grid.minor = element_line(color = "gray80", linewidth = 0.2),
      panel.grid.major = element_blank(),
      panel.grid.minor =  element_blank(),
      axis.text = element_text(size = 8, color = "black"),
      axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7),
      axis.text.y.right = element_text(angle = 0, hjust = 0, size = 7),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.title = element_text(size = 10, color = "black"),
      plot.margin = margin(5, 5, 5, 5),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.margin = margin(0, 0, 0, 10),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.8, "lines"),
      legend.spacing.y = unit(0.2, "lines")
    )
  
  if (!is.null(cdna_processed)) {
    cdna_target_data <- cdna_processed$target
    cdna_query_data <- cdna_processed$query
    
    all_cdna_data <- bind_rows(cdna_target_data, cdna_query_data)
    
    if (!is.null(all_cdna_data) && nrow(all_cdna_data) > 0) {
      vlog("Adding cDNA markers to plot...")
      
      unique_phenotypes <- unique(all_cdna_data$phenotype)
      n_phenotypes <- length(unique_phenotypes)
      
      if (n_phenotypes <= 8) {
        colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
      } else {
        colors <- rainbow(n_phenotypes)
      }
      
      phenotype_colors <- setNames(colors[1:n_phenotypes], unique_phenotypes)
      
      x_range <- range(paf_oriented$adjusted_target_mid, na.rm = TRUE)
      y_range <- range(paf_oriented$adjusted_query_mid, na.rm = TRUE)
      
      x_bottom <- x_range[1] - (x_range[2] - x_range[1]) * 0.05  # 5% to the left
      y_bottom <- y_range[1] - (y_range[2] - y_range[1]) * 0.05  # 5% below the bottom
      
      if (!is.null(cdna_target_data) && nrow(cdna_target_data) > 0) {
        p <- p + 
          geom_point(data = cdna_target_data,
                    aes(x = adjusted_target_mid, y = y_bottom, color = phenotype),
                    size = cdna_size, alpha = cdna_alpha, shape = 17, inherit.aes = FALSE)
      }
      
      if (!is.null(cdna_query_data) && nrow(cdna_query_data) > 0) {
        p <- p + 
          geom_point(data = cdna_query_data,
                    aes(x = x_bottom, y = adjusted_query_mid, color = phenotype),
                    size = cdna_size, alpha = cdna_alpha, shape = 15, inherit.aes = FALSE)  # Square shape
      }
      
      if (!is.null(cdna_target_data) && !is.null(cdna_query_data) && 
          nrow(cdna_target_data) > 0 && nrow(cdna_query_data) > 0) {
        
        shared_cdnas <- inner_join(
          cdna_target_data %>% select(query_name, phenotype, adjusted_target_mid) %>% distinct(),
          cdna_query_data %>% select(query_name, phenotype, adjusted_query_mid) %>% distinct(),
          by = c("query_name", "phenotype"),
          relationship = "many-to-many"
        )
        
        if (nrow(shared_cdnas) > 0) {
          vlog(paste("Adding", nrow(shared_cdnas), "connecting lines for shared cDNAs"))
          
          p <- p + 
            geom_segment(data = shared_cdnas,
                        aes(x = adjusted_target_mid, y = y_bottom,
                            xend = adjusted_target_mid, yend = adjusted_query_mid,
                            color = phenotype),
                        alpha = 0.6, linewidth = 0.5, inherit.aes = FALSE)
          
          p <- p + 
            geom_segment(data = shared_cdnas,
                        aes(x = adjusted_target_mid, y = adjusted_query_mid,
                            xend = x_bottom, yend = adjusted_query_mid,
                            color = phenotype),
                        alpha = 0.6, linewidth = 0.5, inherit.aes = FALSE)
        }
      }
      
      p <- p + 
        scale_color_manual(values = phenotype_colors, name = "cDNA Type") +
        guides(
          size = guide_legend(order = 1, override.aes = list(alpha = 0.7)),
          alpha = guide_legend(order = 2, override.aes = list(size = 2)),
          color = guide_legend(order = 3, override.aes = list(size = 3, alpha = 1))
        )
      
      vlog(paste("Added cDNA markers:",
                ifelse(!is.null(cdna_target_data), nrow(cdna_target_data), 0), "target,",
                ifelse(!is.null(cdna_query_data), nrow(cdna_query_data), 0), "query"))
    }
  }
  
  return(p)
}

##########################################################################################
## Main execution
##########################################################################################

paf_data <- read_paf_flexible(opt$input)

vlog(paste("Applying filters: mapping_quality >", opt$min_quality, ", alignment_length >", opt$min_filter_length))
paf_filtered <- paf_data[paf_data$mapping_quality > opt$min_quality & 
                        paf_data$alignment_length > opt$min_filter_length, ]

vlog(paste("After quality filtering:", nrow(paf_filtered), "alignments"))

if (!is.null(opt$num_chromosomes)) {
  vlog(paste("Selecting top", opt$num_chromosomes, "chromosomes by size"))
  paf_order <- unique(paf_filtered[order(paf_filtered$target_length, decreasing = T), "target_name"])
  paf_filtered <- paf_filtered[paf_filtered$target_name %in% paf_order[1:opt$num_chromosomes], ]
} else {
  paf_order <- unique(paf_filtered[order(paf_filtered$target_length, decreasing = T), "target_name"])
  paf_filtered <- paf_filtered[paf_filtered$target_name %in% paf_order, ]
}

vlog(paste("Final dataset:", nrow(paf_filtered), "alignments"))

synteny_plot <- plot_synteny_dgenies_style(
  paf_filtered, 
  min_length = opt$min_length,
  min_point_size = opt$min_point_size,
  max_point_size = opt$max_point_size,
  x_label = opt$x_label,
  y_label = opt$y_label,
  cdna_size = opt$cdna_size,
  cdna_alpha = opt$cdna_alpha
)

if (opt$title != "") {
  synteny_plot <- synteny_plot + ggtitle(opt$title) +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
}

vlog(paste("Saving plot to:", opt$output))

file_ext <- tolower(tools::file_ext(opt$output))
if (file_ext == "" || file_ext != opt$format) {
  opt$output <- paste0(tools::file_path_sans_ext(opt$output), ".", opt$format)
}

if (opt$format == "pdf") {
  ggsave(opt$output, synteny_plot, width = opt$width, height = opt$height, device = "pdf")
} else if (opt$format == "png") {
  ggsave(opt$output, synteny_plot, width = opt$width, height = opt$height, dpi = opt$dpi, device = "png", bg = "white", scale = 0.8)
} else if (opt$format == "svg") {
  ggsave(opt$output, synteny_plot, width = opt$width, height = opt$height, device = "svg")
} else {
  stop(paste("Unsupported format:", opt$format), call.=FALSE)
}

vlog("Done.")