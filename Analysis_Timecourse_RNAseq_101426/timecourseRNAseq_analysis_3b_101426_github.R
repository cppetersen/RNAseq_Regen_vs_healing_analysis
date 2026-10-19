library(limma)
library(splines)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
##*##
setwd("~/Documents/StudentData/Elias/RNAseq_analysis/analysis101426_timecourseRNAseq_github_test/")
#


head<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis090126/ddv6_HeadR_031225.tsv",sep="\t", header=FALSE)
tail<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis090126/ddv6_TailR_031225.tsv",sep="\t", header=FALSE)
heal<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis090126/ddv6_Heal_031225.tsv",sep="\t", header=FALSE)

colnames(head)<-c("Contig","experiment","Experiment.2","time","log2FC")
colnames(tail)<-c("Contig","experiment","Experiment.2","time","log2FC")
colnames(heal)<-c("Contig","experiment","Experiment.2","time","log2FC")

head.rename<-separate_wider_delim(head,time,delim=" ",names=c("hours","tag"))
tail.rename<-separate_wider_delim(tail,time,delim=" ",names=c("hours","tag"))
heal.rename<-separate_wider_delim(heal,time,delim=" ",names=c("hours","tag"))


#head.rename$hours<-as.numeric(head.rename$hours)

head.wide<-pivot_wider(head.rename,names_from=hours,values_from=c(log2FC))
head.wide<-pivot_wider(head,names_from=time,values_from=c(log2FC))
head.wide<-select(head.wide,1,2,4,5,8,11,13,14,6,7,9,10,12)
colnames(head.wide)<-c("gene","condition","t_0","t_1","t_2","t_4","t_6","t_8","t_12","t_16","t_24","t_36","t_48")
head.long<-head.wide %>% pivot_longer(3:13,names_to="timepoint")
head.long<-head.long %>% separate_wider_delim(timepoint, delim="_", names = c("dum","time"))


tail.wide<-pivot_wider(tail.rename,names_from=hours,values_from=c(log2FC))
tail.wide<-pivot_wider(tail,names_from=time,values_from=c(log2FC))
tail.wide<-select(tail.wide,1,2,4,5,8,11,13,14,6,7,9,10,12)
colnames(tail.wide)<-c("gene","condition","t_0","t_1","t_2","t_4","t_6","t_8","t_12","t_16","t_24","t_36","t_48")
tail.long<-tail.wide %>% pivot_longer(3:13,names_to="timepoint")
tail.long<-tail.long %>% separate_wider_delim(timepoint, delim="_", names = c("dum","time"))


heal.wide<-pivot_wider(heal.rename,names_from=hours,values_from=c(log2FC))
heal.wide<-pivot_wider(heal,names_from=time,values_from=c(log2FC))
heal.wide<-select(heal.wide,1,2,4,5,8,11,13,14,6,7,9,10,12)
colnames(heal.wide)<-c("gene","condition","t_0","t_1","t_2","t_4","t_6","t_8","t_12","t_16","t_24","t_36","t_48")
heal.long<-heal.wide %>% pivot_longer(3:13,names_to="timepoint")
heal.long<-heal.long %>% separate_wider_delim(timepoint, delim="_", names = c("dum","time"))

combined_df_long <- bind_rows(head.long, tail.long,heal.long) %>% select(gene, condition, time, value)
colnames(combined_df_long)<-c("gene","condition","time","log2FC")

Wurtzel.anno<-read.table("~/Documents/StudentData/Elias/RNAseq_analysis/analysis090126/Wurtzel_short_anno.txt",sep="\t",header=TRUE)

injury_genes<-Wurtzel.anno %>% pull(Contig)

replacements <- c("head" = "HeadRegen", "heal" = "Incision", "tail"="TailRegen")
combined_df_long <- combined_df_long %>% 
  mutate(condition = str_replace_all(condition, replacements))

combined_df_long$time<-as.numeric(combined_df_long$time)

df_long<-combined_df_long

df_long$condition <- factor(df_long$condition, levels = c("Incision","HeadRegen","TailRegen"))

## ------------------------------------------------------------
## 0. Expected data format
## ------------------------------------------------------------
## df_long: one row per gene x condition x timepoint, columns:
##   gene      : character/factor, gene identifier
##   condition : factor, "HeadRegen" / "TailRegen" / "Incision"
##   time      : numeric, hours (0, 1, 2, 4, 6, 8, 16, 24, 36, 48, ...)
##   log2FC    : numeric, pre-normalized to that condition's 0h
##
## df_long <- read.csv("your_data.csv")
## df_long$condition <- factor(df_long$condition,
##                              levels = c("Incision","HeadRegen","TailRegen"))

## injury_genes: character vector of your known injury-induced gene
## panel, matching values in df_long$gene
## injury_genes <- c("gene1","gene2", ...)

## ------------------------------------------------------------
## 1. Reshape to a genes x samples matrix (limma expects this
##    "expression matrix" format, same as it would for microarrays)
##    Each "sample" here is one condition-timepoint combination.
## ------------------------------------------------------------
df_long <- df_long %>%
  mutate(sample_id = paste(condition, time, sep = "_t"))

expr_mat <- df_long %>%
  select(gene, sample_id, log2FC) %>%
  pivot_wider(names_from = sample_id, values_from = log2FC) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

sample_info <- df_long %>%
  distinct(sample_id, condition, time) %>%
  arrange(condition, time)

expr_mat <- expr_mat[, sample_info$sample_id]  # ensure column order matches sample_info

# Drop genes with too much missing data to fit a spline reliably
# (require at least 80% of timepoints present)
keep_genes <- rowMeans(!is.na(expr_mat)) >= 0.8
expr_mat <- expr_mat[keep_genes, ]
cat("Genes retained after missingness filter:", nrow(expr_mat), "\n")

## ------------------------------------------------------------
## 2. GOAL 1: which genes are induced, WITHIN each condition
## ------------------------------------------------------------
## Fit a separate natural-spline model per condition with df=4. df=4 is a
## reasonable starting point for ~10 timepoints spanning 0-48h --
## enough flexibility to capture a rise-and-fall injury response,
## but still leaves residual degrees of freedom for the F-test.

spline_df <- 4
#log2FC_cutoff<-0.8

fit_condition_spline <- function(expr_mat, sample_info, cond, spline_df = 4) {
  cols <- sample_info$sample_id[sample_info$condition == cond]
  sub_mat <- expr_mat[, cols]
  sub_time <- sample_info$time[sample_info$condition == cond]
  
  design <- model.matrix(~ ns(sub_time, df = spline_df))
  
  fit <- lmFit(sub_mat, design)
  fit <- eBayes(fit)
  
  # F-test across all spline coefficients (excludes intercept) --
  # tests "is this gene's trajectory non-flat over time"
  tt <- topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = "none")
  
  # Direction/magnitude: fitted curve's peak value (induced = rises
  # above baseline at some point) -- use fitted values, not raw
  # data, since fitted values are the denoised smooth curve
  fitted_vals <- fitted(fit)
  peak_logfc  <- apply(fitted_vals, 1, max)
  trough_logfc <- apply(fitted_vals, 1, min)
  
  tibble::tibble(
    gene = rownames(sub_mat),
    condition = cond,
    F_pvalue = tt$P.Value,
    F_adj_pvalue = tt$adj.P.Val,
    peak_log2FC = peak_logfc,
    trough_log2FC = trough_logfc
  )
}

conditions <- levels(sample_info$condition)

induction_results <- map_dfr(conditions, function(cnd) {
  fit_condition_spline(expr_mat, sample_info, cnd, spline_df = spline_df)
})

## Classify genes: "induced" = significant non-flat trajectory AND
## peak log2FC exceeds a chosen threshold (e.g. 1 = 2-fold change)
induction_results <- induction_results %>%
  mutate(
    induced = F_adj_pvalue < 0.05 & peak_log2FC >= 0.8,
    repressed = F_adj_pvalue < 0.05 & trough_log2FC <= -0.8
  )

## Summary: how many induced/repressed genes per condition
induced_genes_by_condition_table<-induction_results %>%
  group_by(condition) %>%
  summarise(n_induced = sum(induced), n_repressed = sum(repressed), .groups = "drop")
write.table(induced_genes_by_condition_table,"induction_results2.csv",sep=",")

induction_results_wide <- induction_results %>%
  pivot_wider(
    id_cols = gene,
    names_from = condition,
    values_from = 3:8,
    names_glue = "{.value}_{condition}"
  )
#induction_results_wide <-mutate(HeadRegen_vs_Incision_log2FC = )

write.table(induction_results_wide,"induction_results_wide2.csv",sep=",")


## Full list of induced genes per condition, sorted by significance
induced_genes_by_condition <- induction_results %>%
  filter(induced) %>%
  arrange(condition, F_adj_pvalue)

print(induced_genes_by_condition)



#dd_Smed_v6_24180_0_1 is notum
#dd_Smed_v6_28398_0_1 is wnt1
#induction_results_wide %>% filter(gene %in% c("dd_Smed_v6_24180_0_1","dd_Smed_v6_28398_0_1"))

## ------------------------------------------------------------
## 2b. CROSS-CHECK (simple heuristic, no modeling assumptions):
##     require |log2FC| above threshold at >=2 CONSECUTIVE
##     timepoints. A true induction should persist across
##     adjacent timepoints; a single-point spike is more likely
##     measurement noise. This doesn't replace the spline test,
##     but genes that pass BOTH are the most trustworthy calls.
## ------------------------------------------------------------
consecutive_check <- df_long %>%
  filter(gene %in% induced_genes_by_condition$gene) %>%
  arrange(gene, condition, time) %>%
  group_by(gene, condition) %>%
  mutate(above_thresh = log2FC >= 0.7) %>%
  mutate(run_len = with(rle(above_thresh), rep(lengths, lengths))) %>%
  summarise(max_consecutive_run = max(run_len[above_thresh], 0), .groups = "drop")

induced_genes_by_condition <- induced_genes_by_condition %>%
  left_join(consecutive_check, by = c("gene", "condition")) %>%
  mutate(consistent_induction = max_consecutive_run >= 2)
write.table(induced_genes_by_condition,"induced_genes_by_condition2.csv",sep=",")


## Genes passing both the formal spline test AND the consecutive-
## timepoint consistency check -- your highest-confidence induced set
high_confidence_induced <- induced_genes_by_condition %>%
  filter(consistent_induction)

cat("High-confidence induced genes (spline test + consecutive-timepoint check):\n")
high_confidence_induced %>% count(condition) %>% print()


#How many of the Wurtzel genes are IDed as WI in each dataset?
Wurtzel.anno.rename<-Wurtzel.anno
colnames(Wurtzel.anno.rename)<-c("gene","gene.name.short","wurtzel.name","cell.type")
Wurtzel.merge<-left_join(Wurtzel.anno.rename,induction_results_wide,by="gene")

write.table(Wurtzel.merge,"Wurtzel.merge.names.RNAseq.wide.csv")

Wurtzel.merge %>% filter(F_adj_pvalue_Incision<0.05) %>% nrow()
Wurtzel.merge %>% filter(F_adj_pvalue_HeadRegen<0.05) %>% nrow()
Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05) %>% nrow()

Wurtzel.merge %>% filter(F_adj_pvalue_HeadRegen<0.05) %>% nrow()
Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05 & F_adj_pvalue_HeadRegen<0.05) %>% nrow()
Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05 & F_adj_pvalue_Incision<0.05) %>% nrow()

#library(VennDiagram)
library(ggVennDiagram)

my_data <- list(
  Set_1 = Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05 & F_adj_pvalue_HeadRegen<0.05) %>% pull(gene),
  Set_2 = Wurtzel.merge %>% filter(F_adj_pvalue_Incision<0.05) %>% pull(gene)
)

ggVennDiagram(my_data, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFF", high = "#0073C2") +
  theme(legend.position = "none") # Optional: hide the count legend


#install.packages("eulerr")
#install.packages("polylabelr")

library(eulerr)
library(polylabelr)

Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05 | F_adj_pvalue_HeadRegen<0.05) %>% nrow()
Wurtzel.merge %>% filter(F_adj_pvalue_Incision<0.05) %>% nrow()

frac.HeadORTail.induced.in.incision<-(Wurtzel.merge %>% filter(F_adj_pvalue_Incision<0.05) %>% nrow() %>% as.numeric())/(Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05 | F_adj_pvalue_HeadRegen<0.05) %>% nrow() %>% as.numeric())
print(frac.HeadORTail.induced.in.incision)

RNAseq_significant <- list(
  Head_OR_Tail = Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05 | F_adj_pvalue_HeadRegen<0.05) %>% pull(gene),
  Incision = Wurtzel.merge %>% filter(F_adj_pvalue_Incision<0.05) %>% pull(gene)
)

fit_euler <- euler(RNAseq_significant)


plot(fit_euler, 
     quantities = TRUE,           # Displays the counts in each region
     fills = c("skyblue", "pink"), # Custom circle colors
     edges = list(col = "black"))  # Border colors



plot_venn<-plot(fit_euler, 
                quantities = TRUE,           # Displays the counts in each region
                fills = c("skyblue", "pink"), # Custom circle colors
                edges = list(col = "black"))  # Border colors


png("Venn_Wurtzel_HeadTail_vs_Incision.png", width = 800, height = 600, res = 150)
print(plot(plot_venn, quantities=TRUE))
dev.off()

ggsave("Venn_Wurtzel_HeadTail_vs_Incision.png")


RNAseq_significant_Wurtzel_each <- list(
  Head = Wurtzel.merge %>% filter(F_adj_pvalue_HeadRegen<0.05) %>% filter(peak_log2FC_HeadRegen > 0) %>% pull(gene),
  Tail = Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05) %>% filter(peak_log2FC_TailRegen > 0) %>% pull(gene),
  Incision = Wurtzel.merge %>% filter(F_adj_pvalue_Incision<0.05) %>% filter(peak_log2FC_Incision > 0) %>% pull(gene)
)

fit_euler_Wurtzel_each <- euler(RNAseq_significant_Wurtzel_each)

png("euler_diagram_Wurtzel_all.png", width = 800, height = 800, res = 150)
plot_venn_Wurtzel_each<-plot(fit_euler_Wurtzel_each, 
                             quantities = TRUE,           # Displays the counts in each region
                             fills = c("Head"="red", "Tail"="#1E90FF", "Incision"="#4bb254"),
                             edges = list(col = "black"))  # Border colors

print(plot_venn_Wurtzel_each)
dev.off()


fit_euler_Wurtzel_each_classic <- venn(RNAseq_significant_Wurtzel_each)

plot_venn_Wurtzel_each_classic<-plot(fit_euler_Wurtzel_each_classic, 
                                     quantities = TRUE,           # Displays the counts in each region
                                     fills = c("Head"="red", "Tail"="#1E90FF", "Incision"="#4bb254"),
                                     edges = list(col = "black"))  # Border colors
print(plot_venn_Wurtzel_each_classic)



library(eulerr)
library(grid)

RNAseq_significant_Wurtzel_each <- list(
  Head = Wurtzel.merge %>% filter(F_adj_pvalue_HeadRegen < 0.05) %>% filter(peak_log2FC_HeadRegen > 0) %>% pull(gene),
  Tail = Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen < 0.05) %>% filter(peak_log2FC_TailRegen > 0) %>% pull(gene),
  Incision = Wurtzel.merge %>% filter(F_adj_pvalue_Incision < 0.05) %>% filter(peak_log2FC_Incision > 0) %>% pull(gene)
)

fit_euler_Wurtzel_each <- euler(RNAseq_significant_Wurtzel_each)

plot_venn_Wurtzel_each <- plot(fit_euler_Wurtzel_each,
                               quantities = list(fontsize = 16,fontface="bold"),
                               fills = c("Head" = "red", "Tail" = "#1E90FF", "Incision" = "#4bb254"),
                               edges = list(col = "black"),
                               labels = list(fontsize = 16, fontface="bold")
)

png("euler_diagram_Wurtzel_all.png", width = 1000, height = 1000, res = 150)
grid.newpage()
pushViewport(viewport(width = 0.85, height = 0.85))
grid.draw(plot_venn_Wurtzel_each)
popViewport()
dev.off()

##Figure S1B
library(eulerr)
library(polylabelr)

RNAseq_significant <- list(
  Head_OR_Tail = Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05 | F_adj_pvalue_HeadRegen<0.05) %>% pull(gene),
  Incision = Wurtzel.merge %>% filter(F_adj_pvalue_Incision<0.05) %>% pull(gene)
)

fit_euler <- euler(RNAseq_significant)

#"#4bb254", "Amputation" = "#1E90FF"
plot(fit_euler, 
     quantities = TRUE,           # Displays the counts in each region
     fills = c("#1E90FF", "#4bb254"), # Custom circle colors
     edges = list(col = "black"))  # Border colors

png("euler_diagram.png", width = 800, height = 600, res = 150)
venn_plot<-plot(fit_euler, 
                #     quantities = TRUE,           # Displays the counts in each region
                quantities = list(fontsize = 18), 
                fills = c("#1E90FF", "#4bb254"), # Custom circle colors
                #    fills = c("skyblue", "pink"), # Custom circle colors
                edges = list(col = "black"))  # Border colors
print(venn_plot)
dev.off()

##

#collect lists of genes induced in head or tail:

Wurtzel_induced_Head<- Wurtzel.merge %>% filter(F_adj_pvalue_HeadRegen<0.05) %>% filter(peak_log2FC_HeadRegen > 0) %>% pull(gene)
Wurtzel_induced_Tail<- Wurtzel.merge %>% filter(F_adj_pvalue_TailRegen<0.05) %>% filter(peak_log2FC_TailRegen > 0) %>% pull(gene)
Wurtzel_induced_Incision<-Wurtzel.merge %>% filter(F_adj_pvalue_Incision<0.05) %>% filter(peak_log2FC_Incision > 0) %>% pull(gene)
Wurtzel_induced_HorT<- Wurtzel.merge %>% filter(F_adj_pvalue_HeadRegen<0.05 | F_adj_pvalue_TailRegen<0.05) %>% filter(peak_log2FC_HeadRegen > 0) %>% pull(gene)
Wurtzel_induced_HorTnotIncision<-Wurtzel.merge %>% filter(F_adj_pvalue_HeadRegen<0.05 | F_adj_pvalue_TailRegen<0.05) %>% filter(peak_log2FC_HeadRegen > 0) %>% filter(!gene %in% Wurtzel_induced_Incision) %>% pull(gene)

write.table(as.data.frame(Wurtzel_induced_Head),"Wurtzel_induced_Head.csv")
write.table(as.data.frame(Wurtzel_induced_Tail),"Wurtzel_induced_Tail.csv")
write.table(as.data.frame(Wurtzel_induced_Incision),"Wurtzel_induced_Incision.csv")
write.table(as.data.frame(Wurtzel_induced_HorT),"Wurtzel_induced_HorT.csv")
write.table(as.data.frame(Wurtzel_induced_HorTnotIncision),"Wurtzel_induced_HorTnotIncision.csv")


