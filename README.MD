# m6APrediction

Predict m6A modification sites from RNA sequences.

## Install

```r
# install.packages("remotes")
remotes::install_github("Susiekawaiidesu/m6APrediction")


## Usage

```r
library(m6APrediction)

# 单条预测（示例值）
prediction_single(
  seq = "AGCUAGCAUGCAUGCAUGCAUGC",
  gc_content = 0.52,
  RNA_type = "mRNA",
  RNA_region = "3UTR",
  exon_length = 1200,
  distance_to_junction = 150,
  evolutionary_conservation = 0.85,
  DNA_5mer = dna_encoding("AGCUAGCAUGCAUGCAUGC")
)

# 批量预测（data.frame 提前算好列）
df <- data.frame(
  id = c("seq1", "seq2"),
  seq = c("AGCUAGCAUGCAUGC", "GGGGGGGGGGGGGG"),
  gc_content = c(0.52, 0.48),
  RNA_type = c("mRNA", "lncRNA"),
  RNA_region = c("3UTR", "5UTR"),
  exon_length = c(1200, 800),
  distance_to_junction = c(150, 300),
  evolutionary_conservation = c(0.85, 0.22),
  DNA_5mer = I(list(dna_encoding("AGCUAGCAUGCAUGC"),
                    dna_encoding("GGGGGGGGGGGGGG")))
)
prediction_multiple(df)
