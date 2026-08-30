
table(combined_4disease_40GBM$L1_new)
unique(combined_4disease_40GBM$L1_new)

library(Seurat)
library(dplyr)
library(Seurat)
library(dplyr)
library(Matrix)

expr_data <- LayerData(
  combined_4disease_40GBM,
  assay = "RNA",
  layer = "data"
)

# 检查
dim(expr_data)


# ==============================
# 2. 获取细胞类型
# ==============================

celltype <- combined_4disease_40GBM$L1_new

celltype_order <- c(
  "Tumor",
  "Myeloid",
  "T/NK",
  "Fibroblast/Mural",
  "Glial",
  "Endothelial",
  "B cell"
)

celltype <- factor(
  celltype,
  levels = celltype_order
)

table(celltype)

# ==============================
# 3. 计算各 cell type 平均表达
# ==============================

avg_expr <- sapply(
  celltype_order,
  function(ct) {
    
    cells <- which(celltype == ct)
    
    Matrix::rowMeans(
      expr_data[, cells, drop = FALSE]
    )
  }
)

colnames(avg_expr) <- celltype_order

dim(avg_expr)
# ==============================
# 4. 计算相对其他细胞类型的表达优势
# ==============================

top_markers <- lapply(
  celltype_order,
  function(ct) {
    
    other_ct <- setdiff(
      celltype_order,
      ct
    )
    
    specificity <- avg_expr[, ct] -
      rowMeans(
        avg_expr[, other_ct, drop = FALSE]
      )
    
    result <- data.frame(
      gene = rownames(avg_expr),
      avg_expr = avg_expr[, ct],
      specificity = specificity
    )
    
    result %>%
      filter(
        avg_expr > 0,
        is.finite(specificity)
      ) %>%
      arrange(
        desc(specificity)
      ) %>%
      slice_head(n = 10) %>%
      mutate(
        celltype = ct
      )
  }
) %>%
  bind_rows()
top_markers %>%
  select(
    celltype,
    gene,
    avg_expr,
    specificity
  )





marker_genes <- c(
  # Tumor
  "SOX4", "SOX6", "PTPRZ1", "NAALADL2",
  
  # Myeloid
  "TYROBP", "C1QA", "C1QB", "C1QC",
  
  # T/NK
  "CD3E", "CD2", "CCL5", "IL32",
  
  # Fibroblast/Mural
  "COL1A1", "COL1A2", "COL3A1", "CALD1",
  
  # Glial
  "PLP1", "MBP", "PTGDS", "SPOCK1",
  
  # Endothelial
  "VWF", "FLT1", "PTPRG", "CALCRL",
  
  # B cell
  "IGKC", "IGHG1", "CD74", "MZB1"
)

marker_genes[!marker_genes %in% rownames(combined_4disease_40GBM)]














library(Seurat)
library(ggplot2)
library(dplyr)

# ==============================
# 1. 设置 marker genes
# ==============================

marker_genes <- c(
  # Tumor
  "SOX4", "SOX6", "PTPRZ1", "NAALADL2",
  
  # Myeloid
  "TYROBP", "C1QA", "C1QB", "C1QC",
  
  # T/NK
  "CD3E", "CD2", "CCL5", "IL32",
  
  # Fibroblast/Mural
  "COL1A1", "COL1A2", "COL3A1", "CALD1",
  
  # Glial
  "PLP1", "MBP", "PTGDS", "SPOCK1",
  
  # Endothelial
  "VWF", "FLT1", "PTPRG", "CALCRL",
  
  # B cell
  "IGKC", "IGHG1", "CD74", "MZB1"
)


# ==============================
# 2. 设置 cell type 顺序
# ==============================

celltype_order <- c(
  "Tumor",
  "Myeloid",
  "T/NK",
  "Fibroblast/Mural",
  "Glial",
  "Endothelial",
  "B cell"
)

combined_4disease_40GBM$L1_new <- factor(
  combined_4disease_40GBM$L1_new,
  levels = celltype_order
)

Idents(combined_4disease_40GBM) <- "L1_new"


# ==============================
# 3. 绘制 DotPlot
# ==============================
p_figure1C <- DotPlot(
  combined_4disease_40GBM,
  features = marker_genes,
  assay = "RNA",
  dot.scale = 6
) +
  
  # 颜色梯度
  
  scale_color_gradientn(
    colors = c(
      "#3F9BCB",
      "#FFFFFF",
      "#E56895"
    )
  ) +
  
  # 点大小
  scale_size(
    range = c(1, 7)
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    color = "Average expression",
    size = "Percent expressed"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 11,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    
    axis.title = element_blank(),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.border = element_blank(),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

p_figure1C
ggsave(
  "Figure_1C.pdf",
  plot = p_figure1C,
  width = 9.7,
  height = 4,
  units = "in"
)

