#' 质控指标整理（比对区域）
#' 
#' @param data 样本质控结果表格，数据框，必填，具体格式参考QC数据
#' @param xlab 可视化结果x轴，通常为样本，字符型，必填
#' @param var1 堆叠图第一个变量名称，字符型，选填，默认为“exonic”
#' @param var2 堆叠图第二个变量名称，字符型，选填，默认为“intronic”
#' @param var3 堆叠图第三个变量名称，字符型，选填，默认为“intergenic”
#' @param limit_var1 堆叠图第一个变量的阈值，整型，选填，默认为0.4
#' @param limit_var2 堆叠图第二个变量的阈值，整型，选填，默认为0.3
#' @param limit_var3 堆叠图第三个变量的阈值，整型，选填，默认为0
#' @param changelabname 是否需要更改x轴，逻辑型（T or F)，选填，默认为F
#' @param xlabname 可视化结果x轴的新名字，字符型，选填，默认为“”
#' @param ylabname 可视化结果y轴的新名字，字符型，选填，默认为“”
#' @param var1labname 堆叠图第一个变量在表格中的命名，字符型，选填，默认为"Exonic rate"
#' @param var2labname 堆叠图第二个变量在表格中的命名，字符型，选填，默认为"Intronic rate"
#' @param var3labname 堆叠图第三个变量在表格中的命名，字符型，选填，默认为"Intronic rate"
#' @param var3_color 堆叠图第一个变量在图中的颜色，字符型，选填，默认为"#66C2A5"
#' @param var2_color 堆叠图第二个变量在图中的颜色，字符型，选填，默认为"#FC8D62"
#' @param var1_color 堆叠图第三个变量在图中的颜色，字符型，选填，默认为“#8DA0CA”
#' @param title 可视化结果是否需要标题，逻辑型（T or F)，选填，默认为F
#' @param titlename 可视化结果标题，字符型，选填，默认为“”
#' @param baralpha 可视化结果柱状图透明度，整型，选填，默认为0.5
#' @param barwidth 可视化结果柱状图宽度，整型，选填，默认为0.8
#' @param bar_pass_color 可视化结果柱状图中合格样本颜色，字符型，选填，默认为"#ACACAC"（灰色）
#' @param bar_fail_color 可视化结果柱状图中合格样本颜色，字符型，选填，默认为"#D32F2E"（红色）
#' @param facet 可视化结果是否需要分面，逻辑型（T or F)，选填，默认为F
#' @param facet_label 可视化结果分面依据，字符型，如果facet=T，必填（如"batch"）
#' @param showerror 可视化结果下是否需要整理整体结果，逻辑型（T or F)，选填，默认为F
#' @param showsample 可视化结果的x轴内容是否显示，逻辑型（T or F)，选填，默认为F
#' @param axisx_size 可视化结果x轴字体大小，整型，选填，默认为12
#' @param axisy_size 可视化结果y轴字体大小，整型，选填，默认为16
#' @param axis_title_size 可视化结果标题字体大小，整型，选填，默认为16
#' @import ggplot2
#' @import ggthemes
#' @import gridExtra
#' @import ggsci
#' @return
#' @export
#'
#' @examples RawdataQC_stack(data=QC,xlab="library")


RawdataQC_stack <- function(data = data, xlab = xlab,
                                 var1 = "exonic", var2 = "intronic", var3 = "intergenic",
                                 limit_var1 = 0.40, limit_var2 = 0.3, limit_var3 = 0.00,
                                 changelabname = F, xlabname = xlabname, ylabname = ylabname,
                                 var1labname = "Exonic rate", var2labname = "Intronic rate", var3labname = "Intergenic rate",
                                 var3_color = "#66C2A5", var2_color = "#FC8D62", var1_color = "#8DA0CA",
                                 title = F, titlename = titlename,
                                 baralpha = 0.5, barwidth = 0.8,
                                 bar_pass_color = "#ACACAC", bar_fail_color = "#D32F2E",
                                 facet = F, facet_label = "batch",
                                 showerror = F, showsample = F,
                                 axisx_size = 12, axisy_size = 16, axis_title_size = 16) {
  ## Data Prep
  data$pass1 <- "PASS"
  data$pass1[data[, which(colnames(data) == var1)] < limit_var1] <- "FAIL"
  data$pass1 <- factor(data$pass1, levels = c("PASS", "FAIL"))
  data$pass2 <- "PASS"
  data$pass2[data[, which(colnames(data) == var2)] < limit_var2] <- "FAIL"
  data$pass2 <- factor(data$pass2, levels = c("PASS", "FAIL"))
  data$pass3 <- "PASS"
  data$pass3[data[, which(colnames(data) == var3)] < limit_var3] <- "FAIL"
  data$pass3 <- factor(data$pass3, levels = c("PASS", "FAIL"))
  data$pass <- paste(data$pass1, data$pass2, data$pass3, sep = "_")
  data$pass[-which(data$pass == "PASS_PASS_PASS")] <- "FAIL"
  data$pass[which(data$pass == "PASS_PASS_PASS")] <- "PASS"

  ## Table Generation
  SummaryTable_total_reads_faillist <- data.frame(
    Sample = data[which(data$pass == "FAIL"), which(colnames(data) == xlab)],
    Value_Var1 = round(data[which(data$pass == "FAIL"), which(colnames(data) == var1)], 3),
    Value_Var2 = round(data[which(data$pass == "FAIL"), which(colnames(data) == var2)], 3),
    Value_Var3 = round(data[which(data$pass == "FAIL"), which(colnames(data) == var3)], 3)
  )
  names(SummaryTable_total_reads_faillist) <- c(xlab, var1, var2, var3)
  if (changelabname == T) {
    names(SummaryTable_total_reads_faillist) <- c(xlabname, var1labname, var2labname, var3labname)
  }
  data_pass <- as.data.frame(table(data$pass))
  if (length(data_pass$Var1) == 1) {
    if (unique(data_pass$Var1) == "FAIL") {
      SummaryTable_total_reads <- data.frame(Pass = 0, Fail = length(data$pass), Total = dim(data)[1], Pass_Ratio = 0.00)
    } else {
      SummaryTable_total_reads <- data.frame(Pass = length(data$pass), Fail = 0, Total = dim(data)[1], Pass_Ratio = 1.00)
    }
  } else {
    SummaryTable_total_reads <- data.frame(
      Pass = data_pass[which(data_pass$Var1 == "PASS"), which(colnames(data_pass) == "Freq")],
      Fail = data_pass[which(data_pass$Var1 == "FAIL"), which(colnames(data_pass) == "Freq")],
      Total = dim(data)[1],
      Pass_Ratio = round(data_pass[which(data_pass$Var1 == "PASS"), which(colnames(data_pass) == "Freq")] / dim(data)[1], 2)
    )
  }
  names(SummaryTable_total_reads) <- c(
    "Number of passed samples", "Number of failed samples",
    "Number of total samples", "Ratio of passed samples"
  )

  ## Plot Generation
  if (facet == T) {
    data_melt <- data[, c(xlab, facet_label, var1, var2, var3)]
    data_melt <- reshape2::melt(data_melt, id.vars = c(xlab, facet_label), variable.name = "Region")
  } else {
    data_melt <- data[, c(xlab, var1, var2, var3)]
    data_melt <- reshape2::melt(data_melt, id.vars = c(xlab), variable.name = "Region")
  }
  data_melt$Region <- factor(data_melt$Region, levels = c(var3, var2, var1))
  data_melt$library <- factor(data_melt$library, levels = c(data_melt$library
  [order(data_melt[which(data_melt$Region == var1), "value"])]))

  ### Overall Plot
  p <- ggplot(data_melt, aes(x = get(xlab), y = value, fill = Region)) +
    geom_bar(stat = "identity", alpha = baralpha, width = barwidth) +
    theme_few() +
    theme(
      panel.border = element_blank(),
      strip.text.x = element_text(size = 20),
      panel.spacing.x = unit(0.2, "mm")
    ) +
    labs(x = xlab) +
    theme(
      axis.text.y = element_text(size = axisy_size, face = "plain", color = "black"),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = axis_title_size, face = "bold"),
      axis.title.y = element_text(size = axisy_size, face = "plain", color = "black"),
      axis.title.x = element_text(size = axisy_size, face = "plain", color = "black"),
      axis.ticks = element_blank()
    ) +
    scale_fill_manual(values = c(var3_color, var2_color, var1_color))
  if (title == T) {
    p <- p + ggtitle(titlename)
  }
  if (changelabname == T) {
    p <- p + xlab(xlabname) + ylab(ylabname)
  }
  if (showsample == T) {
    p <- p + theme(axis.text.x = element_text(size = axisx_size, angle = 90, face = "plain", color = "black"))
  }
  if (facet == T) {
    p <- p + facet_grid(~ get(facet_label), scales = "free_x", space = "free_x")
  }
  if (showerror == T) {
    tt <- ttheme_default(colhead = list(fg_params = list(parse = TRUE)))
    tbl_1 <- tableGrob(SummaryTable_total_reads, rows = NULL, theme = tt)
    p <- grid.arrange(p, tbl_1, nrow = 2, as.table = TRUE, heights = c(8, 1))
  }

  ### Faillist Plot
  if (dim(SummaryTable_total_reads_faillist)[1] == 0) {
    failplot <- print("There is no fail samples!")
  } else {
    data_f_1 <- data[data[, which(colnames(data) == xlab)] %in% SummaryTable_total_reads_faillist[, 1], ]
    if (facet == T) {
      data_f_2 <- data_f_1[, c(xlab, facet_label, var1, var2, var3)]
      data_f <- reshape2::melt(data_f_2, id.vars = c(xlab, facet_label), variable.name = "Region")
    } else {
      data_f_2 <- data_f_1[, c(xlab, var1, var2, var3)]
      data_f <- reshape2::melt(data_f_2, id.vars = c(xlab), variable.name = "Region")
    }
    data_f$Region <- factor(data_f$Region, levels = c(var3, var2, var1))
    data_f$library <- factor(data_f$library,
      levels = c(data_f$library[order(data_f[which(data_f$Region == var1), "value"])])
    )
    failplot <- ggplot(data_f, aes(x = get(xlab), y = value, fill = Region)) +
      geom_bar(stat = "identity", alpha = baralpha, width = barwidth) +
      theme_few() +
      theme(panel.border = element_blank(), strip.text.x = element_text(size = 20), panel.spacing.x = unit(0.2, "mm")) +
      labs(x = xlab) +
      theme(
        axis.text.y = element_text(size = axisy_size, face = "plain", color = "black"),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = axis_title_size, face = "bold"),
        axis.title.y = element_text(size = axisy_size, face = "plain", color = "black"),
        axis.title.x = element_text(size = axisy_size, face = "plain", color = "black"),
        axis.ticks = element_blank()
      ) +
      scale_fill_manual(values = c(var3_color, var2_color, var1_color))
    if (title == T) {
      failplot <- failplot + ggtitle(titlename)
    }
    if (changelabname == T) {
      failplot <- failplot + xlab(xlabname) + ylab(ylabname)
    }
    if (showsample == T) {
      failplot <- failplot +
        theme(axis.text.x = element_text(size = axisx_size, angle = 70, face = "plain", color = "black", vjust = 0.5, hjust = 0.5))
    }
    if (facet == T) {
      failplot <- failplot + facet_grid(~ get(facet_label), scales = "free_x")
    }
  }
  return(list(
    plot = p,
    faillist = SummaryTable_total_reads_faillist,
    summary = SummaryTable_total_reads,
    failplot = failplot
  ))
}
