#' 质控指标整理
#'
#' @param data 样本质控结果表格，数据框，必填，具体格式参考QC数据
#' @param xlab 可视化结果x轴，通常为样本，字符型，必填
#' @param ylab 可视化结果y轴，通常为待检测的指标，字符型，必填
#' @param changelabname 是否需要更改x轴，逻辑型（T or F)，选填，默认为F
#' @param xlabname 可视化结果x轴的新名字，字符型，选填，默认为“”
#' @param ylabname 可视化结果y轴的新名字，字符型，选填，默认为“”
#' @param title 可视化结果是否需要标题，逻辑型（T or F)，选填，默认为F
#' @param titlename 可视化结果标题，字符型，选填，默认为“”
#' @param baralpha 可视化结果柱状图透明度，整型，选填，默认为0.5
#' @param barwidth 可视化结果柱状图宽度，整型，选填，默认为0.8
#' @param bar_pass_color 可视化结果柱状图中合格样本颜色，字符型，选填，默认为"#ACACAC"（灰色）
#' @param bar_fail_color 可视化结果柱状图中合格样本颜色，字符型，选填，默认为"#D32F2E"（红色）
#' @param facet 可视化结果是否需要分面，逻辑型（T or F)，选填，默认为F
#' @param facet_label 可视化结果分面依据，字符型，如果facet=T，必填（如"batch"）
#' @param showsampletype 可视化结果是否需要区分样本类型，逻辑型（T or F)，选填，默认为F
#' @param sampletype_label 可视化结果是否区分样本类型的依据，字符型，如果showsampletype=T，必填（如"type"）
#' @param showerror 可视化结果下是否需要整理整体结果，逻辑型（T or F)，选填，默认为F
#' @param showsample 可视化结果的x轴内容是否显示，逻辑型（T or F)，选填，默认为F
#' @param limit 指标阈值，整型，选填，默认为0，若不设置，所有样本均合格
#' @param limit_color 指标阈线值颜色，字符型，选填，默认为"red"
#' @param limit_size 指标阈线值宽度，整型，选填，默认为1
#' @param axisx_size 可视化结果x轴字体大小，整型，选填，默认为12
#' @param axisy_size 可视化结果y轴字体大小，整型，选填，默认为16
#' @param axis_title_size 可视化结果标题字体大小，整型，选填，默认为16
#' @param showerrordetail 结果中错误样本是否显示具体数值，逻辑型（T or F)，选填，默认为T
#' @param errordetail_color 结果中错误样本的具体数值颜色，字符型，选填，默认为"black"
#' @param numberdig 结果中错误样本的具体数值保留小数位数，字符型，整型，选填，默认为3
#' @param axisdetail_size 可视化结果y轴字体大小，整型，选填，默认为16
#' @import ggplot2
#' @import ggthemes
#' @import gridExtra
#' @import ggsci
#'
#' @return
#' @export
#'
#' @examples RawdataQC_standard(data=QC,xlab="library",ylab="gc_content",limit=45)
RawdataQC_standard <- function(data = data, xlab = xlab, ylab = ylab,
                                    changelabname = F, xlabname = "", ylabname = "",
                                    title = F, titlename = "",
                                    baralpha = 0.5, barwidth = 0.8,
                                    bar_pass_color = "#ACACAC", bar_fail_color = "#D32F2E",
                                    facet = F, facet_label = "batch",
                                    showsampletype = F, sampletype_label = "type",
                                    showerror = F, showsample = F,
                                    limit = 0, limit_color = "red", limit_size = 1,
                                    axisx_size = 12, axisy_size = 16,axisdetail_size=6,axis_title_size = 16,
                                    showerrordetail = T, errordetail_color = "black", numberdig = 3) {
  ## Data Prep
  data[, which(colnames(data) == xlab)] <- factor(data[, which(colnames(data) == xlab)],
    levels = c(data[, which(colnames(data) == xlab)]
    [order(data[, which(colnames(data) == ylab)], decreasing = T)])
  )
  data[, which(colnames(data) == ylab)] <- as.numeric(data[, which(colnames(data) == ylab)])
  data$pass <- "PASS"
  data$pass[data[, which(colnames(data) == ylab)] < limit] <- "FAIL"
  data$pass <- factor(data$pass, levels = c("PASS", "FAIL"))

  ## Table Generation
  SummaryTable_total_reads_faillist <- data.frame(
    Sample = data[data[, which(colnames(data) == ylab)] < limit, which(colnames(data) == xlab)],
    Value = round(data[data[, which(colnames(data) == ylab)] < limit, which(colnames(data) == ylab)], 3)
  )
  names(SummaryTable_total_reads_faillist) <- c(xlab, ylab)
  if (changelabname == T) {
    names(SummaryTable_total_reads_faillist) <- c(xlabname, ylabname)
  }
  data_pass <- as.data.frame(table(data$pass))

  if (length(data_pass$Var1) == 1) {
    if (unique(data_pass$Var1) == "FAIL") {
      SummaryTable_total_reads <- data.frame(Pass = 0, Fail = length(data$pass), Total = dim(data)[1], Pass_Ratio = 0.00)
    } else {
      SummaryTable_total_reads <- data.frame(
        Pass = length(data$pass), Fail = 0,
        Total = dim(data)[1], Pass_Ratio = 1.00
      )
    }
  } else {
    SummaryTable_total_reads <- data.frame(
      Pass = data_pass[which(data_pass$Var1 == "PASS"), which(colnames(data_pass) == "Freq")],
      Fail = data_pass[which(data_pass$Var1 == "FAIL"), which(colnames(data_pass) == "Freq")],
      Total = dim(data)[1],
      Pass_Ratio = round(data_pass[
        which(data_pass$Var1 == "PASS"),
        which(colnames(data_pass) == "Freq")
      ] / dim(data)[1], 2)
    )
  }
  names(SummaryTable_total_reads) <- c(
    "Number of passed samples", "Number of failed samples",
    "Number of total samples", "Ratio of passed samples"
  )

  ## Plot Generation
  ### Overall Plot
  if (showsampletype == T) {
    p <- ggplot(data) +
      geom_bar(aes(x = get(xlab), y = get(ylab), fill = get(sampletype_label)),
        stat = "identity", alpha = baralpha, width = barwidth
      )
    p <- p + scale_fill_nejm(name = "Type")
  } else {
    p <- ggplot(data, aes(x = get(xlab), y = get(ylab), fill = pass)) +
      geom_bar(stat = "identity", alpha = baralpha, width = barwidth) +
      scale_fill_manual(name = "Data Quality", values = c("PASS" = bar_pass_color, "FAIL" = bar_fail_color))
  }
  p <- p + theme_few() +
    geom_hline(yintercept = limit, linetype = "dotted", color = limit_color, size = limit_size) +
    theme(
      panel.border = element_blank(),
      strip.text.x = element_text(size = 20),
      panel.spacing.x = unit(0.2, "mm")
    ) +
    labs(x = xlab, y = ylab) +
    theme(
      axis.text.y = element_text(size = axisy_size, face = "plain", color = "black"),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = axis_title_size, face = "bold"),
      axis.title.y = element_text(size = axisy_size, face = "plain", color = "black"),
      axis.title.x = element_text(size = axisy_size, face = "plain", color = "black"),
      axis.ticks = element_blank()
    )

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
    data_f <- data[data[, which(colnames(data) == ylab)] < limit, ]
    if (showsampletype == T) {
      failplot <- ggplot(data_f) +
        geom_bar(aes(x = get(xlab), y = get(ylab), fill = get(sampletype_label)), stat = "identity", alpha = baralpha, width = barwidth)
      failplot <- failplot + scale_fill_nejm(name = "Type") +
        theme_few() +
        theme(
          panel.border = element_blank(),
          strip.text.x = element_text(size = 20),
          panel.spacing.x = unit(0.2, "mm")
        )
    } else {
      failplot <- ggplot(data_f, aes(x = get(xlab), y = get(ylab), fill = pass)) +
        geom_bar(stat = "identity", alpha = baralpha, width = barwidth) +
        scale_fill_manual(name = "Data Quality", values = c("PASS" = bar_pass_color, "FAIL" = bar_fail_color)) +
        theme_few() +
        theme(
          panel.border = element_blank(), legend.position = "none",
          strip.text.x = element_text(size = 20), panel.spacing.x = unit(0.2, "mm")
        )
    }
    failplot <- failplot + labs(x = xlab, y = ylab) +
      theme(
        axis.text.y = element_text(size = axisy_size, face = "plain", color = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = axisy_size, face = "plain", color = "black"),
        axis.title.x = element_text(size = axisy_size, face = "plain", color = "black"),
        plot.title = element_text(hjust = 0.5, size = axis_title_size, face = "bold"),
        axis.ticks = element_blank()
      ) +
      theme(axis.text.x = element_text(size = axisx_size, angle = 70, face = "plain", color = "black", vjust = 0.5, hjust = 0.5))
    if (title == T) {
      failplot <- failplot + ggtitle(paste(titlename, "(Failed Samples)"))
    }
    if (changelabname == T) {
      failplot <- failplot + xlab(xlabname) + ylab(ylabname)
    }
    if (facet == T) {
      failplot <- failplot + facet_grid(~ get(facet_label), scales = "free_x", space = "free_x")
    }
    if (showerrordetail == T) {
      failplot <- failplot + geom_text(aes(
        x = get(xlab), y = max(get(ylab)) * 1.1,
        label = sprintf(paste0("%.", numberdig, "f"), get(ylab))
      ),
      size = axisdetail_size,
      color = errordetail_color,
      position = position_dodge(width = 0.9), show.legend = F
      )
    }
  }
  return(list(
    plot = p, faillist = SummaryTable_total_reads_faillist,
    summary = SummaryTable_total_reads, failplot = failplot
  ))
}
