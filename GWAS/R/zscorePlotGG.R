#' Zscore plot
#'
#' This function plots and compares z-scores for each of the methods
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param sib Is a binary for indicating if we look at sibling history
#' @param IncludeTG Is a binary for indicating if we include True genetic liability
#' @keywords ZscorePlot
#' @export
#' @importFrom dplyr %>%
#' @examples
#' zscorePlotGG(total_indiv = 1000, SNP = 1000, h = 0.5, sib=0, includeTG = 0)
zscorePlotGG <- function(total_indiv, SNP, h, sib, includeTG){
  file <- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
  if (file.exists(file)){
    #Read relevant assoc-files
    file_stringer_LTFH <- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
    file_stringer_GWAX<- paste("./data/GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
    file_stringer_case_ctrl<- paste("./data/case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.qassoc", sep="")
    file_stringer_TG<- paste("./data/TG","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.qassoc", sep="")
    ltfh <- data.table::fread(file_stringer_LTFH)
    gwax <- data.table::fread(file_stringer_GWAX)
    case_ctrl <- data.table::fread(file_stringer_case_ctrl)
    TG <- data.table::fread(file_stringer_TG)
    #creating z-scores
    ltfh$Z <- sign(ltfh$BETA)*abs(qnorm(ltfh$P,lower.tail = FALSE))
    case_ctrl$Z <- sign(case_ctrl$BETA)*abs(qnorm(case_ctrl$P,lower.tail = FALSE))
    gwax$Z <- sign(gwax$BETA)*abs(qnorm(gwax$P,lower.tail = FALSE))
    TG$Z <- sign(TG$BETA)*abs(qnorm(TG$P,lower.tail = FALSE))
    # get true causals
    file_stringer_true <- paste("./data/BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.txt", sep="")
    true_beta <- data.table::fread(file_stringer_true)
    true_beta_causals <- true_beta %>%
      dplyr::filter(V2 != 0)
    #Creating subset of data based on true_causals (joined together)
    zscoredata <- true_beta_causals %>%
      dplyr::left_join(select(ltfh,c(SNP,Z)), by = c("V1" = "SNP")) %>%
      dplyr::rename("true_beta"="V2", "Z_ltfh" = "Z") %>%
      dplyr::left_join(select(gwax,c(SNP,Z)), by = c("V1" = "SNP")) %>%
      dplyr::rename("Z_gwax" = "Z") %>%
      dplyr::left_join(select(case_ctrl,c(SNP,Z)), by = c("V1" = "SNP")) %>%
      dplyr::rename("Z_case_ctrl" = "Z")%>%
      dplyr::left_join(select(TG,c(SNP,Z)), by = c("V1" = "SNP")) %>%
      dplyr::rename("Z_TG" = "Z")
    zscoredata <- zscoredata %>%
      dplyr::filter(abs(Z_ltfh) != Inf & abs(Z_gwax) != Inf & abs(Z_case_ctrl) != Inf & abs(Z_TG) != Inf)


    #PLOT OF LTFH VS GWAX
    Lm <- lm(Z_ltfh ~ 0 + Z_gwax, data = zscoredata)
    LmSum <- summary(Lm)
    p1 <- ggplot2::ggplot(data = zscoredata, ggplot2::aes(Z_ltfh, Z_gwax)) +
      ggplot2::geom_point(color="darkgrey") +
      ggplot2::geom_smooth(method = "lm", formula = y~0+x, col = "red") +

      ggplot2::ggtitle(paste("Z-Scores for LT-FH vs. GWAX")) +
      ggplot2::geom_text(x = 0, y = 9, label = paste("Squared slope = ",round(LmSum$coefficients[1]^2,3),sep=""), color = 'black') +

      ggplot2::labs(x=expression(paste(Z[GWAX])),y=expression(paste(Z[LT-FH]))) +

      ggplot2::theme_light() +

      ggplot2::theme(
        legend.position="none",
        panel.border = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank()
      )

    # PLOT OF THE LTFH VS CASE-CTRL
    Lm <- lm(Z_ltfh ~ 0 + Z_case_ctrl, data = zscoredata)
    LmSum <- summary(Lm)
    p2 <- ggplot2::ggplot(data = zscoredata, ggplot2::aes(Z_case_ctrl, Z_ltfh)) +
      ggplot2::geom_point(color="darkgrey") +
      ggplot2::geom_smooth(method = "lm", formula = y~0+x, col = "red") +

      ggplot2::ggtitle(paste("Z-Scores for LT-FH vs CASE-CTRL")) +
      ggplot2::geom_text(x = 0, y = 9, label = paste("Squared slope = ",round(LmSum$coefficients[1]^2,3),sep=""), color = 'black') +

      ggplot2::labs(x=expression(paste(Z[CASE_CTRL])),y=expression(paste(Z[LT-FH]))) +

      ggplot2::theme_light() +

      ggplot2::theme(
        legend.position="none",
        panel.border = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank()
      )

    # PLOT OF THE GWAX VS CASE-CTRL
    Lm <- lm(Z_gwax ~ 0 + Z_case_ctrl, data = zscoredata)
    LmSum <- summary(Lm)
    p3 <- ggplot2::ggplot(data = zscoredata, ggplot2::aes(Z_case_ctrl, Z_gwax)) +
      ggplot2::geom_point(color="darkgrey") +
      ggplot2::geom_smooth(method = "lm", formula = y~0+x, col = "red") +

      ggplot2::ggtitle(paste("Z-Scores for GWAX vs CASE-CTRL")) +
      ggplot2::geom_text(x = 0, y = 9, label = paste("Squared slope = ",round(LmSum$coefficients[1]^2,3),sep=""), color = 'black') +

      ggplot2::labs(x=expression(paste(Z[CASE_CTRL])),y=expression(paste(Z[GWAX]))) +

      ggplot2::theme_light() +

      ggplot2::theme(
        legend.position="none",
        panel.border = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank()
      )

    if (includeTG == 0){return(cowplot::plot_grid(p2,p1,p3,nrow=1))}
    else{
      #PLOT OF TRUE VS LTFH
      Lm <- lm(Z_TG ~ 0 + Z_ltfh, data = zscoredata)
      LmSum <- summary(Lm)
      p4 <- ggplot2::ggplot(data = zscoredata, ggplot2::aes(Z_TG, Z_ltfh)) +
        ggplot2::geom_point(color="darkgrey") +
        ggplot2::geom_smooth(method = "lm", formula = y~0+x, col = "red") +

        ggplot2::ggtitle(paste("Z-Scores for True Genetic Liability vs. LT-FH")) +
        ggplot2::geom_text(x = 0, y = 9, label = paste("Squared slope = ",round(LmSum$coefficients[1]^2,3),sep=""), color = 'black') +

        ggplot2::labs(x=expression(paste(Z[LT-FH])),y=expression(paste(Z[True]))) +

        ggplot2::theme_light() +

        ggplot2::theme(
          legend.position="none",
          panel.border = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank()
        )

      # PLOT OF THE TRUE VS GWAX
      Lm <- lm(Z_TG ~ 0 + Z_gwax, data = zscoredata)
      LmSum <- summary(Lm)
      p5 <- ggplot2::ggplot(data = zscoredata, ggplot2::aes(Z_TG, Z_gwax)) +
        ggplot2::geom_point(color="darkgrey") +
        ggplot2::geom_smooth(method = "lm", formula = y~0+x, col = "red") +

        ggplot2::ggtitle(paste("Z-Scores for True Genetic Liability vs GWAX")) +
        ggplot2::geom_text(x = 0, y = 9, label = paste("Squared slope = ",round(LmSum$coefficients[1]^2,3),sep=""), color = 'black') +

        ggplot2::labs(x=expression(paste(Z[GWAX])),y=expression(paste(Z[True]))) +

        ggplot2::theme_light() +

        ggplot2::theme(
          legend.position="none",
          panel.border = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank()
        )

      # PLOT OF THE TG VS CASE-CTRL
      Lm <- lm(Z_TG ~ 0 + Z_case_ctrl, data = zscoredata)
      LmSum <- summary(Lm)
      p6 <- ggplot2::ggplot(data = zscoredata, ggplot2::aes(Z_TG, Z_gwax)) +
        ggplot2::geom_point(color="darkgrey") +
        ggplot2::geom_smooth(method = "lm", formula = y~0+x, col = "red") +

        ggplot2::ggtitle(paste("Z-Scores for True Genetic Liability vs CASE-CTRL")) +
        ggplot2::geom_text(x = 0, y = 9, label = paste("Squared slope = ",round(LmSum$coefficients[1]^2,3),sep=""), color = 'black') +

        ggplot2::labs(x=expression(paste(Z[CASE_CTRL])),y=expression(paste(Z[True]))) +

        ggplot2::theme_light() +

        ggplot2::theme(
          legend.position="none",
          panel.border = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank()
        )
      return(cowplot::plot_grid(p2,p1,p3,p4,p5,p6,nrow=2))

    }
  } else {
    print("No data exists! - Try another predefined or run our simulation")
  }
}
