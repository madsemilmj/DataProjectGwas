#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(devtools)

library(shiny)
library(knitr)
library(markdown)
library(DT)

library(tidyverse)
library(knitr)
library(data.table)
library(vroom)
library(future) 
library(flock)
library(cowplot)
library(ggplot2)
library(plotly)
library(dplyr)
library(shinythemes)
library(GWAS)
library(purrr)
library(stats)
library(graphics)
library(usethis)
library(here)
#rmdfiles <- c("markdowns/skabelon_til_shiny.rmd", "markdowns/gibs_ltfh_markdown.rmd")
#sapply(rmdfiles, knit, quiet = T)




# Define UI for application that draws a histogram
shinyUI(  navbarPage(title = "GWAS", 
                     
  
                   theme = shinytheme("yeti"),

  navbarMenu(title = "DATA", 
             

             
        
      tabPanel("Welcome", 
               withMathJax(), tags$head(
                 tags$style(
                   HTML(
                     ".MathJax {
                       font-size: 2 em;
                     }"
                   )
                 )
               ),
               
                 
                 
                 
          fluidRow(column(12),
                   
          
  
                 
                img(src = "frontpage_dna.jpg", width="100%")
            
          
          ),
          hr(),
          fluidRow(column(10, offset= 1, 
          (helpText(h4("Welcome to our shiny app, this application is developed to help the user to get a better understanding of Genome Wide Association Studies.
          Before you begin to explore its implementations and visualizations, we encourage you to read the following introduction of the application. 
          The shiny-app is a product with which you can employ genetic data in order to determine significant gene mutations associated with a particular disease.
          Moreover, it is a tool that can compare the methods that go into such analysis."))), br(), 
          
          strong("Biological Terms"),br(), 
          helpText(h4("As we work with biological data, an understanding for several biological concepts are needed in order to gain a complete perspective on the mechanics of the application.", br(), 
          "In general, genotype is an organism's full hereditary information - that is, the passing on of traits from parents to their offspring; whereas phenotype is an organism's actual observed characteristic. 
          This distinction is important in understanding how diseases function in a biological context, entailing that genes contribute to a trait. Thus, the phenotype is the observable expressions of the genes.")),
          helpText(h4("Researchers take blood samples from each participant in each group and harvest DNA from these samples. Machinery enables researchers to survey each person’s complete set of DNA, or genome, from which 
          they strategically select markers of genetic variation. These markers are named single nucleotide polymorphisms, or SNP’s; for an individual, these traits are often referred to as an individual's genotype.")), 
          br(), 
          strong("Our Data"), br(), 
          helpText(h4("Single-Nucleotide Polymorphisms (SNPs) are base-pair values and can obtain the values \\(A,C,G,T\\). In this regard, SNPs are the deviation from a 'reference' genome and are the target of our analysis. If the majority of a population has \\(C\\) with a subgroup that has an \\(A\\), then \\(C\\) is the major allele and \\(A\\) is the minor allele.
          The minor allele is usually considered the ”risk allele”, which is also the case in this project.
          Thus, our data consists of counts on how many minor alleles are present in our population sample - in practice, this implies that a \\(0\\) denotes zero minor alleles, a \\(1\\) denotes one minor allele and \\(2\\) denotes two minor alleles present at a position \\(j\\).In practice, however, we normalize these with their respective frequencies such that they obtain zero mean and standard deviation one")),
          br(), 
          strong("Our simulation"),br(),
          helpText(h4("As we are not able to access real-life genetic data because of privacy protocols, we have simulated our data. This Shiny-app enables its user to simulate \\(N\\) individuals and \\(M\\) SNPs, with \\(C\\) causal ( \\(\\beta \\neq 0 \\))  SNPs. 
                      In the simulation, it is assumed that the SNPs follows a binomial distribution for position \\(1 \\leq j \\leq M\\):
                      $$x_j \\sim binom(2,MAF_j) $$ Here, MAF denotes the minor allele frequency which is the frequency of the allele with the fewest occurrences (considered to be the 'risk' allele). We assume that the MAF follows a uniform distribution for position 
                      \\(1 \\leq j \\leq M\\): $$MAF_j \\sim unif(0.01, 0.49)$$ The causal SNPs constitute only a subset of the whole set of SNPs and are the ones we assume to enjoy a causal relationship with the disease; ideally, the causal SNPs are the SNPs we wish to identify through our analysis. 
In our simulation, we assume that the effect sizes of the causal SNPs follow a normal distribution such that: $$\\beta_j \\sim N \\left(0, \\frac{h^2}{c} \\right)
$$ Here, \\(h^2 \\in (0,1)\\) denotes the heritability of a trait. With regard to the non-causal SNPs, or null SNPs, their effect sizes are imposed such that \\(\\beta=0\\)."
                      
                      
                      )))),
          fluidRow(column(7, offset = 1, 
          strong("Phenotype"), br(), 
          helpText(h4("When we assign phenotypes, we utilize a liability-threshold model in which an individual is a disease case if and only if an underlying continuous-valued liability lies above a threshold \\(T\\). In our case, this liability is a sum of the genetic liability and the enviromental liabilty
          such that: $$ l=l_g+l_e \\sim N(0,1), $$ where \\(l_g\\sim N(0,h^2)\\) and \\(l_e \\sim N(0,1-h^2) \\).
            In other words, we have a case \\((y=1)\\) if \\(l \\geq T\\) with \\(P(l \\geq T)=K\\); (\\(y=0\\), otherwise). Here, \\(T\\) denotes the threshold of disease liability, and \\(K\\) denotes  prevalence of trait."))
          ,br(),
          strong("Analysis"), br(),
          helpText(h4("This application is the result of a genome wide association study (GWAS). In general, GWAS-studies are observational studies of genetic data of different individuals to see if any gentic variations are associated with a trait. In order to explore these associations, our analysis centers on the linear regression: $$ y_i = x_{ij} \\cdot \\beta_j +e_i $$  Here, \\( y_i \\) is the phenotype for individual  \\(i \\), \\(x_{ij} \\) is the SNP for individual \\(i \\) at position \\(j\\) (we normalize these during the anlysis), \\(\\beta_j\\) is the effect size of the SNP at position \\(j\\) and \\(e_i \\) 
          is the individual error term - we assume that the intercept is zero in our analysis. We employ three different outcomes in our analysis: GWAS, GWAX and LT-FH. In short, GWAS uses binary case–control status, ignoring family history; GWAX uses binary proxy case–control status, merging controls who have a family history of disease with disease cases; and LT-FH uses continuous-valued posterior mean genetic liability, 
                      appropriately differentiating all case–control and family history configurations. Figure 1a shows how the different phenotypes vary: "))), column(3,  br(),  br(), img(src = "figure1_ltfh.jpg", width="100%"), br(), "Figure 1a 
                                                                                                                                                                        Hujoel et al., 2019, Combining case-control status and family history of disease increases association power. BioRxiv, 722645. https://doi.org/10.1101/722645")
          
          ), fluidRow(column(10, offset=1, 
          strong("Structure of Shiny App"),br(), 
          helpText(h4("The overall structure of the application is the six tabs where the user can explore the different methods employed in the analysis and the corresponding summary statistics for each method. The first tab 'data' is the section dedicated to the raw data. In this section, you can choose between pre-simulated data, simulate your own data, or upload your own data. Regardless of your choice, the application will perform the analysis and present the results in the ensuing tabs.", br(), br(),
                      
                      "In the tabs 'Case-Control', 'GWAX', 'LTFH', you can modify the analysis in relation to which of the three methods you wish to employ. The app will then present the result in an interactive data-table, enabling the user to click on whichever SNP from the data and receive the position in the Manhattan plot along with a plot presenting the effect size for the selected SNP. In the 'LTFH'-tab, you can also click on the 'Gibbs sampler'-subsection for information and code-elaboration on the Gibbs sampler-algorithm with which we have used to perform the LT-FH analysis. You can also click on the 'visualization'-subsection, where you will find a plot of the posterior mean genetic liabilities for the different family history configurations. ",
                      br(), br(), 
                      "If you click on the tab 'Comparing models', you will find a subsection 'Statistic strength' dedicated to comparing the results and strengths of the three different methods. Here, a \\(Z\\)-score plot displays the interrelated strength between the methods, whereas three boxplots compare ten simulations of the same data. A subsection 'Manhattan Plots' is also available in which you can explore an interactive Manhattan plot of all SNPs and their corresponding \\(p\\)-values and causality information. In conjunction with the Manhattan plot, a confusion matrix provides additional information on how many causal SNPs are tested as significant, etc. ",br(), br(),
"Finally, the 'about' tab presents information on the authors of the application along with a link to the source code." )),
                 
            
          
          br(), br(),br(), br(), br(), br()
          
          
          
          
          )
          
          
          
),#end of fluidrow
    
    
    
    
        
    
   
),


tabPanel("Choose Data", fluidRow( titlePanel(column(3,offset = 6, "Data"))),
    
          
         fluidRow(column(10, offset=1,   br(), br(), br(),
                         h4(HTML("<p>In this tab you can choose the data that goes in to the analysis in the app. <b>Be aware</b> that the choice you make on this page goes in to the entire analysis. If you choose to use the predifined data,
                                 the data is already simulated and you can freely choose between the different options of presimulated data in the analysis tabs. 
                                 If you decide to simulate your own data remember to press the button 'run simulation'. This choice might take a while, depending on the size of data you want to simulate. Lastly, you can also choose to upload your data, in this section it is very important to be very meticulous since the analysis only works with data that aligns perfectly with the format that  <a href='https://zzz.bwh.harvard.edu/plink/data.shtml'> Plink</a> uses.</p>")
                         
                         ))), 
         
         
         fluidRow(
    
    
    
    
    column(10, offset = 1 ,wellPanel(style = "background: #ddf2f3", 
      
        verbatimTextOutput("predifined"),
    radioButtons("datatype", label = h3("Choose data option"),
                 choices = c("Use Predifined Data" = 1, "simulate new data (This may take a while)" = 2, "Upload your own data"=3), 
                 selected = 1))),
          ),

    fluidRow(
        
    # Display this only if the density is shown
    conditionalPanel(condition = "input.datatype == 2", column(5, offset = 1,
                     sliderInput(inputId = "Individuals_sim",
                                 label = "Individuals:",
                                 min = 1000, max = 150000, value = 10000, step = 1000),
                     sliderInput(inputId = "SNPs_sim",
                                 label = "SNPs:",
                                 min = 1000, max = 150000, value = 10000, step = 1000)
    
    ),
    column(5,
                    sliderInput(inputId = "H2_sim",
                                label = "Heritability:",
                                min = 0.01, max = 0.99, value = 0.5, step = 0.01),
                    radioButtons(inputId = "sib_sim",
                                label = "Include sibling history:", choices = c("Yes"=1, "No"=0 )),
                    actionButton("run_sim", "Run simulation"))
   
                                         
                                         
                                         
    ), #conditionalpanel end
   
      
   conditionalPanel(condition = "input.datatype == 1", column(10, offset = 1,
                                                              h3("This choice is recomended - you can now just go through the different pages to investigate the different possibilities of GWAS"))
   ),   
    
    
    
   
   conditionalPanel(condition = "input.datatype == 3",  column(2, offset =2, wellPanel( 
                                                              
                                                              fileInput("upload", h4("PED-file")) , br(), textOutput("uploader") 
                                                      
                                                              
                                                              )),column(2,  wellPanel(
                                                                       fileInput("upload1", h4("Phenotype file")), br(), textOutput("uploader1") 
                                                                       
                                                                       
                                                              )),column(2,  wellPanel(
                                                                       
                                                                       fileInput("upload2", h4("True file")), br(), textOutput("uploader2") 
                                                                      
                                                              )), column(2, wellPanel(
                                                                        
                                                                        fileInput("upload3", h4("MAP-file")), br(), textOutput("uploader3") 
                                                              
                                                                      
                                                              )), 
                    
   ) 
   
  
   
                                                              
  
    
   ), 
    #fluid row end
   hr(),
fluidRow(conditionalPanel(condition= "input.run_sim", column(10, offset = 1,
                          
                          h3("Simulation has started, Done! will be printed when datasimulation is completed"),
                          
                          textOutput("simulator") ))
         
           ),#fluid row end     

         
         
), #tabpanel end 


tabPanel("Simulation",    fluidRow( titlePanel(column(4,offset = 5, "Simulation"))),
         hr(),
         fluidRow(column(10, offset = 1, 
                         helpText(h4(br(),br(),"First, as shown above, we set a seed such that our data is reproducable.",
                                     br(), br(),"This code presents how we simulated our data. The simulation is a prerequisite for our analysis, as we do not have acess to real-life data.", 
                                     br(), "Thus, the simulation is the groundwork for the entire project. Throughout the simulation, time performance has been the key of our optimization process. Thus, our coding choices reflect the fact we are working with big datasets.",
                                     br(), "Before we go into the actual simulation of data, we need to define some functions that become handy later.",
                                     br(), br(),"We use this function such that we are able to normalize our data:", 
                                     br(), br())),
                         
                         verbatimTextOutput("KODE1"), helpText(h4(br(),
                                                                  br(), 
                                                                  "This function simulate the siblings. In our project, an individual can have up to seven siblings. ", br(),
                                                                  br(), "We have looked at stats from families in Denmark here: https://www.dst.dk/da/Statistik/emner/befolkning-og-valg/husstande-familier-boern/boern",
                                                                  br(), br(),
                                                                  "From this, we have obtained the following distribution over families with up to seven siblings: \\(sib_b \\sim Binom(7,0.16) \\)", br(),
                                                                                          br(), br(),
                                                                   "We use this function to convert our denotation of the risk alleles (0,1,2) into (AA, AC, CC), accommodating PLINK's notation for risk alleles. C is the risk allele, i.e. CC = 2 and AA = 0.",
                                                                   br(), br())),
                         
                         verbatimTextOutput("KODE2"), helpText(h4( br(), br(),
                                                                   "We are able adjust our simulation of our data with regard to the number of individuals and SNPs we desire.
                       Further, we can also adjust the heritability, the number of causal SNPs, our liability-threshold and the number of workers with which we parallelize the simulation.",
                                                                   br(), br(),
                                                                   "As we are able to adjust the parameters, we can visualize the effects of such adjustments - look for the sections regarding these aspects.",
                                                                   br(), br())),
                         
                         verbatimTextOutput("KODE3"), helpText(h4(
                           br(), br(),
                           "Our T denotes the threshold of disease liablity such that an individual is a case \\(y=1\\) if \\(l \\geq T \\) with \\(P(l \\geq T)=k \\).", br(),
                           br(), "Thus, k denotes the prevalence of a trait, meaning that a k percentage of a population is assumed to have the trait.", br(), br())),
                         
                         verbatimTextOutput("KODE4"), helpText(h4( br(), br(),
                                                                   "Minor allele frequency (MAF) is the frequency of the allele with the fewest occurences (this allele is also denoted as the risk allele).
                       The risk allele is the allele which we assumed to have a relationsip with the disease.", br(), br(),
                                                                   "We assume how frequent the risk alleles are at each position j (SNP) and we assume that risk alleles are uniformly distributed such that \\( MAF_j\\sim unif(0.01,049) \\). The maximum value of 0.49 implies that we are making assumptions about the minorly frequent allele.", br(), br())),
                         
                         verbatimTextOutput("KODE5"), helpText(h4( br(),
                                                                   "We create a MAP-file that contains four columns which are: \\(\\textit{chromosome}\\), (in our case assigned to 1 for all rows); \\(\\textit{rs#}\\) (which were just  numbered from 1:100.000); \\(\\textit{Genetic distance}\\) (which is a column of zeros) and; \\(\\textit{Base-pair position}\\) (which is also numbered from 1:100.000).
                       The MAP-file and the MAF-file is created with fwrite such that it is readable by PLINK. ", br())),
                         
                         verbatimTextOutput("KODE6"), helpText(h4( br(), br(),
                                                                   "Now, we simulate our causal SNPs. These are the SNPs which are assumed to have a true causal effect. 
                       These causal SNPs is just a fraction of the total number of SNPs. All non-causal SNPs have an effect size of zero.", br(), br(),
                                                                   "The effect sizes of the causal SNPs are simulated as such:\\(  \\beta_j \\sim N(0,\\frac{h^2}{c}) \\)", br(), 
                                                                   "Here, \\( h \\in(0,1) \\) and is the heritability of a trait.", br(), br())),
                         
                         verbatimTextOutput("KODE7"), helpText(h4( br(), br(),
                                                                   "We simulate the genotypes of mom and dad. We simulate the genotypes such that \\( x_j \\sim binom(2,MAF_j) \\)", br(), br())),
                         
                         verbatimTextOutput("KODE8"), helpText(h4( br(), br(), 
                                                                   "Here, for each SNP for the child, we take the sum of both parents and divide it by two, implying that a child's genotype is the mean of its parents' genotypes.", br(), br(),
                                                                   "Further, we simulate the genotypes of seven siblings in exactly the same way as the child.", br(), 
                                                                   "It is assumed that the full liability is the sum of the enviromental liability and the genetic liability.", br(), br(),
                                                                   "In our simulation we have that: \\(l=l_g+l_e \\sim N(0,1) \\), where \\( l_g \\sim N(0,h^2) \\).", br(), br())),
                         
                         verbatimTextOutput("KODE9"), helpText(h4( br(), br(),
                                                                    "We normalize our data with our defined function.", br(), 
                                                                    "We overwrite the files as they preceeding files are no longer needed.", br(), br())),
                         
                         verbatimTextOutput("KODE10"), helpText(h4( br(), br(),
                                                                    "We employ the causal SNPs on each child such that the effect sizes of the causal SNPs are imposed if and only if a child has that particular risk allele; otherwise, it is zero.", br(), br())),
                         
                         verbatimTextOutput("KODE11"), helpText(h4( br(), br(), 
                                                                    "We assume that the environmental liablity follows a normal distribution: \\( l_e \\sim N(0,1-h^2) \\).", br(), br())),
                         
                         verbatimTextOutput("KODE12"), helpText(h4( br(), br(),
                                                                    "As stated above, the full liability is the sum of the genetic liability and the enviromental liability.", br(),
                                                                    "Further, we create the phenotype for each child.", br(), br())),
                         
                         verbatimTextOutput("KODE13"), helpText(h4( br(), br(), 
                                                                    "Following the same process, we create the phenotypes for both mom and dad:", br(), br())),
                         
                         verbatimTextOutput("KODE14"), helpText(h4( br(), br(), 
                                                                    "Now, we bind all of our values together such that we include all availabe information for each individual. 
                       This also includes information on siblings. 
                       We then write txt-files for our data such that we can employ PLINK.
                       PLINK has a number pre-fixed formats, and these are all met in our data.", br(),
                                                                    "The pre-fixed formats are all found here: https://zzz.bwh.harvard.edu/plink/", br(), br(),
                                                                    "The first chunk of data is now written to their respective txt.files:", br(), br() )),
                         
                         verbatimTextOutput("KODE15"), helpText(h4(br(), br(), 
                                                                   "The abobe code takes care of the first chunk of individuals. Now, we only need the rest of the individuals.", br(), br(),
                                                                   "For this, we parallelize our simulation. This implies that our calculations are computed simultaneously such that we can execute our code more effectively with higher performance. 
                      This is particularly important as we simulate huge amounts of data, which is highly time consuming.", br(), br(),
                                                                   "If we did not parallelize our simulation, we would not be able to execute the code, as we do not have access to suffcient processing power for the gravity of such computation. Hence, parallelization is a necessity.", br(), br(),
                                                                   "When parallelizing, it is important that we make sure that we never write on the same data. Therefore, we initiate by making a temporary file to make sure we are able to write to a file simultaneous when using parallel computing.
                       The future_apply/future_lapply are apply/lapply in the context of the 'future' framework, which is an R package that supports parallelization. 
                       As mentioned before, it is paramount that each worker does not work on (simulate) the same individual. 
                       Therefore, the 'id_fra'  and 'id_til' make sure that we are always working with only an individual at the particular chunk of individuals. 'id_fra' is the first individual in the chunk for each iteration. 
                       Subsuquently, 'id_til' is the last individual in the chunk for each iteration.", br(), br(),
                                                                   "Future.seed is set to true such that the data is reproducable.", br(), br())),
                         
                         verbatimTextOutput("KODE16"), helpText(h4( br(), br(),
                                                                    
                         )
                         ))
                  
                  
                  
                  
                  
                  
                  
         ))

         
         
         

           ), #DATA end tab end
     
# CASECTRLCASECASECTRLCTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRLCASECTRL


    
    
    tabPanel("Case-Control", 
             
             fluidRow( titlePanel(column(10, "GWAS with Case/Control"))),
             hr(),
             fluidRow(column(12, HTML("<p>This section presents the summary statistics for the regression run by <a href='https://zzz.bwh.harvard.edu/plink/data.shtml'> Plink</a> in the context of case-control. The analysis is run directly from R using linear regression by calling the command prompt:</p>"),   verbatimTextOutput("GWASplink"), "This implies that we have a binary phenotype as our dependent variable, indicating whether or not the offspring is a case. (See GWAS in figure 1a in the 'welcome tab').
           In the data table, you can click on whichever SNP you desire and obtain the p-value; even further, a manhattan plot appears from which you can observe all SNPs and their respective LOD-scores (a LOD-score: \\(-\\log_{10}(P-value)\\)) along with an arrow that indicates which SNP you clicked on. Additionally, a plot with the estimated beta coefficient, or the effect size, for the selected SNP, holding everything else constant will appear
The slope of the graph indicates whether the selected SNP is estimated to have a positive or negative effect on the phenotype.
In practice, we narrow our analysis down to a single SNP \\(j\\):$$y_i=x_{ij} \\cdot \\beta_j+ e_i$$
Here, \\(y_i\\) is the phenotype for individual \\(i\\), \\(x_{ij}\\) is the SNP for individual \\(i\\) at position \\(j\\), \\(\\beta_j\\) is the effect size of the SNP at posistion \\(j\\), and \\(e_i\\) is the individual error term.
Moreover, the true \\(\\beta\\) (causal/non-causal) for the selected SNP is added along with SE-bands so that we can observe the interplay of these elements.")),
             br(),
             fluidRow(conditionalPanel(condition = "input.datatype == 1",
             column(4, wellPanel(br(), br(),
               radioButtons(inputId = "Individuals_cc",        label = h3("Individuals"),
                            choices = c("10.000" = format(as.numeric(10000), scientific = F), "50.000" = as.numeric(50000), "100.000"= format(as.numeric(100000),scientific = F)), 
                            selected = format(as.numeric(100000),scientific = F), inline = T),br(),
               radioButtons("SNPs_cc",        label = h3("SNPs"),
                            choices = c("50.000" = as.numeric(50000), "75.000" = 75000, "100.000"=format(as.numeric(100000),scientific = F)), 
                            selected = format(as.numeric(100000),scientific = F) ,inline = T), br(),
               radioButtons("h2_cc",        label = h3("liability"),
                            choices = c("0.25" = as.numeric(25), "0.5" = as.numeric(50), "0.75"=as.numeric(75)), 
                            selected = 50 ,inline = T), br(), br(),br() )),
              
                                    
               column(8, DTOutput("casectrl"))),
             
             conditionalPanel(condition = "input.datatype == 2", 
                              column(4, wellPanel("You have  chosen to simulate the below data configuration, if no datatable appears on the right side of the screen it means you forgot to press 'run simulation'. You can change this in 'Data' \\(\\rightarrow \\) 'Choose data'", tableOutput("table_of_sim"))),
                              
                              
                              column(8, DTOutput("casectrlsim"))   ),
             conditionalPanel(condition = "input.datatype == 3", column(10, offset = 1, wellPanel( h3("You have chosen to upload your own data, yet you did not upload a valid dataset. Please return to 'Data' \\(\\rightarrow \\) 'Choose data' to upload valid data or choose a different data configuration." 
             
             ))))),
             
             fluidRow(conditionalPanel(condition = "input.datatype == 1",column(6, plotOutput("manhattanCC")),
                                       column(6, plotOutput("RegressionCC"))),
                      conditionalPanel(condition = "input.datatype == 2",column(6, plotOutput("manhattanCCsim")),
                                       column(6, plotOutput("RegressionCCsim"))),
                      
                      
             br(), br(),br(), br(),br(), br(),br(), br() ) # fluid row end 
             
            
             ), #tabpanel end

    

          
          
# GWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAXGWAX          

          

           tabPanel("GWAX", 
                    
                    
                    fluidRow( titlePanel(column(10, "Genome Wide Association with family history as a proXy"))),
                    hr(),
                    fluidRow(column(12, HTML("<p>This section presents the summary statistics for the regression run by <a href='https://zzz.bwh.harvard.edu/plink/data.shtml'> Plink</a> in the context of GWAX. The analysis is run directly from R using linear regression by calling the command prompt:</p>"),   verbatimTextOutput("GWAXplink"),"This implies that we have a binary phenotype as our dependent variable, indicating whether or not the offspring or any family member is a case. (See GWAX in figure 1a in the 'welcome tab').
In the data table, you can click on whichever SNP you desire and obtain the p-value; even further, a manhattan plot appears from which you can observe all SNPs and their respective LOD-scores (a LOD-score: \\(-\\log_{10}(P-value)\\)) along with an arrow that indicates which SNP you clicked on. Additionally, a plot with the estimated beta coefficient, or the effect size, for the selected SNP, holding everything else constant will appear.
The slope of the graph indicates whether the selected SNP is estimated to have a positive or negative effect on the phenotype.
In practice, we narrow our analysis down to a single SNP \\(j\\):$$y_i=x_{ij} \\cdot \\beta_j+e_i$$
Here, \\(y_i\\) is the phenotype for individual \\(i\\), \\(x_{ij}\\) is the SNP for individual \\(i\\)  at position \\(j\\), \\(\\beta_j\\) is the effect size of the SNP at posistion \\(j\\), and \\(e_i\\) is the individual error term.
Moreover, the true \\(\\beta\\) (causal/non-causal) for the selected SNP is added along with SE-bands so that we can observe the interplay of these elements.")),
                    br(),
                    fluidRow(
                    conditionalPanel(condition = "input.datatype == 1",
                    column(4, wellPanel(
                      radioButtons(inputId = "Individuals_GW", h3("Individuals"), c("10.000" = format(as.numeric(10000), scientific = F), "50.000" = as.numeric(50000), "100.000"= format(as.numeric(100000),scientific = F)), 
                                   selected = format(as.numeric(100000),scientific = F),inline = T),
                      radioButtons("SNPs_GW",        label = h3("SNPs"),
                                   choices = c("50.000" = 50000, "75.000" = 75000, "100.000"=format(as.numeric(100000),scientific = F)), 
                                   selected = format(as.numeric(100000),scientific = F),inline = T),
                      radioButtons("h2_GW",        label = h3("liability"),
                                   choices = c("0.25" = 25, "0.5" = 50, "0.75"=75), 
                                   selected = 50, inline = T),
                      radioButtons("Sib_GW",        label = h3("Include sibling history"),
                                   choices = c("Yes" = 1, "No" = 0), 
                                   selected = 1 ,inline = T))),
                      column(8, DTOutput("GWAXhead")
                      
                    )),conditionalPanel(condition = "input.datatype==2", 
                    column(4, wellPanel("You have  chosen to simulate the below data configuration, if no datatable appears on the right side of the screen it means you forgot to press 'run simulation'. You can change this in 'Data' \\(\\rightarrow \\) 'Choose data'", 
                           
                           tableOutput("table_of_sim1"), radioButtons("Sib_GW_sim",        label = h4("Include sibling history"),
                                                                      choices = c("Yes" = 1, "No" = 0), 
                                                                      selected = 1,inline = T))), 
                    column(8, DTOutput("GWAXheadsim"))
                    ),
                    conditionalPanel(condition = "input.datatype == 3", column(10, offset = 1, wellPanel( h3("You have chosen to upload your own data, yet you did not upload a valid dataset. Please return to 'Data' \\(\\rightarrow \\) 'Choose data' to upload valid data or choose a different data configuration." 
                                                                                                             
                    ))))),
                      
             
                    
                    fluidRow(conditionalPanel(condition = "input.datatype == 1",column(6, plotOutput("manhattanGWAX")),
                                              column(6, plotOutput("RegressionGWAX"))),
                             conditionalPanel(condition = "input.datatype == 2",column(6, plotOutput("manhattanGWAXsim")),
                                              column(6, plotOutput("RegressionGWAXsim")))
                             
                             
                    )
                    

                     
                    ), #tabpanel end
                    
                    
                    
                    
      
                    
                    
                   
                    
                    
                    
                    
                    
                  
           
      #navbarMenu end

          
# LTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFHLTFH
          
navbarMenu(title = "LTFH", 
           tabPanel("Gibbs Sampling", withMathJax(), titlePanel(fluidRow(column(8, offset =4,  "Convergence Test Gibbs Sampling"))), 
                    fluidRow(column(10, offset = 1,
                                    helpText(h3("How and why do we use Gibbs Sampling")),
                                    helpText(h4('\\(\\textbf{Gibbs sampling}\\) is a sampling algorithm designed to draw samples from \\( \\textbf{posterior distributions} \\)
                    (however it can be used to draw samples from any distribution like in our case the multivariate normal) 
                    when sampling from those distributions are not straight-forward. When we know the conditional distribution 
                    of each of the parameters from the target distribution (posterior distribution), and hence 
                    it is possible to draw samples from them - we use the Gibbs sampling algorithm (compared to for instance other 
                    more general algorithms like \\(\\textbf{Metropolis Hastings}\\)).')), br(),
                                    helpText(h4('The algorithm reduces the difficult problem of sampling from a multivariate distribution 
                             down to a sequence of simpler univariate problems.Therefore in order to use the Gibbs Sampling 
                             algorithm we need to know (and be able to draw samples from) the full conditional distributions for 
                             each parameter.')), br(),
                                    helpText(h4("In our case, we have the following assumption about the distribution of liabilities: $$
                                  \\left(\\begin{matrix}
                                  \\epsilon_{o,e} \\\\
                                  \\epsilon_{o,g} \\\\
                                  \\epsilon_{p_1} \\\\
                                  \\epsilon_{p_2} \\\\
                                  \\epsilon_s
                                  \\end{matrix}\\right) \\sim MVN_5\\left(\\left(\\begin{matrix}
                                                                                                                0 \\\\
                                                                                                                0 \\\\
                                                                                                                0 \\\\
                                                                                                                0 \\\\
                                                                                                                0
                                                                                                                \\end{matrix}\\right), \\left(\\begin{matrix}
                                               1-h^2 & 0 & 0 & 0 & 0 \\\\
                                               0 & h^2 & \\frac{h^2}{2} & \\frac{h^2}{2} & \\frac{h^2}{2} \\\\
                                               0 & \\frac{h^2}{2} & 1 & 0 & \\frac{h^2}{2} \\\\
                                               0 & \\frac{h^2}{2} & 0 & 1 & \\frac{h^2}{2} \\\\
                                               0 & \\frac{h^2}{2} & \\frac{h^2}{2} & \\frac{h^2}{2} & 1 \\\\
                                               \\end{matrix}\\right)\\right)
                                  $$")), br(),
                                    helpText(h4('Where: \\(\\epsilon_{o,e} \\) is the environmental liability for the offspring, \\(\\epsilon_{o,g}\\) is the genetic liability for the offspring, while \\(\\left(\\epsilon_{p_1},\\epsilon_{p_2},\\epsilon_{s} \\right) \\) are the full liability for the mom, dad and sibling(s) respectively.')), br(),
                                    helpText(h4("Since we in our case are interested in the posterior mean genetic liability 
                                \\(E[\\epsilon_{o,g}|z_o,z_{p_1},z_{p_2},\\overset{\\rightarrow}{z_s}]\\) 
                                for each individual given the case-control status of the genotyped individual 
                                \\((z_o)\\), both parents \\((z_{p_1}, z_{p_2})\\), and sibling(s) \\((\\overset{\\rightarrow}{z_s})\\)
                                - where an individual is a case \\((z=1)\\) if and only if \\(\\epsilon \\geq T\\), otherwise we have \\((z=0)\\) - 
                                we are interested in drawing samples from a truncated version of a linear transformation of the above multivariate normal, 
                                where the bounds depends on the specific configuration. The linear transformation of multivariate normal distribution is defined by the matrix \\(A\\) satisfying: 
                                $$
                                A\\left(\\begin{matrix}
                                  \\epsilon_{o,e} \\\\
                                  \\epsilon_{o,g} \\\\
                                  \\epsilon_{p_1} \\\\
                                  \\epsilon_{p_2} \\\\
                                  \\epsilon_s
                                  \\end{matrix}\\right) = \\left(\\begin{matrix}
                                  \\epsilon_{o,g} \\\\
                                  \\epsilon_{o} \\\\
                                  \\epsilon_{p_1} \\\\
                                  \\epsilon_{p_2} \\\\
                                  \\epsilon_s
                                  \\end{matrix}\\right)
                                $$
                                where \\(\\epsilon_o = \\epsilon_{o,e}+\\epsilon_{o,g} \\), and because the linear transformed random variable then follows a \\(MVN_5(A\\boldsymbol{\\mu}, A\\Sigma A^T)\\), where \\(\\boldsymbol{\\mu}\\) 
                                and \\(\\Sigma\\) is defined as above we get:
                                
                                $$
                                \\left(\\begin{matrix}
                                  \\epsilon_{o,g} \\\\
                                  \\epsilon_{o} \\\\
                                  \\epsilon_{p_1} \\\\
                                  \\epsilon_{p_2} \\\\
                                  \\epsilon_s
                                  \\end{matrix}\\right) \\sim MVN_5\\left(\\left(\\begin{matrix}
                                                                                                                0 \\\\
                                                                                                                0 \\\\
                                                                                                                0 \\\\
                                                                                                                0 \\\\
                                                                                                                0
                                                                                                                \\end{matrix}\\right), \\left(\\begin{matrix}
h^2 & h^2 & \\frac{h^2}{2} & \\frac{h^2}{2} & \\frac{h^2}{2} \\\\
h^2 & 1 & \\frac{h^2}{2} & \\frac{h^2}{2} & \\frac{h^2}{2} \\\\
\\frac{h^2}{2} & \\frac{h^2}{2} & 1 & 0 & \\frac{h^2}{2} \\\\
\\frac{h^2}{2} & \\frac{h^2}{2} & 0 & 1 & \\frac{h^2}{2} \\\\
\\frac{h^2}{2} & \\frac{h^2}{2} & \\frac{h^2}{2} & \\frac{h^2}{2} & 1
\\end{matrix}\\right)\\right)
                                $$
                                
                                
                                Choosing to draw from a truncated multivariate 
                                normal optimize the time taken to estimate the posterior mean genetic liability 
                                \\(E[\\epsilon_{o,g}|z_o,z_{p_1},z_{p_2},\\overset{\\rightarrow}{z_s}]\\), since we would need a 
                                s.e.m (Standard Error of the Mean) for the estimated paramter to be less than \\(0.01\\) 
                                compared with the alternative of just drawing samples from the multivariate normal 
                                (where we would need a huge amount of samples such that s.e.m is less than \\(0.01\\) for all configurations). 
                                The fact that we are interested in drawing samples from the truncated multivariate normal fits nicely with 
                                the Gibbs-sampling algorithm where we just draw samples from the truncated versions of the  conditional 
                                distributions.")), br(),
                                    
                                    helpText(h4("In this case we know the univariate distribution of the full conditional distribution for each parameter, 
                                which will be normal distributions and can be calculated using:")),
                                    
                    )),
                    fluidRow(column(4,offset=4,
                                    img(src = "ConDist.PNG", width="100%"), br(), "                                       Source= Wikipedia"
                    )), 
                    fluidRow(column(10, offset = 1,
                                    
                                    helpText(h4("Once we are able to draw samples from the conditional distributions, the general Gibbs Sampling algorithm is:"))
                                    
                    )),
                    fluidRow(column(4,offset=4,
                                    img(src = "GibsAlgo.PNG", width="100%"), br(), "                                       Source= Bayesian Statistical Methods, Reich and Ghosh"
                    )), 
                    fluidRow(column(10, offset = 1,
                                    helpText(h4("The way we implemented the Gibbs-Sampling algorithm, where we repeatedly draw samples from the full 
                                conditional distributions (truncated normal distributions) is by the use of both the Cummulative Distriution 
                                Function (CDF) and the inverse CDF (quantile function) for the normal distribution, as well as the Uniform 
                                distribution.")), br(),
                                    helpText(h4("We first define the lowerBound \\((a)\\) and upperBound \\((b)\\) for the truncated distribution, then we use the 
                                CDF of the full conditional distributions (Normal Distributions) to find the probabilities \\(P(X\\leq a)\\) and 
                                \\(P(X\\leq b)\\), where \\(X\\sim Normal(\\mu,\\sigma^2)\\). (see below figure) This gives an interval of probabilites 
                                such that \\(p = F^{-1}(z)\\) where \\(z\\in (P(X\\leq a),P(X\\leq b))\\) satisfies \\(a\\leq p \\leq b\\). We then draw a 
                                single sample/value from the Uniform distribution with minimum \\(=P(X\\leq a)\\) and maximum \\(=P(X\\leq b)\\), and 
                                the value of that sample is then converted using the invers CDF function \\(F^{-1}\\) for the normal distribution 
                                to give us the actual sample.")), br())),
                                    fluidRow(column(4,offset=4,
                                    img(src = "cdfbillede.PNG", width="100%"), br(), "Source= https://www.r-bloggers.com/2020/08/generating-data-from-a-truncated-distribution/?fbclid=IwAR2GaH-Dk0mEfwG0Pu0cua6x7vqyfgKD64efzPpN5ZIU9I7Pdz3lHQqtCbQ"
                                    )), 
                    fluidRow(column(10, offset = 1,
                                    helpText(h4("Here is our implementation of the Gibbs Sampling algorithm. In the below code, we do not specify any specific
                                upper and lower bound (they are set to \\(\\infty\\) and \\(-\\infty\\) respectively), which is used to draw samples 
                                from the truncated version of the multivariate normal. In the subsequent sections we will explore different 
                                configurations resulting in different truncated multivariate normal distributions.")), br(),
                                    verbatimTextOutput("gibbscode1"), br(),
                                    helpText(h4("We will now look at different configurations as well as different initial values and explore the 
                                behaviour/convergence of our samples. In each configuration and stat-values we will plot traceplots 
                                for each of the parameters in the multivariate normal, as well as our defined burn-in used in our 
                                analysis (red vertical line). We will also plot two of the parameters (that are supposed to be normal 
                                distributed) with their respective bounds (red vertical and horizontal lines) with contour-lines for 
                                their correspoding bivariate normal-distribution superimposed.")), br(),
                                    helpText(h3(strong("All family-members is a case with start values at \\(10\\)"))), br(),
                                    img(src = "Gibsplot1.svg", width="100%"), br(),
                                    helpText(h3(strong("No family-members is a case with start values at \\(10\\)"))), br(),
                                    img(src = "Gibsplot2.svg", width="100%")
                                    
                                    
                                    
                                    
                    )),
                    
                    
                    
           
            
                    
                    
                    
                    
           ),
           
           tabPanel("Method", 
                    fluidRow( titlePanel(column(10, "Posterior Plots"))),
                    hr(),
                    fluidRow(column(12, "This plot visualizes the",  em("posterior mean genetic liabilities for all case-control and family history configurations"),".", span( "The blue plot", style = "color:blue"), "displays a \\(Normal(0,h^2)\\) distribution (initial genetic liability).", span("The red plot", style = "color:red"), "is the distribution of the configuration, with the dotted line representing the posterior mean genetic liability. 
                             The user is encouraged to adjust the parameters such that it becomes evident how the posterior mean genetic liability is affected by the different configurations.")),  br(),
                    strong("Be attentive:"), "the plot does not appear if the given configuration does, ", strong("NOT"), "appear in the chosen data.",
                    
                    
                    hr() ,
                    
                    fluidRow(column(4, wellPanel(
                      checkboxGroupInput(inputId = "Plotsterior",        label = h3("Choose case history"),
                                         choices= c("Child is a case", "Mother is a case","Father is a case" ), choiceValues = c(1,2,3)), 
                      conditionalPanel(condition = "input.Sib_ltfh == 1", numericInput("Sibs", "Siblings:", 2, min = 0, max = 7),
                                       numericInput("Sibs_case", "Siblings that are cases", 1, min=0, max=7),
                      ))),  
                      column(8, plotOutput("plotsteriot_pre"))
                      
                      
                    )
           ),  # tabpanel     
           tabPanel("Summary Statistics",
                    fluidRow( titlePanel(column(10, "Liability Threshold Family History"))),
                    hr(),
                    fluidRow(column(12, HTML("<p>This section presents the summary statistics for the regression run by <a href='https://zzz.bwh.harvard.edu/plink/data.shtml'> Plink</a> in the context of LTFH.The analysis is run directly from R using linear regression by calling the command prompt:</p>"),   verbatimTextOutput("LTFHplink"),"Here the association statistics are computed via regression of genotypes and posterior mean genetic liabilities. (See LTFH in figure 1a in the 'welcome tab').
In the data table, you can click on whichever SNP you desire and obtain the p-value; even further, a manhattan plot appears from which you can observe all SNPs and their respective LOD-scores (a LOD-score: (a LOD-score: \\(-\\log_{10}(P-value)\\)) along with an arrow 
that indicates which SNP you clicked on. Additionally, a plot with the estimated beta coefficient, or the effect size, for the selected SNP, holding everything else constant will appear
The slope of the graph indicates whether the selected SNP is estimated to have a positive or negative effect on the phenotype.
In practice, we narrow our analysis down to a single SNP \\(j\\):$$y_i=x_{ij} \\cdot \\beta_j+e_i$$
Here, \\(y_i\\) is the phenotype for individual \\(i\\), \\(x_{ij}\\) is the SNP for individual \\(i\\) at position \\(j\\), \\(\\beta_j\\) is the effect size of the SNP at posistion \\(j\\) and \\(e_i\\) is the individual error term.
Moreover, the true \\(\\beta\\) (causal/non-causal) for the selected SNP is added along with SE-bands so that we can observe the interplay of these elements.
                                    
                                    ")),
                    br(),
                    fluidRow(conditionalPanel(condition = "input.datatype == 1",
                    
                    
                    
                    column(4, wellPanel(
                      radioButtons(inputId = "Individuals_lfth",        label = h3("Individuals"),
                                   choices = c("10.000" = 10000 , "50.000" = 50000, "100.000"= format(as.numeric(100000),scientific = F)), 
                                   selected = format(as.numeric(100000),scientific = F) ,inline = T ),
                      radioButtons("SNPs_ltfh",        label = h3("SNPs"),
                                   choices = c("50.000" = 50000, "75.000" = 75000, "100.000"=format(as.numeric(100000),scientific = F)), 
                                   selected = format(as.numeric(100000),scientific = F),inline = T),
                      radioButtons("h2_ltfh",        label = h3("liability"),
                                   choices = c("0.25" = 25, "0.5" = 50, "0.75"=75), 
                                   selected = 50,inline = T),
                      radioButtons("Sib_ltfh",        label = h3("Include sibling history"),
                                   choices = c("Yes" = 1, "No" = 0), 
                                   selected = 1,inline = T)
                      
                    )),
                    column(8, DTOutput("LTFHhead"))), 
                    
                    conditionalPanel(condition = "input.datatype == 2",   
                                     column(4, wellPanel("You have  chosen to simulate the below data configuration, if no datatable appears on the right side of the screen it means you forgot to press 'run simulation'. You can change this in 'Data' \\(\\rightarrow \\) 'Choose data'",
                                       tableOutput("table_of_sim2"), 
                                       radioButtons("Sib_ltfh_sim", label = h4("Include sibling history"), choices = c("Yes" = 1, "No" = 0), selected = 1,inline = T))), 
                      column(8, DTOutput("LTFHheadsim"))
                                     
                    
              
                    ),
                    conditionalPanel(condition = "input.datatype == 3", column(10, offset = 1, wellPanel( h3("You have chosen to upload your own data, yet you did not upload a valid dataset. Please return to 'Data' \\(\\rightarrow \\) 'Choose data' to upload valid data or choose a different data configuration." 
                                                                                                             
                    ))))),
                    
                    
                    
                    fluidRow(conditionalPanel(condition = "input.datatype == 1",column(6, plotOutput("manhattanLTFH")),
                             column(6, plotOutput("Regressionltfh"))),
                             conditionalPanel(condition = "input.datatype == 2",column(6, plotOutput("manhattanLTFHsim")),
                                              column(6, plotOutput("Regressionltfhsim")))
                            
                             
                    ) #fluidrow end
                    
                    ) #tabpanel end 
                    
               
                    
                    
                    
                    
                    

           
     
           
           
           
           
           
), #navbarMenu end

#comparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodelscomparingmodels
          
          navbarMenu(title = "Comparing models", 
                     
                     tabPanel("Statistic Strength",
                              fluidRow( titlePanel(column(10, " Statistic Strength "))),
                              hr(),
                              fluidRow(column(12, helpText(h4("The value of the \\(Z\\)-score represents how many standard deviations a value is away from the mean. In the plot,
                                              the \\(Z\\)-statistics are plotted against each other for the different methods, where the intercept is zero and the slope is free.
                                              The slope squared is an expression of the strength of one method versus another. Below we plot each of the methods against eachother, 
                                                              on the left you have the option to choose between different data-configurations, as well as the option to include sibling history
                                                              (only relevant for GWAX and LT-FH) or if you want to include the method having the true genetic liabilities as phenotypes." )),
                                        br(),
                                        helpText(h4("The \\(\\textbf{box plot}\\) below first compares the average \\(\\chi^2\\) for both the true null SNP's and the causals SNPs for \\(10\\) simulations 
                                          all having the same configuration, with \\(100000\\) individuals, \\(100000\\) SNPs, liability of \\(0.5\\) and with a prevalence of \\(0.05\\).
                                          This plot will not change based on your choise on the left menu. This plot is used both as a check for everything is how it should be (we espect
                                          the mean of the \\( \\chi^2\\) to be around \\(1\\) under the null, but also as a way to compare the strength of the methods.
                                          The term \\(\\textit{power}\\) is calculated as the number of causal SNPs found in the analysis as a fraction of the total causal SNPs in the data. 
                                          The data simulation all contained \\(1000\\) causal SNPs.")), br(),
                                        helpText(h4(br(),strong("Disclaimer"), "These plots take a while to load."))
                                        
                                        
                                        
                              )),
                              br(),
                              fluidRow(conditionalPanel(condition = "input.datatype == 1",
                                column(3, wellPanel(
                                radioButtons(inputId = "Individuals_co",        label = h3("Individuals"),
                                             choices = c("10.000" = 10000 , "50.000" = 50000, "100.000"= format(as.numeric(100000),scientific = F)), 
                                             selected = format(as.numeric(100000),scientific = F) ,inline = T ),
                                radioButtons("SNPs_co",        label = h3("SNPs"),
                                             choices = c("50.000" = 50000, "75.000" = 75000, "100.000"=format(as.numeric(100000),scientific = F)), 
                                             selected = format(as.numeric(100000),scientific = F),inline = T),
                                radioButtons("h2_co",        label = h3("liability"),
                                             choices = c("0.25" = 25, "0.5" = 50, "0.75"= 75), 
                                             selected = 50 ,inline = T),
                                radioButtons("Sib_box1",        label = h3("Include sibling history"),
                                             choices = c("Yes" = format(as.numeric(1),scientific = F), "No" = format(as.numeric(0),scientific = F)), 
                                             selected = format(as.numeric(1),scientific = F) ,inline = T),
                                radioButtons("Truegen",   label= h3("Inlcude true genentic liability"), choices =c("Yes" = format(as.numeric(1),scientific = F), "No" = format(as.numeric(0),scientific = F)), 
                                selected = format(as.numeric(0),scientific = F) ,inline = T)
                                
                              )
                              
                              ), column(9, plotOutput("zscoreGG")),
                              
                              column(9, offset = 15,  plotOutput("boxcompare"))),
                              
                              conditionalPanel(condition = "input.datatype == 2",  
                                               column(3, wellPanel( "You have  chosen to simulate the below data configuration, if no plots appears after 15 seconds it means you forgot to press 'run simulation'. You can change this in 'Data' \\(\\rightarrow \\) 'Choose data'", tableOutput("table_of_sim3"), 
                                                            radioButtons("Sib_z_sim",        label = h4("Include sibling history"),
                                                                           choices = c("Yes" = 1, "No" = 0),                      selected = 1,inline = T), 
                                radioButtons("Truegen1",   label= h3("Inlcude true genentic liability"), choices =c("Yes" = format(as.numeric(1), scientific = F), "No" = format(as.numeric(0),scientific = F)), 
                                             selected = format(as.numeric(1),scientific = F) ,inline = T))), 
                                column(9, plotOutput("zscoreGG_sim")),
                                
                                column(9, offset = 15,  plotOutput("boxcompare1"))),
                              
                              
                              
                              conditionalPanel(condition = "input.datatype == 3", column(10, offset = 1, wellPanel( h3("You have chosen to upload your own data, yet you did not upload a valid dataset. Please return to 'Data' \\(\\rightarrow \\) 'Choose data' to upload valid data or choose a different data configuration." 
                                                                                                                       
                              ))))),
                              
                                
                                
                                
                                
                                
                               #fluidrow end
                            fluidRow(
                              br(), br(),br(), br(), br(), br(),br(), br(),br(), br()
                            )
                              
                              
                              
                     ),
                     
                     tabPanel("Manhattan Plots",
                              titlePanel("Manhattan Plots"),
                              fluidRow(column(7, helpText(h4("The Manhattan plot is a big scatter plot of the LOD-score 
                                                             against the SNP-number. The LOD score is the negative logarithm 
                                                             with base ten of the \\(p\\)-value; thus, the higher the LOD-score, the 
                                                             lower the \\(p\\)-value. All SNPs with LOD-scores under \\(1\\), corresponding to all \\(p\\)-values 
                                                             over \\(0.1\\), are excluded from the plot for time complexity reasons. If the user holds the 
                                                             mouse over a given SNP, then the plot provides information on the LOD-score for each SNP 
                                                             along with the true causality status."), br(), 
                                                          helpText(h4(span("The blue SNPs",style = "color:cornflowerblue"),"are the correctly predicted causal SNPs,", 
                                                                      span("the orange SNPs",style = "color:orange"), "are the causal SNPs we did not predict,", 
                                                                      span("the grey SNPs", style ="color:darkgrey"),"are the correctly predicted not-causal SNPs, 
                                                                                                                        and", span("the red SNPs", style="color:red"), "are not-causal SNPs the analysis predicted as causal. 
                                                                                                                        Further, a confusion matrix provides information on how many causal SNPs that the analysis 
                                                                                                                        detected under the selected method. The blue SNPs are the detected causal SNPs."))),
                                              helpText(h4(br(), br(),strong("Disclaimer"), "These plots take a while to load."))
                                              
                              ), column(5,  conditionalPanel(condition = "input.datatype == 1",  plotOutput("ConfussionM")) ,
                                                        conditionalPanel(condition = "input.datatype == 2",  plotOutput("ConfussionMsim")) )  ),
                              
                              fluidRow( column(12,offset = 5,  radioButtons(inputId = "Choose_type",        label = h3("Choose manhattan plot:"),
                                                          choices= c("Case-control"=1, "GWAX"=2,"LTFH"=3, "True Genetic liability"=4 ), selected = 1, inline = T))) ,
                              
                              fluidRow(conditionalPanel(condition = "input.datatype == 1",
                                       column(3, wellPanel(
                                radioButtons(inputId = "Individuals_ma",        label = h3("Individuals"),
                                             choices = c("10,000" = as.numeric(10000) , "50,000" = as.numeric(50000), "100,000"= format(as.numeric(100000),scientific = F)), 
                                             selected = format(as.numeric(100000),scientific = F) ,inline = T ),
                                radioButtons("SNPs_ma",        label = h3("SNPs"),
                                             choices = c("50,000" = as.numeric(50000), "75,000" = as.numeric(75000), "100,000"=format(as.numeric(100000),scientific = F)), 
                                             selected = format(as.numeric(100000),scientific = F),inline = T),
                                radioButtons("h2_ma",        label = h3("liability"),
                                             choices = c("0.25" = 25, "0.5" = 50, "0.75"= 75), 
                                             selected = 50 ,inline = T), 
                                radioButtons("Sib_box", label = h3("Include sibling history"), choices = c("Yes" = 1, "No" = 0), selected = 1,inline = T),
                                radioButtons(inputId = "Bonferronicoefficient",    label = h3("Bonferroni coefficient:"), choices = c("1" = as.numeric(1), "5,000" = as.numeric(5000), "100,000"= format(as.numeric(100000),scientific = F), "1,000,000"= format(as.numeric(1000000),scientific = F)), 
                                             selected = format(as.numeric(100000),scientific = F) ,inline = T)
                                
                              ))
                                
                                
                              
                              , 
                              
                                       column(9, plotlyOutput("bigmanhattancc")) 
                              
                      
                              ), 
                              conditionalPanel(condition = "input.datatype == 2",  
                                               column(3, wellPanel("You have  chosen to simulate the below data configuration, if no plot appears after a 15 seconds it means you forgot to press 'run simulation'. You can change this in 'Data' \\(\\rightarrow \\) 'Choose data'",  tableOutput("table_of_sim4"), 
                                                                     radioButtons("Sib_ma_sim",        label = h4("Include sibling history"),
                                                                                  choices = c("Yes" = 1, "No" = 0),                      selected = 1,inline = T), 
                                                                     radioButtons(inputId = "Bonferronicoefficientsim",    label = h3("Bonferroni coefficient:"), choices = c("1" = 1, "5,000" = 5000, "100,000"= format(as.numeric(100000),scientific = F), "1,000,000"= format(as.numeric(1000000),scientific = F)), 
                                                                                  selected = format(as.numeric(100000),scientific = F) ,inline = T)
                              
                                                                     
                                               
                                                                     )),
                                               
                                               
                                               column(9, plotlyOutput("bigmanhattansim"))
                                               ),
                              
                              conditionalPanel(condition = "input.datatype == 3", column(10, offset = 1, wellPanel( h3("You have chosen to upload your own data, yet you did not upload a valid dataset. Please return to 'Data' \\(\\rightarrow \\) 'Choose data' to upload valid data or choose a different data configuration." 
                                                                                                                       
                              )))),
                              
                              
                              br(), br(), br(), br(), br(), br(),  br(), br(), br() )
                            
                              
                              
                              )
                     
                     

                     
          ), #navbarMenu end
          

tabPanel("About us",                     fluidRow( titlePanel(column(4, offset = 6, "About us"))),
         hr(),
         
         fluidRow(column(10, offset=1, h4("The authors of this app are", strong("Jeppe Dalgaard Rode, Mads Emil Marker Jungersen,"), "and", strong("Thomas Lykke Rasmussen."), 
         "We are all students currently enrolled in the Data Science program at the department of mathematics at Aarhus University. 
         This Shiny app is created in connection with a dataproject-course durring our fourth semester where we had to work on a project with relation to Bayesian statistics. 
         Because we did not have access to real-life data, we have simulated all data such that we were able to employ methods such as GWAS. 
         This app serves to visualize our results in the exam, but hopefully it can also bring some value to other students who wish to engage in a similar project. 
         We have in the process made the package", code("GWAS"),  "which can be accesed through devtools via link:"),
         
         br(),
         h4(HTML("<p> By Clickling here <a href='https://github.com/madsemilmj/DataProjectGwas/tree/main/GWAS/R'> Github GWAS</a>.</p>")), 
         br(),
         h4("Or in R by loading: ", br())),
         fluidRow(column(10, offset=1, wellPanel(code("library(devtools)"), br(), 
                                       code("install_github('madsemilmj/DataProjectGwas/GWAS')"), br(),
                                            code("library(GWAS)")
                                   )))),
         
         
         fluidRow(column( 10, offset=1,strong(h3("Packages")), 
           verbatimTextOutput("sessionInfo")
                         )),

        fluidRow(column(10, offset=1, h4("If one were to experience any issus please contact us on email:"), br(), br(),  )),
        fluidRow(br(), column(3, offset = 1, img(src = "billede_Jeppe.png", width="50%"),br(), hr(),
                              h4("Jeppe Dalgaard Rode"), br(), 
                              h4("E-mail: jeppe.dalgaardrode@gmail.com"),
                              h4("Phone: +45 53 31 31 00")),
                 column(3, offset = 1, img(src = "billede_MadsEmil.png", width="50%"), br(), hr(),
                        h4("Mads Emil Marker Jungersen"), br(), 
                        h4("E-mail:  madsemil.mj@gmail.com"),
                        h4("Phone: +45 28 45 60 86")
                        
                        ), 
                 column(3,  offset = 1, img(src = "billede_Thomas.png", width="50%"), br(), hr(),
                        h4("Thomas Lykke Rasmussen"), br(), 
                        h4("E-mail:  thomaslykke@outlook.com"),
                        h4("Phone: +45 26 17 28 28", br(), br(), br(), br(),br(), br(), br(),br())
                        
                        ),
                               
                               
                )
                        
                        
                 )
         
         
 

    


)

    
)

