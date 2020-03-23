#' A function to simulate bulk or single cell RNA sequencing data
#'
#' @description This function simulates (bulk/single cell) RNA-seq dataset from semi-parametrically estimated distributions of gene expression levels in a given real source RNA-seq dataset
#'
#' @param n.sim an integer for the number of simulations to be generated
#' @param s.data a real source dataset (a SingleCellExperiment class or a matrix/data.frame of counts with genes in rows and samples in columns)
#' @param batch NULL or a vector containing batch indicators for each sample/cell in the source data
#' @param group NULL or a vector containing group indicators for each sample/cell in the source data
#' @param n.genes a numeric value for the total number of genes to be simulated
#' @param cand.DE.genes a list object contatining canidiate null and non-null (DE/predictor) genes.
#' If NULL (the default), an internal function determines candidate genes based on log-fold-change and other statistics.
#' The user can also pass a list of canidate null and non-null genes (they must be disjoint). In particular, the list should contain
#' two character vectors (for the name of the features/genes in the source data) with names 'null.genes' and 'nonnull.genes'.
#' For example, cand.DE.genes=list(null.genes=c('A', 'B'), nonnull.genes=c('C', 'D')).
#' @param pDE a numeric value between 0 and 1 indicating the desired fraction of DE genes in the simulated data
#' @param batch.config a numerical vector containing fractions for the composition of samples/cells per batch. The fractions must sum to 1. The number of batches to be simulated is equal to the size of the vector. (Example, batch.config=c(0.6, 0.4)   means simulate 2 batches with 60\% of the simulated samples/cells in batch 1 and the rest 40\% in the second batch. Another example,  batch.config=c(0.3, 0.35, 0.25) means simulate 3 batches with the first, second and third batches contain 30\%, 35\% and 25\% of the samples/cells, respectively).
#' @param group.config a numerical vector containing fractions for the composition of samples/cells per group. Similar definition to 'batch.config'. The number of groups to be simulated is equal to the size of the vector. The fractions must sum to 1.
#' @param lfc.thrld a positive numeric value for the minimum absolute log-fold-change for selecting
#' candidate DE genes in the source data (when group is not NULL, pDE>0 and cand.DE.genes is NULL)
#' @param t.thrld a positive numeric value for the minimum absolute t-test statistic for the
#' log-fold-changes of genes for selecting candidate DE genes in the source data (when group
#' is not NULL, pDE>0 and cand.DE.genes is NULL)
#' @param llStat.thrld a positive numeric value for the minimum squared test statistics from
#' the log-linear model to select candidate DE genes in the source data (when group is not NULL,
#' pDE>0 and cand.DE.genes is NULL)
#' @param model.zero.prob a logical value whether to model the zero expression probability separately (suitable for simulating of single-cell RNA-seq data or zero-inflated data)
#' @param tot.samples a numerical value for the total number of samples/cells to be simulated.
#' @param result.format a character value for the type of format for the output. Choice can be 'SCE' for SingleCellExperiment class or "list" for a list object that contains the simulated count, column information and row information.
#' @param return.details a logical value. If TRUE, detailed results including estimated parameters and densities will be returned
#' @param genewiseCor a logical value, if TRUE (default) the simulation will retain the gene-to-gene correlation structure of the source data using Gausian-copulas . Note that if it is TRUE, the program will be slow or it may fail for a limited memory size.
#' @param prior.count a positive constant to be added to the CPM before log transformation, to avoid log(0). The default is 1.
#' @param log.CPM.transform a logical value. If TRUE, the source data will be transformed into log-(CPM+const) before estimating the probability distributions
#' @param lib.size.params NULL or a named numerical vector containing parameters for simulating library sizes from log-normal distribution. If lib.size.params =NULL (default), 
#' then the package will fit a log-normal distribution for the library sizes in the source data to simulate new library sizes. 
#' If the user would like to specify the parameters of the log-normal distribution for the desired library sizes, 
#' then the log-mean and log-SD params of rlnorm() functions can be passed using this argument. 
#' Example, lib.size.params = c(meanlog=10, sdlog=0.2). See also ?rlnorm.
#' @param variable.lib.size a logical value. If FALSE (default), then the expected library sizes are simulated once and remains the same for every replication (if n.sim>1).
#' @param verbose a logical value, if TRUE a message about the status of the simulation will be printed on the console
#' @param w see ?hist
#' 
#' @return a list of SingleCellExperiment/list objects each containing simulated counts (not normalized), smple/cell level information in colData, and gene/feature level information in rowData.
#'
#' @details This function uses a specially designed exponential family for density estimation
#' to constructs the distribution of gene expression levels from a given real gene expression data
#' (e.g. single-cell or bulk sequencing data), and subsequently, simulates a new
#' from the estimated distributions.#'
#' For simulation of single-cell RNA-seq data (or any zero inflated gene expression data), the programm
#' involves an additional step to explicitly account for the high abundance of zero counts (if required).
#' This step models the probability of zero counts as a function the mean
#' expression of the gene and the library size of the cell (both in log scale) to add excess zeros.
#' This can be done by using \emph{model.zero.prob=TRUE}. Note that,
#' for extremly large size data, it is recomended to use a random sample of cells to
#' reduce computation time. To enable this, add the argument \emph{subset.data=TRUE} and you
#' can specify the number of cells to be used using \emph{n.samples} argument.
#' For example \emph{n.samples=400}.
#' Given known groups of samples/cells in the source data, DGE is simulated by independently
#' sampling data from distributions constructed for each group seprately. In particular, this procedure is
#' applied on a set of genes with absolute log-fold-change in the source data more than a given threshold (\emph{lfc.thrld}).
#' Moreover, when  the source dataset involves samples/cells processed in different batches, our
#' simulation procedure incorporates this batch effect in the simulated data, if required.
#' Different experimental designs can be simulated using the group and batch configuration arguments to
#' simulate biologica/experimental conditions and batchs, respectively.
#' Also, it is important to filter the source data so that genes with suffient expression will be used to
#' estimate the probability distributions.
#'
#' @references
#' \itemize{
#' \item Assefa, A. T., Vandesompele, J., & Thas, O. (2020). SPsimSeq: semi-parametric simulation of bulk and single cell RNA sequencing data. \emph{Bioinformatics}, doi: https://doi.org/10.1093/bioinformatics/btaa105.
#' \item Efron, B., & Tibshirani, R. (1996). Using specially designed exponential families for density estimation. \emph{The Annals of Statistics}, 24(6), 2431-2461.
#' }
#'
#' @examples
#' #----------------------------------------------------------------
#' # Example 1: simulating bulk RNA-seq
#'
#' # load the Zhang bulk RNA-seq data (availabl with the package)
#' data("zhang.data.sub")
#'
#' # filter genes with sufficient expression (important step to avoid bugs)
#' zhang.counts <- zhang.data.sub$counts[rowSums(zhang.data.sub$counts > 0)>=5, ]
#' MYCN.status  <- zhang.data.sub$MYCN.status
#'
#' # We simulate only a single data (n.sim = 1) with the following property
#' # - 2000 genes ( n.genes = 2000)
#' # - 20 samples (tot.samples = 20)
#' # - the samples are equally divided into 2 groups each with 90 samples
#' #   (group.config = c(0.5, 0.5))
#' # - all samples are from a single batch (batch = NULL, batch.config = 1)
#' # - we add 10% DE genes (pDE = 0.1)
#' # - the DE genes have a log-fold-change of at least 0.5 in
#' #   the source data (lfc.thrld = 0.5)
#' # - we do not model the zeroes separately, they are the part of density
#' #    estimation (model.zero.prob = FALSE)
#'
#' # simulate data
#' set.seed(6452)
#' zhang.counts2 <- zhang.counts[sample(nrow(zhang.counts), 2000), ]
#' sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts2, batch = NULL,
#'                           group = MYCN.status, n.genes = 2000, batch.config = 1,
#'                           group.config = c(0.5, 0.5), tot.samples = 20,
#'                           pDE = 0.1, lfc.thrld = 0.5, result.format = "list")
#'
#'
#' head(sim.data.bulk$counts[[1]][, seq_len(5)])  # count data
#' head(sim.data.bulk$colData)        # sample info
#' head(sim.data.bulk$rowData)        # gene info
#'
#'
#' #----------------------------------------------------------------
#' # Example 2: simulating single cell RNA-seq from a single batch (read-counts)
#' # we simulate only a single scRNA-seq data (n.sim = 1) with the following property
#' # - 2000 genes (n.genes = 2000)
#' # - 100 cells (tot.samples = 100)
#' # - the cells are equally divided into 2 groups each with 50 cells
#' #   (group.config = c(0.5, 0.5))
#' # - all cells are from a single batch (batch = NULL, batch.config = 1)
#' # - we add 10% DE genes (pDE = 0.1)
#' # - the DE genes have a log-fold-change of at least 0.5
#' #' # - we model the zeroes separately (model.zero.prob = TRUE)
#' # - the ouput will be in SingleCellExperiment class object (result.format = "SCE")
#'
#' library(SingleCellExperiment)
#'
#' # load the NGP nutlin data (availabl with the package, processed with
#' # SMARTer/C1 protocol, and contains read-counts)
#' data("scNGP.data")
#'
#' # filter genes with sufficient expression (important step to avoid bugs)
#' scNGP.data2 <- scNGP.data[rowSums(counts(scNGP.data) > 0)>=5, ]
#' treatment <- ifelse(scNGP.data2$characteristics..treatment=="nutlin",2,1)
#'
#' set.seed(654321)
#' scNGP.data2 <- scNGP.data2[sample(nrow(scNGP.data2), 2000), ]
#'
#' # simulate data (we simulate here only a single data, n.sim = 1)
#' #sim.data.sc <- SPsimSeq(n.sim = 1, s.data = scNGP.data2, batch = NULL,
#'  #                       group = treatment, n.genes = 2000, batch.config = 1,
#'   #                      group.config = c(0.5, 0.5), tot.samples = 100,
#'    #                     pDE = 0.1, lfc.thrld = 0.5, model.zero.prob = TRUE,
#'     #                    result.format = "SCE")
#'
#' #sim.data.sc1 <- sim.data.sc[[1]]
#' #class(sim.data.sc1)
#' #head(counts(sim.data.sc1)[, seq_len(5)])
#' #colData(sim.data.sc1)
#' #rowData(sim.data.sc1)
#'
#' @export
#' @importFrom stats pnorm dnorm runif rbinom predict approx quantile glm rlnorm lm sd var setNames
#' @importFrom plyr rbind.fill
SPsimSeq <- function(n.sim = 1, s.data, batch = rep(1, ncol(s.data)), 
                     group = rep(1, ncol(s.data)), 
                     n.genes = 1000, batch.config = 1, group.config = 1, 
                     pDE = 0.1, cand.DE.genes = NULL, lfc.thrld = 0.5, 
                     t.thrld = 2.5, llStat.thrld = 5, tot.samples = ncol(s.data), 
                     model.zero.prob = FALSE, genewiseCor = TRUE,
                     log.CPM.transform = TRUE, lib.size.params = NULL, 
                     variable.lib.size = FALSE,
                     result.format = "SCE", return.details = FALSE, 
                     verbose = TRUE, prior.count = 1, const.mult = 1e6)
{
  #DATA EXTRACTION
  s.data = extractMat(s.data)
  #INPUT CHECKS
  checkInputs <- checkInputValidity(s.data = s.data, group = group, batch = batch,
                                    group.config = group.config, batch.config = batch.config, 
                                    w = w, log.CPM.transform = log.CPM.transform, 
                                    prior.count = prior.count, lib.size.params = lib.size.params)
  #CPM TRANSFORM
  cpm.data <- if(log.CPM.transform){
    # calculate log CPM 
     calculateCPM(s.data, prior.count = prior.count, 
                             const.mult = const.mult)
  } else s.data
  #PARAMETER ESTIMATION
  # Estimate library size distributions
  LS = colSums(s.data) #Library sizes
  if(variable.lib.size && is.null(lib.size.params)){
    if(verbose){message("Fitting library size distirbution ...")}
    lib.size.params <- estLibSizeDistr(LS = LS, batch = batch)
  }
  # Estimate correlation matrices
  if(genewiseCor){
    if(verbose) {message("Estimating featurewise correlations ...")}
    corMats.batch <- obtCorMatsBatch(cpm.data = cpm.data, batch = batch)
  }
  #Find candidate DE genes
  if(is.null(cand.DE.genes)){
    if(verbose) {message("Selecting candidate DE genes ...")}
    cand.DE.genes = if(length(unique(group))>1){ 
      chooseCandGenes(cpm.data = cpm.data, group = group, const = prior.count, 
                      lfc.thrld = lfc.thrld, t.thrld = t.thrld, 
                      llStat.thrld = llStat.thrld, pDE = pDE, n.genes = n.genes)
    } else {
      list(null.genes = rownames(s.data)) 
    } 
  }
  #Selected genes
  null.genes0    <- cand.DE.genes$null.genes
  nonnull.genes0 <- cand.DE.genes$nonnull.genes
  # Fit logistic regression for the probability of zeros
  if(all(s.data!=0)) model.zero.prob = FALSE
  if(model.zero.prob){
    if(verbose) {message("Fitting zero probability model ...")}
    fracZero.logit.list <- fracZeroLogitModel(s.data = s.data, batch = batch, 
                                              cpm.data = cpm.data)
  }
  # Estimate batch specific densities
  if(verbose) {message("Estimating densities ...")}
  densList <- lapply(c(null.genes0, nonnull.genes0), function(gene){ 
    geneParmEst(cpm.data.i = cpm.data[gene, ], batch = batch, group = group,
                de.ind = gene %in% nonnull.genes0,
                model.zero.prob = model.zero.prob, min.val = log(prior.count), w = w)
  })
  densList = setNames(densList, c(null.genes0, nonnull.genes0))
  #SIMULATION
  ## EXPERIMENT CONFIGURATION
  if(verbose) {message("Configuring design ...")}
  exprmt.design <- configExperiment(batch.config = batch.config, group.config = group.config,
                                    tot.samples = tot.samples)
  ## DATA GENERATION
  if(verbose) {message("Simulating data ...")}
  sim.data.list <- lapply(seq_len(n.sim), function(h){
    if(verbose) {message(" ...", h, " of ", n.sim)}
    #Sample libray sizes
    samLS <- if(variable.lib.size & log.CPM.transform){
       genLibSizes(fit.ln = lib.size.params, n.batch = n.batch, 
                       n.group = n.group, config.mat = config.mat)
    } else LS
    #Sample copula
    copSam = genCopula(corMats.batch, n.batch = n.batch, batch = batch)
    
    # sample DE and null genes
    selctGenes <- selectGenes(pDE = pDE, group = group, n.genes = n.genes,
                                 null.genes0 = null.genes0, nonnull.genes0 = nonnull.genes0,
                                 group.config = group.config)
    #Generate data
    sim.dat <- lapply(selctGenes, function(gene){ 
      SPsimPerGene(cpm.data = cpm.data, densList.ii = densList[[gene]],
                   DE.ind.ii = gene %in% nonnull.genes0, sel.genes.ii = gene,
                   n.batch = n.batch, n.group = n.group, group = group, batch=batch,
                   log.CPM.transform = log.CPM.transform, config.mat=config.mat,
                   null.group=null.group, LL = if(variable.lib.size) samLS else LL, copulas.batch=copulas.batch,
                   const = const, min.val = min.val, model.zero.prob=model.zero.prob,
                   tot.samples=tot.samples, fracZero.logit.list = fracZero.logit.list)
    })
    sim.data.h <- prepareSPsimOutputs(sim.dat=sim.dat, n.batch=n.batch, n.group=n.group,
                        DE.ind = selctGenes %in% nonnull.genes0, 
                        sel.genes = selctGenes, LL=ELS, 
                        result.format=result.format, log.CPM.transform=log.CPM.transform)
    return(sim.data.h)
  })
  if(return.details){
    list("sim.data.list" = sim.data.list, "detailed.results" = densList)
  }else{
    sim.data.list
  }
  
}