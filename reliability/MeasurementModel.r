
## ########################################### ##
##         SOFTWARE LICENSE AGREEMENT          ##
## ########################################### ##

# This work is licensed under a Creative Commons Attribution 4.0 International License.
# See http://creativecommons.org/licenses/by/4.0/ for details on license restrictions.
# Please attribute any work that makes use or is derived from this work
#   by refering to https://github.com/jalmar/openmx-models


## ########################################### ##
##             DESCRIPTION OF FILE             ##
## ########################################### ##

# In this R script, a measurement model with two latent factors (nf=2) is provided.
# This measurement model can be used to obtain a estimate of the association between
# the reliable components of two constructs that is free of random measurement error.
# In this example, there are two constructs with a total of three observed variables (nv=3).
# Two parallel scores "M1_H1" and "M1_H2" are available for the first construct "M1".
# The first latent factor "F1" is derived from these two parallel scores.
# Since the two parallel scores are measurements of the same construct,
# their expected means are constrained to be equal.
# The second latent factor "F1" is derived from the single score of the second construct "M2".
# If there are repeated measurements of the second construct "M2", the measurement model
# can easily be expanded for reliablity modelling of the second construct "M2".
# The quality of the fit is assessed with the Comparative Fit Index (CFI),
# the Tucker-Lewis Index (TLI), and the Root Mean Square Error of the Approximation (RMSEA).
# These fit indices are based on the comparison between the measurement model to
# a saturated model and an independence (or null) model. The saturated model is a model that
# estimates all means, variances and covariances freely for all variables (i.e. unconstrained).
# The independence model estimates only the means and variances with covariances fixed to zero.
# Finally, a summary of the most relevant parameters is printed.
# For more details, see https://github.com/jalmar/openmx-models/reliability


## ############################################################################


## ########################################### ##
##       LOAD LIBRARIES AND SET OPTIONS        ##
## ########################################### ##

## load OpenMx library
require(OpenMx) # SEM
options(warn=1)


## ########################################### ##
##    INPUT DATASETS AND SELECTED VARIABLES    ##
## ########################################### ##

## INPUT: a data frame data_df with the parallel and single scores measurements
data_df <- readRDS("ExampleData.rds")
selVars <- c("FC_Z_H1", "FC_Z_H2", "BHV")


## ########################################### ##
##      CALCULATE MEANS AND (CO)VARIANCES      ##
## ########################################### ##

## calculate initial mean and variance of measures
init_means <- apply(data_df[,selVars, drop=F], 2, mean, na.rm=T)
init_vars <- apply(data_df[,selVars, drop=F], 2, var, na.rm=T)
init_sds <- apply(data_df[,selVars, drop=F], 2, sd, na.rm=T)

## calculate initial covariance matrices between measures
init_covars <- var(data_df[,selVars, drop=F], na.rm=T)
init_corrs <- cor(data_df[,selVars, drop=F], method="pearson", use="pairwise.complete.obs")


## ########################################### ##
##        MEASUREMENT MODEL PARAMETERS         ##
## ########################################### ##

## number of observed variables
nf=2 # number of factors
nv=3 # number of measurements

## initial guess at factor variances [nf x nf matrices]; NOTE: standardized symmetric matrix with variances on diagnoal fixed to 1.0!
# initialize with rough estimate of expected correlation; e.g. corr(M1_H1, M2)
init_factor_variance_values = c(1.0, 0.5,
                                0.5, 1.0)
init_factor_variance_paths = c(F, T,
                               T, F)
init_factor_variance_labels = c("varF1", "covF12",
                                "covF12", "varF2")

## loading of factors on the measures [nf x nv matrices]:
# factor 1 loads on the two parallel scores of measure 1
# factor 2 loads on the single score of measure 2
# initialize with rough estimate of factor loadings; e.g. cov(M1_H1, M1_H2))
init_factor_loading_values = c(0.70, 0.70, 0.00, # factor 1
                               0.00, 0.00, 1.00) # factor 2
init_factor_loading_paths = c(T, T, F, # factor 1
                              F, F, F) # factor 2; NOTE: third position fixed for single score measurement
init_factor_loading_labels = c("f1", "f1", NA, # factor 1
                                NA,   NA, "f2") # factor 2

## loadings of residual variances on the measures [1 x nv matrices]
# initialize with rough estimate of residual variance; e.g. 1 - sqrt(corr(M1_H1, M1_H2)) or 1 - factor loading
init_es_loading_values = c(0.30, 0.30, 0.00)
init_es_loading_paths = c(T, T, F) # NOTE: third position fixed for single score measurement
init_es_loading_labels = paste0("es", 1:nv)

## restrict the means of variables to be equal through shared labels [1 x nv matrices]
# Here, the means for parallel scores of the same construct are contrained to be equal
init_means_labels = c("mean_FC_Z", "mean_FC_Z", "mean_BHV")


## ########################################### ##
##  OPENMX SPECIFICATION OF MEASUREMENT MODEL  ##
## ########################################### ##

## measurement model
MeasurementModel <- mxModel("mm",
    
    # variance matrix Vf to store total variance for latent phenotypic factor; NOTE: standardized with variance on diagonals fixed to 1.0
    mxMatrix(type="Full", nrow=nf, ncol=nf, free=init_factor_variance_paths, values=init_factor_variance_values, labels=init_factor_variance_labels, name="Vf", lbound=-1, ubound=1),
    
    # matrices for correlations between latent factors (only nf>1)
    mxMatrix(type="Iden", ncol=nf, nrow=nf, name="Idnf"), # identity matrix of size nf x nf
    mxAlgebra(expression=solve(sqrt(Idnf*Vf)) %*% Vf %*% solve(sqrt(Idnf*Vf)), name="Rphf"),
    
    # matrix fs to store path coefficients for latent phenotypic factors loading on observed variables
    mxMatrix(type="Full", nrow=nv, ncol=nf, free=init_factor_loading_paths, values=init_factor_loading_values*init_vars, labels=init_factor_loading_labels, name="fs", lbound=0),
    
    # matrix es to store path coefficients for measurement-specific factors loading on each measurement; NOTE: fixing off-diagonals to zero assumes no correlation between residual variances of measurements
    mxMatrix(type="Diag", nrow=nv, ncol=nv, free=init_es_loading_paths, values=init_es_loading_values*init_vars, labels=init_es_loading_labels, name="es", lbound=0),
    
    # matrix Es to store variance for measurement-specific factors
    mxAlgebra(expression=es%*%t(es), name="Es"),
    
    # total variance of the measures is the proportion of variance explained by the factor(s) plus residual variances
    mxAlgebra(expression=fs%&%Vf+Es, name="Vm"), # (co)variances between the measures
    mxMatrix(type="Iden", nrow=nv, ncol=nv, free=F, name="Inv"),
    mxAlgebra(expression=solve(sqrt(Inv*Vm)), name="SDm"), # standard deviation of the measures
    
    # phenotypic correlation between measures
    mxAlgebra(expression=solve(sqrt(Inv*Vm)) %*% Vm %*% solve(sqrt(Inv*Vm)), name="Rphm"),
    
    # standardized estimates for loadings on latent factor
    mxAlgebra(expression=SDm%*%fs, "Fs_std"),
    mxAlgebra(expression=Fs_std*Fs_std, "Fs_std2"),
    
    # standardized estimates for measurement-specific factors
    mxAlgebra(expression=SDm%*%es, name="Es_std"),
    mxAlgebra(expression=Es_std*Es_std, name="Es_std2"),
    
    # observed data of subjects
    mxData(observed=data_df, type="raw"),
    
    # expected means vector
    mxMatrix(type="Full", nrow=1, ncol=nv, free=T, values=init_means, labels=init_means_labels, name="expMeans"),
    
    # expected covariance matrix; NOTE: %&% is quadratic product (A %&% B == ABA'); same as Vm
    mxAlgebra(expression=mm.fs %&% mm.Vf + mm.Es, name="expCov"),
        
    # optimization objective
    mxExpectationNormal(means="expMeans", covariance="expCov", dimnames=selVars),
    mxFitFunctionML(),
    
    # calculate confidence intervals of freely estimated parameters
    mxCI(sprintf("mm.Fs_std2[%s]", apply(which(matrix(init_factor_loading_paths, nrow=nv, ncol=nf), arr.ind = TRUE), MARGIN=1, paste0, collapse=",")), interval = 0.95, type="both"), # standardized factor loadings
    mxCI(sprintf("mm.Rphf[%s]", apply(which(lower.tri(matrix(1,nf,nf)), arr.ind = TRUE), MARGIN=1, paste0, collapse=",")), interval = 0.95, type="both"), # correlation between factors
    mxCI(sprintf("mm.Rphm[%s]", apply(which(lower.tri(matrix(1,nv,nv)), arr.ind = TRUE), MARGIN=1, paste0, collapse=",")), interval = 0.95, type="both") # correlation between measures
    
) # END of mxModel MeasurementModel

## assign first value to parameters with same labels but different starting values
MeasurementModel <- omxAssignFirstParameters(MeasurementModel)


## ########################################### ##
##          FIT THE MEASUREMENT MODEL          ##
## ########################################### ##

MeasurementModelFit <- mxTryHard(MeasurementModel, extraTries=30, checkHess=F, intervals=T, silent=T)


## ########################################### ##
##  EXTRACT PARAMETERS OF INTEREST FROM MODEL  ##
## ########################################### ##

## summary of measurement model model fit
print(summary(MeasurementModelFit))

## corrected/reliable association between the two constructs
message(sprintf("Corrected association between the two constructs = %0.3f [%0.3f; %0.3f]", MeasurementModelFit$output$algebras$mm.Rphf[2,1], MeasurementModelFit$output$confidenceIntervals["mm.Rphf[2,1]","lbound"], MeasurementModelFit$output$confidenceIntervals["mm.Rphf[2,1]","ubound"]))

## standardized factor loading on parallel scores of the first construct
message(sprintf("Standardized factor loading on first parallel score = %2.1f%% [%2.1f%%; %2.1f%%]", MeasurementModelFit$output$algebras$mm.Fs_std2[1,1]*100, MeasurementModelFit$output$confidenceIntervals["mm.Fs_std2[1,1]","lbound"]*100, MeasurementModelFit$output$confidenceIntervals["mm.Fs_std2[1,1]","ubound"]*100))
message(sprintf("Standardized factor loading on second parallel score = %2.1f%% [%2.1f%%; %2.1f%%]", MeasurementModelFit$output$algebras$mm.Fs_std2[2,1]*100, MeasurementModelFit$output$confidenceIntervals["mm.Fs_std2[2,1]","lbound"]*100, MeasurementModelFit$output$confidenceIntervals["mm.Fs_std2[2,1]","ubound"]*100))

## test-retest reliability between the two parallel scores
message(sprintf("Test-retest reliability between the two parallel scores = %0.3f [%0.3f; %0.3f]", MeasurementModelFit$output$algebras$mm.Rphm[2,1], MeasurementModelFit$output$confidenceIntervals["mm.Rphm[2,1]","lbound"], MeasurementModelFit$output$confidenceIntervals["mm.Rphm[2,1]","ubound"]))

## uncorrected associations of the individual parallel scores with the second measurement
message(sprintf("Uncorrected association between the first parallel score and second construct = %0.3f [%0.3f; %0.3f]", MeasurementModelFit$output$algebras$mm.Rphm[3,1], MeasurementModelFit$output$confidenceIntervals["mm.Rphm[3,1]","lbound"], MeasurementModelFit$output$confidenceIntervals["mm.Rphm[3,1]","ubound"]))
message(sprintf("Uncorrected association between the second parallel score and second construct = %0.3f [%0.3f; %0.3f]", MeasurementModelFit$output$algebras$mm.Rphm[3,2], MeasurementModelFit$output$confidenceIntervals["mm.Rphm[3,2]","lbound"], MeasurementModelFit$output$confidenceIntervals["mm.Rphm[3,2]","ubound"]))


## ########################################### ##
##          DETERMINE GOODNESS OF FIT          ##
## ########################################### ##

## specification of the saturated model
SaturatedModel <- mxModel("sat",
    
    # expected (co)variance matrix; unconstrained variances and covariance estimates
    mxMatrix(type="Symm", nrow=nv, ncol=nv, free=T, values=vech(init_covars), name="expCov"),
    
    # observed data of subjects
    mxData(observed=data_df, type="raw"),
    
    # expected means vector; unconstrained means estimates
    mxMatrix(type="Full", nrow=1, ncol=nv, free=T, values=init_means, labels=paste0("mean_m",1:nv), name="expMeans"),
        
    # optimization objective
    mxExpectationNormal(means="expMeans", covariance="expCov", dimnames=selVars),
    mxFitFunctionML()
    
) # END of mxModel SaturatedModel

SaturatedModel <- omxAssignFirstParameters(SaturatedModel)

## specification of the independence/null model
IndependenceModel <- mxModel("indep",
    
    # expected (co)variance matrix; unconstrained variances and covariance estimates
    mxMatrix(type="Diag", nrow=nv, ncol=nv, free=T, values=init_vars, name="expCov"),
    
    # observed data of subjects
    mxData(observed=data_df, type="raw"),
    
    # expected means vector; same means constraints as substantive model
    mxMatrix(type="Full", nrow=1, ncol=nv, free=T, values=init_means, labels=init_means_labels, name="expMeans"),
        
    # optimization objective
    mxExpectationNormal(means="expMeans", covariance="expCov", dimnames=selVars),
    mxFitFunctionML()
    
) # END of mxModel IndependenceModel

IndependenceModel <- omxAssignFirstParameters(IndependenceModel)

## fit saturated and independence model
SaturatedModelFit <- mxTryHard(SaturatedModel, extraTries=10, checkHess=F, intervals=F)
IndependenceModelFit <- mxTryHard(IndependenceModel, extraTries=10, checkHess=F, intervals=F)

## print goodness of fit indices
tmp_summary <- summary(MeasurementModelFit, refModels=list(Saturated=SaturatedModelFit, Independence=IndependenceModelFit))
message(sprintf("Goodness of fit indices are CFI=%0.3f, TLI=%0.3f, and RMSEA=%0.3f", tmp_summary$CFI, tmp_summary$TLI, tmp_summary$RMSEA))

