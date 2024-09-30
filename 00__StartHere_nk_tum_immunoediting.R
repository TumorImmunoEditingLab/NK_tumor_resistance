# Specifying parameters here ----
# remotes::install_github("rstudio/renv@v0.16.0")
if(!("renv" %in% installed.packages())){install.packages("renv", version="0.15.4")}
if(!("import" %in% installed.packages())){renv::install("import")}
if(!("yaml" %in% installed.packages())){renv::install("yaml")}
if(!("knitr" %in% installed.packages())){renv::install("knitr")}

# setting up initial parameters
base_dir <- "/home/rstudio/workspace/"
data_dir <- file.path(base_dir, "datasets") 
dir.create(data_dir)
results_dir <- file.path(base_dir, "results_dir") 
dir.create(results_dir)
setwd(base_dir)

# project setup
# loading config file ----
config_file <- file.path(base_dir, "nk_tum_immunoediting_config.yaml")  
config <- yaml::read_yaml(config_file)

#project setup
project_name <- config$project_name
nthreads <- config$nthreads # e.g. 4 for furrr multisession

# this make take a while...
# https://rstudio.github.io/renv/reference/config.html
# getOption(x, default = NULL)
# renv::settings$use.cache()
# getOption('renv.config.pak.enabled')
# if the project does not automatically activate run:
if(Sys.getenv("RENV_PATHS_CACHE") != "/renv_cache") {Sys.setenv(RENV_PATHS_CACHE = "/renv_cache")}
if(Sys.getenv("RENV_PATHS_LIBRARY") != "/home/rstudio/renv_library") {Sys.setenv(RENV_PATHS_LIBRARY = "/home/rstudio/renv_library")}
#if(Sys.getenv("RENV_CONFIG_PAK_ENABLED") != "TRUE") {Sys.setenv(RENV_CONFIG_PAK_ENABLED = TRUE)} # add pak, targets and benchmarkme to docker!!!
# setting root dir
# 0. renv::activate
renv::activate(project = base_dir)

#IF packages are not loaded properly, use: 
# renv::restore()

# Only for manual installation of the project and recovery from the lock file.
# renv::init(project = base_dir, bare=TRUE, bioconductor = "3.16")
# FOR the first run restore environment and packages from renv.lock 
# 1. restore original environment
# renv::restore(project = base_dir, lockfile = file.path(base_dir, paste0(project_name, "_renv.lock")), prompt = FALSE)
# file.copy(from = file.path(base_dir, paste0(project_name, "_renv.lock")), to = file.path(base_dir, "renv.lock")) # to use with renv::diagnostics()

# for some reason apeglm package cannot be spanshotted to the renv.lock file. so, it must be installed manually.
BiocManager::install("apeglm")



# Load the GTF from ENSEMBL
download.file(url = "https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz",
              destfile = "~/workspace/datasets/annotations/Mus_musculus.GRCm38.102.gtf.gz")
R.utils::gunzip(filename = "~/workspace/datasets/annotations/Mus_musculus.GRCm38.102.gtf.gz")

######### 
# PUT HERE PULLS FROM GEO 


