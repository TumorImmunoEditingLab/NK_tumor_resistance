# importing only key functions that are actually used - not to polute namespace!
import::from(.from = readr, read_csv, cols)
import::from(magrittr, "%>%")
import::from(dplyr, mutate, select, filter, rename, arrange, desc, group_by, summarise, ungroup)  # dplyr_mutate = mutate
import::from(purrr, map)
import::from(future, plan, multisession, sequential)
import::from(furrr, furrr_options, future_map2)
import::from(ggplot2, .all=TRUE) # importing all as there is too many
import::from(grid, gpar) # needed in complexheatmap
import::from(kableExtra, kable_styling, kbl)

import::from(.from = GenomicFeatures, makeTxDbFromGFF)
import::from(.from = AnnotationDbi, annot_db_keys = keys, annot_db_select = select)
import::from(.from = DESeq2, .all=TRUE)
import::from(.from = tximport, tximport)

import::from(.from = here::here("utils/filterDatasets.R"), "filterDatasets", .character_only=TRUE) # used for filtering
import::from(.from = here::here("utils/generateEnsemblAnnotation.R"), "generateEnsemblAnnotation", .character_only=TRUE) # used for filtering