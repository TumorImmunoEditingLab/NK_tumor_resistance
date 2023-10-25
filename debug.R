GDCprepare <- function(
    query,
    save = FALSE,
    save.filename,
    directory = "GDCdata",
    summarizedExperiment = TRUE,
    remove.files.prepared = FALSE,
    add.gistic2.mut = NULL,
    mut.pipeline = "mutect2",
    mutant_variant_classification = c(
      "Frame_Shift_Del",
      "Frame_Shift_Ins",
      "Missense_Mutation",
      "Nonsense_Mutation",
      "Splice_Site",
      "In_Frame_Del",
      "In_Frame_Ins",
      "Translation_Start_Site",
      "Nonstop_Mutation"
    )
){
  
  isServeOK()
  if(missing(query)) stop("Please set query parameter")
  
  test.duplicated.cases <-
    (
      any(
        duplicated(query$results[[1]]$cases)) & # any duplicated
        !(query$data.type %in% c(
          "Clinical data",
          "Protein expression quantification",
          "Raw intensities",
          "Masked Intensities",
          "Clinical Supplement",
          "Masked Somatic Mutation",
          "Biospecimen Supplement"
        )
        )
    )
  
  
  if(test.duplicated.cases) {
    dup <- query$results[[1]]$cases[duplicated(query$results[[1]]$cases)]
    cols <- c("tags","cases","experimental_strategy","analysis_workflow_type")
    cols <- cols[cols %in% colnames(query$results[[1]])]
    dup <- query$results[[1]][query$results[[1]]$cases %in% dup,cols]
    dup <- dup[order(dup$cases),]
    print(knitr::kable(dup))
    stop("There are samples duplicated. We will not be able to prepare it")
  }
  
  if (!save & remove.files.prepared) {
    stop("To remove the files, please set save to TRUE. Otherwise, the data will be lost")
  }
  # We save the files in project/source/data.category/data.type/file_id/file_name
  source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(
    query$results[[1]]$project, source,
    gsub(" ","_",query$results[[1]]$data_category),
    gsub(" ","_",query$results[[1]]$data_type),
    gsub(" ","_",query$results[[1]]$file_id),
    gsub(" ","_",query$results[[1]]$file_name)
  )
  
  files <- file.path(directory, files)
  
  # For IDAT prepare since we need to put all IDATs in the same folder the code below will not work
  # a second run
  if (!all(file.exists(files))) {
    # We have to check we movedthe files
    if (query$data.type == "Masked Intensities" | query$data.category == "Raw microarray data"){
      files.idat <- file.path(
        query$results[[1]]$project, source,
        gsub(" ","_",query$results[[1]]$data_category),
        gsub(" ","_",query$results[[1]]$data_type),
        gsub(" ","_",query$results[[1]]$file_name)
      )
      files.idat <- file.path(directory, files.idat)
      if (!all(file.exists(files) | file.exists(files.idat))) {
        stop(
          paste0(
            "I couldn't find all the files from the query. ",
            "Please check if the directory parameter is right ",
            "or `GDCdownload` downloaded the samples."
          )
        )
      }
    } else {
      stop(
        paste0(
          "I couldn't find all the files from the query. ",
          "Please check if the directory parameter is right ",
          "or `GDCdownload` downloaded the samples."
        )
      )
    }
  }
  
  cases <- ifelse(grepl("TCGA|TARGET",query$results[[1]]$project %>% unlist()),query$results[[1]]$cases,query$results[[1]]$sample.submitter_id)
  
  if (grepl("Transcriptome Profiling", query$data.category, ignore.case = TRUE)){
    
    if(unique(query$results[[1]]$experimental_strategy) == "scRNA-Seq"){
      #if (grepl("Single Cell Analysis", unique(query$results[[1]]$data_type), ignore.case = TRUE)){
      
      data <- readSingleCellAnalysis(
        files = files,
        data_format = unique(query$results[[1]]$data_format),
        workflow.type = unique(query$results[[1]]$analysis_workflow_type),
        cases = cases
      )
      return(data)
    } else {
      data <- readTranscriptomeProfiling(
        files = files,
        data.type = ifelse(!is.na(query$data.type),  as.character(query$data.type),  unique(query$results[[1]]$data_type)),
        workflow.type = unique(query$results[[1]]$analysis_workflow_type),
        cases = cases,
        summarizedExperiment
      )
    }
  } else if(grepl("Copy Number Variation",query$data.category,ignore.case = TRUE)) {
    if (unique(query$results[[1]]$data_type) == "Gene Level Copy Number Scores") {
      data <- readGISTIC(files, query$results[[1]]$cases)
    } else if (unique(query$results[[1]]$data_type) == "Gene Level Copy Number") {
      data <- readGeneLevelCopyNumber(files, query$results[[1]]$cases,summarizedExperiment = summarizedExperiment)
    } else {
      data <- readCopyNumberVariation(files, query$results[[1]]$cases)
    }
  }  else if (grepl("Methylation Beta Value",query$data.type, ignore.case = TRUE)) {
    data <- readDNAmethylation(files, cases = cases, summarizedExperiment, unique(query$platform))
  }  else if (grepl("Raw intensities|Masked Intensities",query$data.type, ignore.case = TRUE)) {
    # preparing IDAT files
    data <- readIDATDNAmethylation(files, barcode = cases, summarizedExperiment, unique(query$platform), query$legacy)
  }  else if (grepl("Proteome Profiling",query$data.category,ignore.case = TRUE)) {
    data <- readProteomeProfiling(files, cases = cases)
  }  else if (grepl("Protein expression",query$data.category,ignore.case = TRUE)) {
    data <- readProteinExpression(files, cases = cases)
    if(summarizedExperiment) message("SummarizedExperiment not implemented, if you need samples metadata use the function TCGAbiolinks:::colDataPrepare")
  }  else if (grepl("Simple Nucleotide Variation",query$data.category,ignore.case = TRUE)) {
    if(grepl("Masked Somatic Mutation",query$results[[1]]$data_type[1],ignore.case = TRUE) | source == "legacy")
      suppressWarnings(data <- readSimpleNucleotideVariationMaf(files))
  }  else if (grepl("Clinical|Biospecimen", query$data.category, ignore.case = TRUE)){
    data <- readClinical(files, query$data.type, cases = cases)
    summarizedExperiment <- FALSE
  } else if (grepl("Gene expression",query$data.category,ignore.case = TRUE)) {
    if (query$data.type == "Gene expression quantification")
      data <- readGeneExpressionQuantification(
        files = files,
        cases = cases,
        summarizedExperiment = summarizedExperiment,
        genome = ifelse(query$legacy,"hg19","hg38"),
        experimental.strategy = unique(query$results[[1]]$experimental_strategy)
      )
    
    if (query$data.type == "miRNA gene quantification")
      data <- readGeneExpressionQuantification(
        files = files,
        cases = cases,
        summarizedExperiment = FALSE,
        genome = ifelse(query$legacy,"hg19","hg38"),
        experimental.strategy = unique(query$results[[1]]$experimental_strategy)
      )
    
    if (query$data.type == "miRNA isoform quantification")
      data <- readmiRNAIsoformQuantification(
        files = files,
        cases = query$results[[1]]$cases
      )
    
    if (query$data.type == "Isoform expression quantification")
      data <- readIsoformExpressionQuantification(files = files, cases = cases)
    
    if (query$data.type == "Exon quantification")
      data <- readExonQuantification(
        files = files,
        cases = cases,
        summarizedExperiment = summarizedExperiment
      )
    
  }
  # Add data release to object
  if (summarizedExperiment & !is.data.frame(data)) {
    metadata(data) <- list("data_release" = getGDCInfo()$data_release)
  }
  
  
  if ((!is.null(add.gistic2.mut)) & summarizedExperiment) {
    message("=> Adding GISTIC2 and mutation information....")
    genes <- tolower(levels(EAGenes$Gene))
    if (!all(tolower(add.gistic2.mut) %in% genes)) {
      message(
        paste("These genes were not found:\n",
              paste(add.gistic2.mut[! tolower(add.gistic2.mut) %in% genes],collapse = "\n=> ")
        )
      )
    }
    add.gistic2.mut <- add.gistic2.mut[tolower(add.gistic2.mut) %in% tolower(genes)]
    if (length(add.gistic2.mut) > 0){
      info <- colData(data)
      for(i in unlist(query$project)){
        info <- get.mut.gistc.information(
          info,
          i,
          add.gistic2.mut,
          mut.pipeline = mut.pipeline,
          mutant_variant_classification = mutant_variant_classification
        )
      }
      colData(data) <- info
    }
  }
  if("samples" %in% colnames(data)){
    if(any(duplicated(data$sample))) {
      message("Replicates found.")
      if(any(data$is_ffpe)) message("FFPE should be removed. You can modify the data with the following command:\ndata <- data[,!data$is_ffpe]")
      print(as.data.frame(colData(data)[data$sample %in% data$sample[duplicated(data$sample)],c("is_ffpe"),drop = FALSE]))
    }
  }
  
  
  if(save){
    if(missing(save.filename) & !missing(query)) save.filename <- paste0(query$project,gsub(" ","_", query$data.category),gsub(" ","_",date()),".RData")
    message(paste0("=> Saving file: ",save.filename))
    save(data, file = save.filename)
    message("=> File saved")
    
    # save is true, due to the check in the beggining of the code
    if(remove.files.prepared){
      # removes files and empty directories
      remove.files.recursively(files)
    }
  }
  
  return(data)
}


readClinical <- function(files, data.type, cases){
  if(data.type == "Clinical data"){
    suppressMessages({
      ret <- plyr::alply(files,.margins = 1,readr::read_tsv, .progress = "text")
    })
    names(ret) <- gsub("nationwidechildrens.org_","",gsub(".txt","",basename(files)))
  } else if(data.type %in% c("Clinical Supplement","Biospecimen Supplement")){
    ret <- plyr::alply(files,.margins = 1,function(f) {
      readr::read_tsv(f,col_types = readr::cols())
    }, .progress = "text")
    names(ret) <- gsub("nationwidechildrens.org_","",gsub(".txt","",basename(files)))
  }
  return(ret)
}

readGeneExpressionQuantification <- function(
    files,
    cases,
    genome = "hg19",
    summarizedExperiment = TRUE,
    experimental.strategy,
    platform
){
  skip <- unique((ifelse(experimental.strategy == "Gene expression array",1,0)))
  
  if(length(skip) > 1) stop("It is not possible to handle those different platforms together")
  
  print.header(paste0("Reading ", length(files)," files"),"subsection")
  ret <- plyr::alply(
    .data = seq_along(files),
    .margins = 1,
    .fun = function(i,cases){
      
      data <- fread(
        input = files[i],
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        skip = skip
      )
      
      if(!missing(cases)) {
        assay.list <<- gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)])
        # We will use this because there might be more than one col for each samples
        setnames(
          data,
          colnames(data)[2:ncol(data)],
          paste0(gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)]),"_",cases[i])
        )
      }
      data
    },.progress = "time",cases = cases)
  
  print.header(paste0("Merging ", length(files)," files"),"subsection")
  
  # Just check if the data is in the same order, since we will not merge
  # the data frames to save memory
  stopifnot(all(unlist(ret %>% map(function(y){all(y[,1] ==  ret[[1]][,1])}) )))
  
  # need to check if it works in all cases
  df <- ret %>% map( ~ (.x %>% dplyr::select(-1))) %>% bind_cols()
  df <- bind_cols(ret[[1]][,1],df)
  
  if (summarizedExperiment) {
    df <- makeSEfromGeneExpressionQuantification(df, assay.list, genome = genome)
  } else {
    rownames(df) <- df$gene_id
    df$gene_id <- NULL
  }
  return(df)
}


makeSEfromGeneExpressionQuantification <- function(
    df,
    assay.list,
    genome = "hg19"
){
  
  # Access genome information to create SE
  gene.location <- get.GRCh.bioMart(genome)
  
  if(all(grepl("\\|",df[[1]]))){
    aux <- strsplit(df$gene_id,"\\|")
    GeneID <- unlist(lapply(aux,function(x) x[2]))
    df$entrezgene_id <- as.numeric(GeneID)
    gene.location <- gene.location[!duplicated(gene.location$entrezgene_id),]
    df <- merge(df, gene.location, by = "entrezgene_id")
  } else {
    df$external_gene_name <- as.character(df[[1]])
    df <- merge(df, gene.location, by = "external_gene_name")
  }
  
  
  if("transcript_id" %in% assay.list){
    rowRanges <- GRanges(
      seqnames = paste0("chr", df$chromosome_name),
      ranges = IRanges(start = df$start_position,
                       end = df$end_position),
      strand = df$strand,
      gene_id = df$external_gene_name,
      entrezgene = df$entrezgene_id,
      ensembl_gene_id = df$ensembl_gene_id,
      transcript_id = subset(df, select = 5)
    )
    names(rowRanges) <- as.character(df$gene_id)
    assay.list <- assay.list[which(assay.list != "transcript_id")]
  } else {
    rowRanges <- GRanges(
      seqnames = paste0("chr", df$chromosome_name),
      ranges = IRanges(
        start = df$start_position,
        end = df$end_position
      ),
      strand = df$strand,
      gene_id = df$external_gene_name,
      entrezgene = df$entrezgene_id,
      ensembl_gene_id = df$ensembl_gene_id
    )
    names(rowRanges) <- as.character(df$external_gene_name)
  }
  
  suppressWarnings({
    assays <- lapply(assay.list, function(x) {
      return(
        data.matrix(
          subset(df, select = grep(x,colnames(df),ignore.case = TRUE))
        )
      )
    })
  })
  
  names(assays) <- assay.list
  regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                  "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
  samples <- na.omit(unique(str_match(colnames(df),regex)[,1]))
  colData <-  colDataPrepare(samples)
  
  assays <- lapply(assays, function(x){
    colnames(x) <- NULL
    rownames(x) <- NULL
    return(x)
  })
  
  message("Available assays in SummarizedExperiment : \n  => ",paste(names(assays), collapse = "\n  => "))
  rse <- SummarizedExperiment(
    assays = assays,
    rowRanges = rowRanges,
    colData = colData
  )
  return(rse)
}


colDataPrepareTARGET <- function(barcode){
  message("Adding description to TARGET samples")
  tissue.code <- c(
    '01',
    '02',
    '03',
    '04',
    '05',
    '06',
    '07',
    '08',
    '09',
    '10',
    '11',
    '12',
    '13',
    '14',
    '15',
    '16',
    '17',
    '20',
    '40',
    '41',
    '42',
    '50',
    '60',
    '61',
    '99'
  )
  
  definition <- c(
    "Primary solid Tumor", # 01
    "Recurrent Solid Tumor", # 02
    "Primary Blood Derived Cancer - Peripheral Blood", # 03
    "Recurrent Blood Derived Cancer - Bone Marrow", # 04
    "Additional - New Primary", # 05
    "Metastatic", # 06
    "Additional Metastatic", # 07
    "Tissue disease-specific post-adjuvant therapy", # 08
    "Primary Blood Derived Cancer - Bone Marrow", # 09
    "Blood Derived Normal", # 10
    "Solid Tissue Normal",  # 11
    "Buccal Cell Normal",   # 12
    "EBV Immortalized Normal", # 13
    "Bone Marrow Normal", # 14
    "Fibroblasts from Bone Marrow Normal", # 15
    "Mononuclear Cells from Bone Marrow Normal", # 16
    "Lymphatic Tissue Normal (including centroblasts)", # 17
    "Control Analyte", # 20
    "Recurrent Blood Derived Cancer - Peripheral Blood", # 40
    "Blood Derived Cancer- Bone Marrow, Post-treatment", # 41
    "Blood Derived Cancer- Peripheral Blood, Post-treatment", # 42
    "Cell line from patient tumor", # 50
    "Xenograft from patient not grown as intermediate on plastic tissue culture dish", # 60
    "Xenograft grown in mice from established cell lines", #61
    "Granulocytes after a Ficoll separation"
  ) # 99
  aux <- DataFrame(tissue.code = tissue.code,definition)
  
  # in case multiple equal barcode
  regex <- paste0("[:alnum:]{6}-[:alnum:]{2}-[:alnum:]{6}",
                  "-[:alnum:]{3}-[:alnum:]{3}")
  samples <- str_match(barcode,regex)[,1]
  
  samples.df <- barcode %>% as.data.frame %>%  tidyr::separate(".",1:5 %>% as.character)
  
  ret <- DataFrame(
    barcode = barcode,
    tumor.code = samples.df[,2],
    sample = paste0(samples.df[,1],"-",samples.df[,2],"-",samples.df[,3],"-",samples.df[,4]),
    patient = paste0(samples.df[,1],"-",samples.df[,2],"-",samples.df[,3]),
    case.unique.id = paste0(samples.df[,3]),
    tissue.code = substr(samples.df[,4], 1, 2),
    nucleic.acid.code = substr(samples.df[,5], 3, 3)
  )
  
  ret <- merge(ret,aux, by = "tissue.code", sort = FALSE)
  
  tumor.code <- c('00','01','02','03','04','10','15','20','21','30','40',
                  '41','50','51','52','60','61','62','63','64','65','70','71','80','81')
  
  tumor.definition <- c(
    "Non-cancerous tissue", # 00
    "Diffuse Large B-Cell Lymphoma (DLBCL)", # 01
    "Lung Cancer (all types)", # 02
    "Cervical Cancer (all types)", # 03
    "Anal Cancer (all types)", # 04
    "Acute lymphoblastic leukemia (ALL)", # 10
    "Mixed phenotype acute leukemia (MPAL)", # 15
    "Acute myeloid leukemia (AML)", # 20
    "Induction Failure AML (AML-IF)", # 21
    "Neuroblastoma (NBL)", # 30
    "Osteosarcoma (OS)",  # 40
    "Ewing sarcoma",   # 41
    "Wilms tumor (WT)", # 50
    "Clear cell sarcoma of the kidney (CCSK)", # 51
    "Rhabdoid tumor (kidney) (RT)", # 52
    "CNS, ependymoma", # 60
    "CNS, glioblastoma (GBM)", # 61
    "CNS, rhabdoid tumor", # 62
    "CNS, low grade glioma (LGG)", # 63
    "CNS, medulloblastoma", # 64
    "CNS, other", # 65
    "NHL, anaplastic large cell lymphoma", # 70
    "NHL, Burkitt lymphoma (BL)", # 71
    "Rhabdomyosarcoma", #80
    "Soft tissue sarcoma, non-rhabdomyosarcoma"
  ) # 81
  aux <- DataFrame(tumor.code = tumor.code,tumor.definition)
  ret <- merge(ret,aux, by = "tumor.code", sort = FALSE)
  
  nucleic.acid.code <- c('D','E','W','X','Y','R','S')
  nucleic.acid.description <-  c(
    "DNA, unamplified, from the first isolation of a tissue",
    "DNA, unamplified, from the first isolation of a tissue embedded in FFPE",
    "DNA, whole genome amplified by Qiagen (one independent reaction)",
    "DNA, whole genome amplified by Qiagen (a second, separate independent reaction)",
    "DNA, whole genome amplified by Qiagen (pool of 'W' and 'X' aliquots)",
    "RNA, from the first isolation of a tissue",
    "RNA, from the first isolation of a tissue embedded in FFPE"
  )
  aux <- DataFrame(nucleic.acid.code = nucleic.acid.code,nucleic.acid.description)
  ret <- merge(ret,aux, by = "nucleic.acid.code", sort = FALSE)
  
  
  ret <- ret[match(barcode,ret$barcode),]
  rownames(ret) <- gsub("\\.","-",make.names(ret$barcode,unique=TRUE))
  ret$code <- NULL
  return(DataFrame(ret))
}

colDataPrepareBLGSP <- function(barcode){
  message("Adding description to BLGSP samples")
  tissue.code <- c(
    '01',
    '02',
    '03',
    '04',
    '05',
    '06',
    '07',
    '08',
    '09',
    '10',
    '11',
    '12',
    '13',
    '14',
    '15',
    '16',
    '17',
    '20',
    '40',
    '41',
    '42',
    '50',
    '60',
    '61',
    '99'
  )
  
  definition <- c(
    "Primary solid Tumor", # 01
    "Recurrent Solid Tumor", # 02
    "Primary Blood Derived Cancer - Peripheral Blood", # 03
    "Recurrent Blood Derived Cancer - Bone Marrow", # 04
    "Additional - New Primary", # 05
    "Metastatic", # 06
    "Additional Metastatic", # 07
    "Tissue disease-specific post-adjuvant therapy", # 08
    "Primary Blood Derived Cancer - Bone Marrow", # 09
    "Blood Derived Normal", # 10
    "Solid Tissue Normal",  # 11
    "Buccal Cell Normal",   # 12
    "EBV Immortalized Normal", # 13
    "Bone Marrow Normal", # 14
    "Fibroblasts from Bone Marrow Normal", # 15
    "Mononuclear Cells from Bone Marrow Normal", # 16
    "Lymphatic Tissue Normal (including centroblasts)", # 17
    "Control Analyte", # 20
    "Recurrent Blood Derived Cancer - Peripheral Blood", # 40
    "Blood Derived Cancer- Bone Marrow, Post-treatment", # 41
    "Blood Derived Cancer- Peripheral Blood, Post-treatment", # 42
    "Cell line from patient tumor", # 50
    "Xenograft from patient not grown as intermediate on plastic tissue culture dish", # 60
    "Xenograft grown in mice from established cell lines", #61
    "Granulocytes after a Ficoll separation"
  ) # 99
  aux <- DataFrame(tissue.code = tissue.code,definition)
  
  # in case multiple equal barcode
  regex <- paste0("[:alnum:]{5}-[:alnum:]{2}-[:alnum:]{2}",
                  "-[:alnum:]{5}-[:alnum:]{3}")
  samples <- str_match(barcode,regex)[,1]
  
  samples.df <- barcode %>% as.data.frame %>%  tidyr::separate(".",1:5 %>% as.character)
  
  ret <- DataFrame(
    barcode = barcode,
    tumor.code = samples.df[,2],
    sample = paste0(samples.df[,1],"-",samples.df[,2],"-",samples.df[,3],"-",samples.df[,4]),
    patient = paste0(samples.df[,1],"-",samples.df[,2],"-",samples.df[,3]),
    case.unique.id = paste0(samples.df[,3]),
    tissue.code = substr(samples.df[,4], 1, 2),
    nucleic.acid.code = substr(samples.df[,5], 3, 3)
  )
  
  ret <- merge(ret,aux, by = "tissue.code", sort = FALSE)
  
  tumor.code <- c('00','01','02','03','04','10','15','20','21','30','40',
                  '41','50','51','52','60','61','62','63','64','65','70','71','80','81')
  
  tumor.definition <- c(
    "Non-cancerous tissue", # 00
    "Diffuse Large B-Cell Lymphoma (DLBCL)", # 01
    "Lung Cancer (all types)", # 02
    "Cervical Cancer (all types)", # 03
    "Anal Cancer (all types)", # 04
    "Acute lymphoblastic leukemia (ALL)", # 10
    "Mixed phenotype acute leukemia (MPAL)", # 15
    "Acute myeloid leukemia (AML)", # 20
    "Induction Failure AML (AML-IF)", # 21
    "Neuroblastoma (NBL)", # 30
    "Osteosarcoma (OS)",  # 40
    "Ewing sarcoma",   # 41
    "Wilms tumor (WT)", # 50
    "Clear cell sarcoma of the kidney (CCSK)", # 51
    "Rhabdoid tumor (kidney) (RT)", # 52
    "CNS, ependymoma", # 60
    "CNS, glioblastoma (GBM)", # 61
    "CNS, rhabdoid tumor", # 62
    "CNS, low grade glioma (LGG)", # 63
    "CNS, medulloblastoma", # 64
    "CNS, other", # 65
    "NHL, anaplastic large cell lymphoma", # 70
    "NHL, Burkitt lymphoma (BL)", # 71
    "Rhabdomyosarcoma", #80
    "Soft tissue sarcoma, non-rhabdomyosarcoma"
  ) # 81
  aux <- DataFrame(tumor.code = tumor.code,tumor.definition)
  ret <- merge(ret,aux, by = "tumor.code", sort = FALSE)
  
  # nucleic.acid.code <- c('D','E','W','X','Y','R','S')
  # nucleic.acid.description <-  c(
  #   "DNA, unamplified, from the first isolation of a tissue",
  #   "DNA, unamplified, from the first isolation of a tissue embedded in FFPE",
  #   "DNA, whole genome amplified by Qiagen (one independent reaction)",
  #   "DNA, whole genome amplified by Qiagen (a second, separate independent reaction)",
  #   "DNA, whole genome amplified by Qiagen (pool of 'W' and 'X' aliquots)",
  #   "RNA, from the first isolation of a tissue",
  #   "RNA, from the first isolation of a tissue embedded in FFPE"
  # )
  # aux <- DataFrame(nucleic.acid.code = nucleic.acid.code,nucleic.acid.description)
  # ret <- merge(ret,aux, by = "nucleic.acid.code", sort = FALSE)
  
  
  ret <- ret[match(barcode,ret$barcode),]
  rownames(ret) <- gsub("\\.","-",make.names(ret$barcode,unique=TRUE))
  ret$code <- NULL
  return(DataFrame(ret))
}

colDataPrepareOHSU <- function(barcode){
  message("Adding description to OHSU samples")
  
  browser()
  samples <- str_replace_all(string = barcode, pattern = "R", replacement = "-R")
  
  samples.df <- samples %>% as.data.frame %>%  tidyr::separate(".",1:2 %>% as.character)
  
  ret <- DataFrame(
    barcode = barcode,
    sample = paste0(samples.df[,1])
  )
  
  
  ret <- ret[match(barcode,ret$barcode),]
  rownames(ret) <- ret$barcode
  ret$code <- NULL
  return(DataFrame(ret))
}

colDataPrepare <- function(barcode){
  # For the moment this will work only for TCGA Data
  # We should search what TARGET data means
  message("Starting to add information to samples")
  ret <- NULL
  browser()
  if(all(grepl("TARGET",barcode))) ret <- colDataPrepareTARGET(barcode)
  if(all(grepl("BLGSP",barcode))) ret <- colDataPrepareBLGSP(barcode)
  if(all(grepl("\\d{4}R",barcode))) ret <- colDataPrepareOHSU(barcode)
  if(all(grepl("TCGA",barcode))) ret <- colDataPrepareTCGA(barcode)
  if(all(grepl("MMRF",barcode))) ret <- colDataPrepareMMRF(barcode)
  
  # How to deal with mixed samples "C3N-02003-01;C3N-02003-021" ?
  # Check if this breaks the pacakge
  if(any(grepl("C3N-|C3L-",barcode))) {
    ret <- data.frame(
      sample =  sapply(barcode, function(x) stringr::str_split(x,";") %>% unlist()) %>%
        unlist %>% unique,stringsAsFactors = FALSE
    )
  }
  if(is.null(ret)) {
    ret <- data.frame(
      sample = barcode %>% unique,
      stringsAsFactors = FALSE
    )
  }
  
#modified code - Clinical data addition REMOVED - doesn't work for TARGET-AML
  
  # message(" => Add clinical information to samples")
  # # There is a limitation on the size of the string, so this step will be splited in cases of 100
  # patient.info <- NULL
  # patient.info <- splitAPICall(
  #   FUN = getBarcodeInfo,
  #   step = 10,
  #   items = ret$sample
  # )
  # 
  # if(!is.null(patient.info)) {
  #   ret$sample_submitter_id <- ret$sample %>% as.character()
  #   ret <- left_join(ret %>% as.data.frame, patient.info %>% unique, by = "sample_submitter_id")
  # }
  # ret$bcr_patient_barcode <- ret$sample %>% as.character()
  # ret$sample_submitter_id <- ret$sample %>% as.character()
  # 
  # if(!"project_id" %in% colnames(ret)) {
  #   if("disease_type" %in% colnames(ret)){
  #     aux <- getGDCprojects()[,c(5,7)]
  #     aux <- aux[aux$disease_type == unique(ret$disease_type),2]
  #     ret$project_id <- as.character(aux)
  #   }
  # }
  # # There is no subtype info for target, return as it is
  # if(any(grepl("TCGA",barcode))) {
  #   ret <- addSubtypeInfo(ret)
  # }
  # 
  # # na.omit should not be here, exceptional case
  # if(is.null(ret)) return(data.frame(row.names = barcode, barcode,stringsAsFactors = FALSE))
  # 
  # # Add purity information from http://www.nature.com/articles/ncomms9971
  # # purity  <- getPurityinfo()
  # # ret <- merge(ret, purity, by = "sample", all.x = TRUE, sort = FALSE)
  # 
  # # Put data in the right order
  # ret <- ret[!duplicated(ret$bcr_patient_barcode),]
  # 
  # 
  # idx <- sapply(substr(barcode,1,min(str_length(ret$bcr_patient_barcode))), function(x) {
  #   grep(x,ret$bcr_patient_barcode)
  # })
  # 
  # # the code above does not work, since the projects have different string lengths
  # if(all(ret$project_id == "TARGET-ALL-P3")) {
  #   idx <- sapply(gsub("-[[:alnum:]]{3}$","",barcode), function(x) {
  #     grep(x,ret$bcr_patient_barcode)
  #   })
  # }
  # 
  # if(any(ret$project_id == "CPTAC-3",na.rm = T)) {
  #   idx <- sapply(gsub("-[[:alnum:]]{3}$","",barcode), function(x) {
  #     if(grepl(";",x = x)) x <- stringr::str_split(barcode[1],";")[[1]][1] # mixed samples
  #     grep(x,ret$bcr_patient_barcode)
  #   })
  # }
  # 
  # if(any(ret$project_id == "CMI-MBC",na.rm = T)) {
  #   idx <- match(barcode,ret$bcr_patient_barcode)
  # }
  # 
  # ret <- ret[idx,]
  # 
  if("barcode" %in% colnames(ret)) ret$barcode <- barcode
  
  rownames(ret) <- barcode
  return(ret)
}



makeSEfromTranscriptomeProfilingSTAR <- function(data, cases, assay.list){
  # How many cases do we have?
  # We wil consider col 1 is the ensemble gene id, other ones are data
  size <- ncol(data)
  # Prepare Patient table
  colData <-  colDataPrepare(cases)
  
  # one ensemblID can be mapped to multiple entrezgene ID
  gene.location <- get.GRCh.bioMart("hg38",as.granges = TRUE)
  
  metrics <- subset(data, !grepl("ENSG", data$gene_id))
  data <- subset(data, grepl("ENSG", data$gene_id))
  found.genes <- table(data$gene_id %in% gene.location$gene_id)
  
  if("FALSE" %in% names(found.genes))
    message(paste0("From the ", nrow(data), " genes we couldn't map ", found.genes[["FALSE"]]))
  
  # Prepare data table
  # Remove the version from the ensembl gene id
  assays <- list(
    data.matrix(data[,grep("^unstranded",colnames(data))]),
    data.matrix(data[,grep("stranded_first",colnames(data))]),
    data.matrix(data[,grep("stranded_second",colnames(data))]),
    data.matrix(data[,grep("fpkm_unstrand",colnames(data))]),
    data.matrix(data[,grep("tpm_unstrand",colnames(data))])
  )
  names(assays) <- c("unstranded","stranded_first","stranded_second","fpkm_unstrand","tpm_unstrand")
  assays <- lapply(assays, function(x){
    colnames(x) <- NULL
    rownames(x) <- NULL
    return(x)
  })
  
  # Prepare rowRanges
  rowRanges <- gene.location[match(data$gene_id, gene.location$gene_id),]
  names(rowRanges) <- as.character(data$gene_id)
  
  
  message("Available assays in SummarizedExperiment : \n  => ",paste(names(assays), collapse = "\n  => "))
  rse <- SummarizedExperiment(
    assays = assays,
    rowRanges = rowRanges,
    colData = colData
  )
  metadata(rse) <- metrics
  return(rse)
}


makeSEfromTranscriptomeProfiling <- function(data, cases, assay.list){
  # How many cases do we have?
  # We wil consider col 1 is the ensemble gene id, other ones are data
  size <- ncol(data)
  # Prepare Patient table
  colData <-  colDataPrepare(cases)
  
  # one ensemblID can be mapped to multiple entrezgene ID
  gene.location <- get.GRCh.bioMart("hg38")
  gene.location <- gene.location[!duplicated(gene.location$ensembl_gene_id),]
  
  data$ensembl_gene_id <- as.character(gsub("\\.[0-9]*","",data$X1))
  data <- subset(data, grepl("ENSG", data$ensembl_gene_id))
  found.genes <- table(data$ensembl_gene_id %in% gene.location$ensembl_gene_id)
  
  if("FALSE" %in% names(found.genes))
    message(paste0("From the ", nrow(data), " genes we couldn't map ", found.genes[["FALSE"]]))
  
  data <- merge(data, gene.location, by = "ensembl_gene_id")
  
  # Prepare data table
  # Remove the version from the ensembl gene id
  assays <- list(data.matrix(data[,2:size+1]))
  names(assays) <- assay.list
  assays <- lapply(assays, function(x){
    colnames(x) <- NULL
    rownames(x) <- NULL
    return(x)
  })
  
  # Prepare rowRanges
  rowRanges <- GRanges(
    seqnames = paste0("chr", data$chromosome_name),
    ranges = IRanges(
      start = data$start_position,
      end = data$end_position
    ),
    strand = data$strand,
    ensembl_gene_id = data$ensembl_gene_id,
    external_gene_name = data$external_gene_name,
    original_ensembl_gene_id = data$X1
  )
  names(rowRanges) <- as.character(data$ensembl_gene_id)
  message("Available assays in SummarizedExperiment : \n  => ",paste(names(assays), collapse = "\n  => "))
  rse <- SummarizedExperiment(
    assays = assays,
    rowRanges = rowRanges,
    colData = colData
  )
  
  return(rse)
}


#' @importFrom dplyr left_join
#' @importFrom plyr alply join_all
#' @importFrom purrr map_dfc map map_df
readTranscriptomeProfiling <- function(
    files,
    data.type,
    workflow.type,
    cases,
    summarizedExperiment
) {
  
  if(grepl("Gene Expression Quantification", data.type, ignore.case = TRUE)){
    # Status working for:
    #  - htseq
    #  - FPKM
    #  - FPKM-UQ
    if(grepl("HTSeq",workflow.type)){
      
      x <- plyr::alply(files,1, function(f) {
        readr::read_tsv(
          file = f,
          col_names = FALSE,
          progress = FALSE,
          col_types = c("cd")
        )
      }, .progress = "time")
      
      # Just check if the data is in the same order, since we will not merge
      # the data frames to save memory
      stopifnot(all(unlist(x %>% map(function(y){all(y[,1] ==  x[[1]][,1])}) )))
      
      # need to check if it works in all cases
      # tested for HTSeq - FPKM-UQ and Counts only
      df <- bind_cols(x[[1]][,1],x %>%  map_df(2))
      if(!missing(cases))  colnames(df)[-1] <- cases
      if(summarizedExperiment) df <- makeSEfromTranscriptomeProfiling(df,cases,workflow.type)
    } else  if(grepl("STAR",workflow.type)){
      
      # read files that has 4 not necessary rows, and has several columns
      # gene_id gene_name gene_type
      # unstranded stranded_first stranded_second tpm_unstranded fpkm_unstranded
      x <- plyr::alply(files,1, function(f) {
        readr::read_tsv(
          file = f,
          col_names = TRUE,
          progress = FALSE,
          show_col_types = FALSE,
          skip = 1
        )
      }, .progress = "time")
      
      # bind all counts and then add the gene metadata
      suppressMessages({
        df <- x %>%  map_dfc(.f = function(y) y[,4:8])
        df <- bind_cols(x[[1]][,1:3],df)
      })
      
      # Adding barcodes to columns names, if user wants a dataframe
      if(!missing(cases))  {
        colnames(df)[-c(1:3)] <- sapply(cases, function(x){
          stringr::str_c(c("unstranded_","stranded_first_", "stranded_second_","tpm_unstranded_","fpkm_unstranded_"),x)}
        ) %>% as.character()
      }
      if(summarizedExperiment) df <- makeSEfromTranscriptomeProfilingSTAR(df,cases,workflow.type)
    }
  } else if(grepl("miRNA", workflow.type, ignore.case = TRUE) & grepl("miRNA", data.type, ignore.case = TRUE)) {
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for (i in seq_along(files)) {
      data <- read_tsv(file = files[i], col_names = TRUE,col_types = "cidc")
      if(!missing(cases))
        colnames(data)[2:ncol(data)] <- paste0(colnames(data)[2:ncol(data)],"_",cases[i])
      
      if(i == 1) df <- data
      if(i != 1) df <- merge(df, data, by=colnames(df)[1],all = TRUE)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
  } else if(grepl("Isoform Expression Quantification", data.type, ignore.case = TRUE)){
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for (i in seq_along(files)) {
      data <- read_tsv(file = files[i], col_names = TRUE, col_types = c("ccidcc"))
      if(!missing(cases)) data$barcode <- cases[i] else data$file <- i
      if(i == 1) df <- data
      if(i != 1) df <- rbind(df,data)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  return(df)
}

readGeneLevelCopyNumber <- function(files, cases, summarizedExperiment = FALSE){
  message("Reading Gene Level Copy Number files")
  gistic.df <- NULL
  gistic.list <- plyr::alply(files,1,.fun = function(file) {
    #message("Reading file: ", file)
    data <- read_tsv(
      file = file,
      col_names = TRUE,
      progress = TRUE,
      col_types = readr::cols()
    )
    colnames(data)[-c(1:5)] <- paste0(cases[match(file,files)],"_",  colnames(data)[-c(1:5)])
    
    return(data)
  })
  
  # Check if the data is in the same order
  stopifnot(all(unlist(gistic.list %>% map(function(y){all(y[,1:5] ==  gistic.list[[1]][,1:5])}) )))
  
  # need to check if it works in all cases
  # tested for HTSeq - FPKM-UQ and Counts only
  df <- bind_cols(
    gistic.list[[1]][,1:5], # gene info
    gistic.list %>%  map_dfc(.f = function(x) x[,6:8]) # copy number, min,max
  )
  
  if(summarizedExperiment) {
    se <- makeSEfromGeneLevelCopyNumber(df, cases)
    return(se)
  }
  return(df)
}

makeSEfromGeneLevelCopyNumber <- function(df, cases){
  message("Creating a SummarizedExperiment object")
  rowRanges <- GRanges(
    seqnames = df$chromosome,
    ranges = IRanges(start = df$start, end = df$end),
    gene_id = df$gene_id,
    gene_name = df$gene_name
  )
  
  names(rowRanges) <- as.character(df$gene_id)
  colData <-  colDataPrepare(cases)
  
  # We have three data columns for each files
  assays <- lapply(
    c("[^max|min]_copy_number","min_copy_number", "max_copy_number"),
    function(x){
      ret <- data.matrix(subset(df,select = grep(x,colnames(df))))
      colnames(ret) <- cases
      rownames(ret) <- NULL
      ret
    })
  names(assays) <- c("copy_number","min_copy_number", "max_copy_number")
  
  message("Available assays in SummarizedExperiment : \n  => ",paste(names(assays), collapse = "\n  => "))
  rse <- SummarizedExperiment(assays = assays, rowRanges = rowRanges, colData = colData)
  return(rse)
}


readGISTIC <- function(files, cases){
  message("Reading GISTIC file")
  gistic.df <- NULL
  gistic.list <- plyr::alply(files,1,.fun = function(file) {
    message("Reading file: ", file)
    data <- read_tsv(
      file = file,
      col_names = TRUE,
      progress = TRUE,
      col_types = readr::cols()
    )
    
    aliquot <- colnames(data)[-c(1:3)]
    info <- splitAPICall(
      FUN = getBarcodefromAliquot,
      step = 20,
      items = aliquot
    )
    
    barcode <- as.character(info$submitter_id)[match(aliquot,as.character(info$aliquot_id))]
    colnames(data)[-c(1:3)] <- barcode
    return(data)
  })
  gistic.df <- gistic.list %>%
    join_all(by =  c("Gene Symbol","Gene ID","Cytoband"), type='full')
  
  return(gistic.df)
}

# Reads Copy Number Variation files to a data frame, basically it will rbind it
#' @importFrom purrr map2_dfr
readCopyNumberVariation <- function(files, cases){
  message("Reading copy number variation files")
  
  col_types <- ifelse(any(grepl('ascat2', files)),"ccnnnnn","ccnnnd")
  
  purrr::map2_dfr(
    .x = files,
    .y = cases,
    .f = function(file,case) {
      data <- readr::read_tsv(file, col_names = TRUE, col_types = col_types);
      if(!missing(case)) data$Sample <- case
      data
    })
}

# getBarcodeInfo(c("TCGA-A6-6650-01B"))
addFFPE <- function(df) {
  message("Add FFPE information. More information at: \n=> https://cancergenome.nih.gov/cancersselected/biospeccriteria \n=> http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html")
  barcode <- df$barcode
  patient <- df$patient
  
  ffpe.info <- NULL
  ffpe.info <- splitAPICall(FUN = getFFPE,
                            step = 20,
                            items = df$patient)
  
  df <- merge(df, ffpe.info,by.x = "sample", by.y = "submitter_id")
  df <- df[match(barcode,df$barcode),]
  return(df)
}

# getFFPE("TCGA-A6-6650")
#' @importFrom plyr rbind.fill
#' @importFrom httr content
getFFPE <- function(patient){
  baseURL <- "https://api.gdc.cancer.gov/cases/?"
  options.pretty <- "pretty=true"
  options.expand <- "expand=samples"
  option.size <- paste0("size=",length(patient))
  options.filter <- paste0("filters=",
                           URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.submitter_id","value":['),
                           paste0('"',paste(patient,collapse = '","')),
                           URLencode('"]}}]}'))
  url <- paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&"))
  json  <- tryCatch(
    getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
    error = function(e) {
      message(paste("Error: ", e, sep = " "))
      message("We will retry to access GDC again! URL:")
      #message(url)
      fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
    }
  )
  
  results <- json$data$hits
  results <- rbind.fill(results$samples)[,c("submitter_id","is_ffpe")]
  return(results)
}

getAliquot_ids <- function(barcode){
  baseURL <- "https://api.gdc.cancer.gov/cases/?"
  options.fields <- "fields=samples.portions.analytes.aliquots.aliquot_id,samples.portions.analytes.aliquots.submitter_id"
  options.pretty <- "pretty=true"
  option.size <- paste0("size=",length(barcode))
  #message(paste(barcode,collapse = '","'))
  #message(paste0('"',paste(barcode,collapse = '","')))
  options.filter <- paste0("filters=",
                           URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.submitter_id","value":['),
                           paste0('"',paste(barcode,collapse = '","')),
                           URLencode('"]}}]}'))
  #message(paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&")))
  url <- paste0(baseURL,paste(options.pretty,options.fields, option.size, options.filter, sep = "&"))
  json  <- tryCatch(
    getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
    error = function(e) {
      message(paste("Error: ", e, sep = " "))
      message("We will retry to access GDC again! URL:")
      #message(url)
      fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
    }
  )
  results <- unlist(json$data$hits$samples)
  results.barcode <- grep("TCGA",results,value = TRUE)
  results.aliquot <- grep("TCGA",results,value = TRUE,invert = TRUE)
  
  df <- data.frame(results.aliquot,results.barcode)
  colnames(df) <- c("aliquot_id","barcode")
  return(df)
}


# getBarcodeInfo(c("TCGA-OR-A5K3-01A","C3N-00321-01"))
# barcode is: sample_submitter_id
#' @importFrom dplyr bind_cols slice row_number
#'
getBarcodeInfo <- function(barcode) {
  baseURL <- "https://api.gdc.cancer.gov/cases/?"
  options.pretty <- "pretty=true"
  options.expand <- "expand=samples,project,diagnoses,diagnoses.treatments,annotations,family_histories,demographic,exposures"
  option.size <- paste0("size=",length(barcode))
  options.filter <- paste0("filters=",
                           URLencode('{"op":"or","content":[{"op":"in","content":{"field":"cases.submitter_id","value":['),
                           paste0('"',paste(barcode,collapse = '","')),
                           URLencode('"]}},'),
                           URLencode('{"op":"in","content":{"field":"submitter_sample_ids","value":['),
                           paste0('"',paste(barcode,collapse = '","')),
                           URLencode('"]}},'),
                           URLencode('{"op":"in","content":{"field":"submitter_aliquot_ids","value":['),
                           paste0('"',paste(barcode,collapse = '","')),
                           URLencode('"]}}'),
                           URLencode(']}')
  )
  url <- paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&"))
  #message(url)
  json  <- tryCatch(
    getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
    error = function(e) {
      message(paste("Error: ", e, sep = " "))
      message("We will retry to access GDC again! URL:")
      #message(url)
      fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
    }
  )
  
  results <- json$data$hits
  
  # no results
  if(length(results) == 0){
    return(data.frame(barcode,stringsAsFactors = FALSE))
  }
  
  submitter_id <- results$submitter_id
  submitter_aliquot_ids <- results$submitter_aliquot_ids
  
  if(!is.null(results$samples)) {
    samples <- rbindlist(results$samples, fill = TRUE)
    samples <- samples[match(barcode,samples$submitter_id),]
    samples$sample_submitter_id <- str_extract_all(samples$submitter_id,paste(barcode,collapse = "|")) %>%
      unlist %>% as.character
    
    tryCatch({
      samples$submitter_id <-
        str_extract_all(samples$submitter_id,
                        paste(c(submitter_id,barcode), collapse = "|"),
                        simplify = TRUE) %>% as.character
    }, error = function(e){
      samples$submitter_id <- submitter_id
    })
    df <- samples[!is.na(samples$submitter_id),]
    suppressWarnings({
      df[,c("updated_datetime","created_datetime")] <- NULL
    })
  }
  
  
  # We dont have the same cols for TCGA and TARGET so we need to check them
  if(!is.null(results$diagnoses)) {
    diagnoses <- rbindlist(lapply(results$diagnoses, function(x) if(is.null(x)) data.frame(NA) else x),fill = TRUE)
    diagnoses[,c("updated_datetime","created_datetime","state")] <- NULL
    if(any(grepl("submitter_id", colnames(diagnoses)))) {
      diagnoses$submitter_id <- gsub("_diagnosis.*|-DIAG|diag-","", diagnoses$submitter_id)
    }  else {
      diagnoses$submitter_id <- submitter_id
    }
    # this is required since the sample might not have a diagnosis
    if(!any(df$submitter_id %in% diagnoses$submitter_id)){
      diagnoses$submitter_id <- NULL
      # The sample migth have different sample types
      # The diagnosis the same for each one of these samples
      # in that case we will have a 1 to mapping and binding will
      # not work. We need then to replicate diagnosis to each sample
      # and not each patient
      # Cases can be replicated with getBarcodeInfo(c("BA2691R","BA2577R","BA2748R"))
      if(nrow(diagnoses) <  nrow(df)){
        diagnoses <- plyr::ldply(
          1:length(results$submitter_sample_ids),
          .fun = function(x){
            diagnoses[x] %>% # replicate diagnoses the number of samples
              as.data.frame() %>%
              dplyr::slice(rep(dplyr::row_number(), sum(results$submitter_sample_ids[[x]] %in% barcode)))})
      }
      
      df <- dplyr::bind_cols(df %>% as.data.frame,diagnoses %>% as.data.frame)
    } else {
      df <- left_join(df, diagnoses, by = "submitter_id")
    }
  }
  if(!is.null(results$exposures)) {
    exposures <- rbindlist(lapply(results$exposures, function(x) if(is.null(x)) data.frame(NA) else x),fill = TRUE)
    exposures[,c("updated_datetime","created_datetime","state")] <- NULL
    if(any(grepl("submitter_id", colnames(exposures)))) {
      exposures$submitter_id <- gsub("_exposure.*|-EXP","", exposures$submitter_id)
    }  else {
      exposures$submitter_id <- submitter_id
    }
    if(!any(df$submitter_id %in% exposures$submitter_id)){
      exposures$submitter_id <- NULL
      df <- dplyr::bind_cols(df,exposures)
    } else {
      df <- left_join(df, exposures, by = "submitter_id")
    }
  }
  
  
  if(!is.null(results$demographic)) {
    demographic <- results$demographic
    demographic[,c("updated_datetime","created_datetime","state")] <- NULL
    if(any(grepl("submitter_id", colnames(demographic)))) {
      demographic$submitter_id <- gsub("_demographic.*|-DEMO|demo-","", results$demographic$submitter_id)
    } else {
      demographic$submitter_id <- submitter_id
    }
    
    if(!any(df$submitter_id %in% demographic$submitter_id)){
      demographic$submitter_id <- NULL
      demographic$updated_datetime  <- NULL
      demographic$created_datetime <- NULL
      # The sample migth have different sample types
      # The diagnosis the same for each one of these samples
      # in that case we will have a 1 to mapping and binding will
      # not work. We need then to replicate diagnosis to each sample
      # and not each patient
      # Cases can be replicated with getBarcodeInfo(c("BA2691R","BA2577R","BA2748R"))
      if(nrow(demographic) <  nrow(df)){
        demographic <- plyr::ldply(
          1:length(results$submitter_sample_ids),
          .fun = function(x){
            demographic[x,] %>% # replicate diagnoses the number of samples
              as.data.frame() %>%
              dplyr::slice(rep(dplyr::row_number(), sum(results$submitter_sample_ids[[x]] %in% barcode)))})
      }
      
      df <- dplyr::bind_cols(df  %>% as.data.frame,demographic)
    } else {
      df <- left_join(df,demographic, by = "submitter_id")
    }
  }
  
  df$bcr_patient_barcode <- df$submitter_id %>% as.character()
  projects.info <- results$project
  projects.info <- results$project[,grep("state",colnames(projects.info),invert = TRUE)]
  
  if(any(submitter_id %in% df$submitter_id)){
    projects.info <-  cbind("submitter_id" = submitter_id, projects.info)
    
    suppressWarnings({
      df <- left_join(df,
                      projects.info,
                      by = "submitter_id")
    })
  } else {
    
    if(nrow(projects.info) <  nrow(df)){
      projects.info <- plyr::ldply(
        1:length(results$submitter_sample_ids),
        .fun = function(x){
          projects.info[x,] %>% # replicate diagnoses the number of samples
            as.data.frame() %>%
            dplyr::slice(rep(dplyr::row_number(), sum(results$submitter_sample_ids[[x]] %in% barcode)))})
    }
    
    df <- dplyr::bind_cols(df,projects.info)
  }
  
  # Adding in the same order
  
  
  if(any(substr(barcode,1,str_length(df$submitter_id)) %in% df$submitter_id)){
    df <- df[match(substr(barcode,1,str_length(df$sample_submitter_id)),df$sample_submitter_id),]
    # This line should not exists, but some patients does not have clinical data
    # case: TCGA-R8-A6YH"
    # this has been reported to GDC, waiting answers
    # So we will remove this NA cases
    df <- df[!is.na(df$submitter_id),]
  } else {
    idx <- sapply(substr(barcode,1,str_length(df$submitter_aliquot_ids) %>% max),
                  FUN = function(x){
                    grep(x,df$submitter_aliquot_ids)
                  })
    df <- df[idx,]
  }
  
  # remove empty columns
  df <- df %>% as.data.frame() %>% dplyr::select(which(colSums(is.na(df)) < nrow(df)))
  return(df)
}

#' @title Prepare CEL files into an AffyBatch.
#' @description Prepare CEL files into an AffyBatch.
#' @param ClinData write
#' @param PathFolder write
#' @param TabCel write
#' @examples
#' \dontrun{
#' to add example
#' }
#' @export
#' @return Normalized Expression data from Affy eSets
TCGAprepare_Affy <- function(ClinData, PathFolder, TabCel){
  if (!requireNamespace("affy", quietly = TRUE)) {
    stop("affy package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("affy package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  affy_batch <- affy::ReadAffy(filenames = as.character(paste(TabCel$samples, ".CEL", sep = "")))
  
  eset <- affy::rma(affy_batch)
  
  mat <- Biobase::exprs(eset)
  
  return(mat)
  
}

splitAPICall <- function(FUN, step = 20, items){
  info <- NULL
  info <- tryCatch({
    for(i in 0:(ceiling(length(items)/step) - 1)){
      start <- 1 + step * i
      end <- ifelse(((i + 1) * step) > length(items), length(items),((i + 1) * step))
      if(is.null(info)) {
        info <- FUN(items[start:end])
      } else {
        info <- plyr::rbind.fill(info, FUN(items[start:end]))
      }
    }
    info
  }, error = function(e) {
    step <- 2
    for(i in 0:(ceiling(length(items)/step) - 1)){
      start <- 1 + step * i
      end <- ifelse(((i + 1) * step) > length(items), length(items),((i + 1) * step))
      if(is.null(info)) {
        info <- FUN(items[start:end])
      } else {
        info <- plyr::rbind.fill(info, FUN(items[start:end]))
      }
    }
  })
  info
}



library(purrr)

query <- GDCquery(
  project = "TARGET-AML",
  data.category = c("Transcriptome Profiling"),
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")

black_list <- c("236dda1b-5649-4824-8d32-6e411e84b590","a65bf2da-bcc2-488c-bba2-b955b9f5fe15",
                "2d7fa837-74ff-4f66-840f-1121e3c1fc84","32a7be9a-75de-4286-82a2-74c79e01f733",
                "fe7aea91-81f6-45d5-be51-ce5f1aa355b5","f0859c93-8d7f-404a-a58b-6a380153aff8",
                "abb0a6b6-9c14-40e5-8ffa-81af4d73fc0a","3f2e13af-b20c-4113-81c8-1f719239716f",
                "75a09478-e8db-4941-9ff0-9032121db459","b1154f54-8954-460f-be65-6c62756804b4",
                "116fc74b-510d-4feb-baf4-f9b0a74894fd","8efacaab-8d61-4225-837b-c42a1e4c914f",
                "a372f2c6-ef56-4861-8f56-71b494ba3c82","dddd0b69-d3e0-456b-953d-c1fd05ad9bc8",
                "2910ce5f-cce4-47fc-9cd1-29460c645a31","bcce450e-53c9-401c-9e3c-15726aee6b74",
                "3133ccaa-b4e4-48b3-8430-ea0f082c8a48","adf95118-0837-44a8-8b28-8626cbf227c7",
                "836222f1-6067-4fc6-bba1-a587a67d7e89","a12d6fab-9113-4223-9290-f8d9ac7c8ae9",
                "8e042700-bb72-4b94-90ca-6d67d48a96f0","b8f44d0a-6812-4a0a-b338-acd695283105",
                "9ec79d96-3b36-4f45-9d1a-f5966b9a6352","8f9f21f9-419b-4366-a106-9391a58162d2",
                "d014957c-bf98-4958-b890-413775a73cfa","bb3169af-8965-4b46-a3a4-14c7ee101f20",
                "ea380475-03ab-4f91-84ec-8d93e7371e65","39499910-9b5d-4d5a-9433-9feab1603157",
                "bfd657f8-16a8-419a-9329-ba4da1c075ab","bace435c-275a-411c-9ca3-4920b6970abb",
                "3a4d2e18-f553-4e24-bbc1-7732d19c49b1","ba2aefd4-7716-4857-bcd4-cc64ee981c30",
                "59ed7f0d-4697-4f76-aa7f-2193749fceaf","894b3ae0-af8c-427f-b831-d39b6ec051c3",
                "d7b835a4-43db-4858-ac06-7b40f4de7436","776a4b95-1a90-4db4-b83d-3b7aa862b148",
                "28fce5a1-a731-442c-a2a8-13d29c89e7ca","5991c0a7-375c-43c9-b85b-b936af5aa645",
                "c5a4b24c-1c23-40f7-b13a-ad4f8a70a876","89b9d47c-7b30-43ff-bc05-2010b2f0436e",
                "144ac8dd-9238-42a1-b4f0-5043770f2d27",
                "9c74b223-0bc4-4420-a4bd-cfe22988c5b4","6f7bfe53-6000-4d35-944c-a6be40b64642",
                "dabffad0-831a-40ae-afb5-804b9607f175","1e761c6f-2768-4aaa-afc3-9d08e583b239",
                "ce32c34c-b200-48ec-b2e7-ddfc20ea9d09","e04cf779-6fab-45a0-be8c-646e9c2f29f8",
                "d6cb2c7f-a27c-4a8e-8717-2d0aef8813d4","42fd6adc-1b02-429f-950a-9efe0b0d1695",
                "c0203a07-2fc6-42bc-a051-7bc6d896b156","72fd418b-79ca-4183-8072-a4aba1d3ffc4",
                "eb3a8090-e4d0-4a95-b1fa-122a220889a3","fccb5112-b630-42e5-80a6-661a6c48c66e",
                "335a5c6b-d228-47c6-a076-e16babe020f9","40d78d55-4249-4238-91a9-f928a625f342",
                "20eacea7-b494-4b73-9578-f98f92cf1a09","764f4696-4845-4f5c-8b8e-49ec2c76792f",
                "2f0de586-517e-4e49-adb1-d3293ecb606b","9605582a-eebd-48f1-b035-85d07e2283fc",
                "4fe07def-f108-43fd-b119-7f0537a801ce","7f914771-56e5-49e1-bf1a-31884de58829",
                "46db18ef-e39b-4a22-ae88-d1d3d7e07d70", "f8186338-e939-477f-a140-b27769949253",
                "a2c9c73f-832a-41fe-8ef2-b5b814f8eabb", "e22adaa4-03c1-4719-83c1-7f3d17e49f6e",
                "2b65893d-7d87-4826-875d-8865d28c79fc", "852222c6-52f1-4530-8071-77300b8d25a2",
                "7bcc4806-ec42-4d16-a5f8-108dd6f79d9d",
                "f9f6ed14-b993-4b1c-8bd9-8a8e9cc6348f","830b9273-d289-417b-ae52-2a477dc59cca",
                "13556ca2-7c66-47cd-bacd-7c63ec3f7e18","6da7c2bc-c47f-47db-9020-8f3bb753632f",
                "f2ca7571-9112-4f2c-8591-7028671bc0b9","cb85de6f-132d-4846-b467-11e906cd2041",
                "d5807d3b-3cec-4676-815e-d8f3b0a1b615","a5a831e1-df4b-44ba-a903-e8bcba06f0da",
                "727f5e27-09e9-425f-9757-81da66c054bd","e2f6e032-e16f-44ca-839c-a3e1cace79b7",
                "d9fdf6f2-3071-4b62-a3ec-b04071d13157","c6698c82-2293-48a9-8ba7-68295f027be6",
                "21953a49-b4c2-4c0f-b9e8-b69b7a9dcfbe","91b9cfa6-3b37-4a78-ae84-d56c026e3698",
                "a7dea07f-c24a-4b5c-a80f-f14222f4408d","19789ae0-eb4c-401b-81d5-6d13c38ca3c2",
                "9e294ad8-264d-43a2-b9d2-02aa3307bf19","b1fa9c0a-3f92-4df7-ba73-775e3845c77a",
                "3eb843d9-1a5a-47aa-a1f0-04dc6265f5b1","ef27cbfa-817d-42d1-a55b-f443334b5cb6",
                "8cebd1b6-f0d6-4a3f-b402-a172f6d5c9e8","ef1631f9-826d-40c9-beff-aaa0c79191cd",
                "916b4d0f-f9c3-4d19-b468-71ea6c738561","5f53e6f1-3e51-4f6b-a5af-df1e04f05e6f",
                "d12e2a8f-702c-462a-b3b6-4f3fdb9d0368","eb79c252-1843-4ddd-9cad-46119922e4c3",
                "3c1cdd30-f9cc-4068-9990-defeef13e961","ef015494-c6d6-408b-a1f3-9ab956e7354a",
                "122af8f0-12d3-4dcd-9758-a0d8a6f30713","84486d8d-183f-40cd-8652-2dbc63743fc5",
                "23dc3697-03da-4ede-a435-d6d4f33d6971",
                "27e4394b-ad21-44b9-aded-f539bc7d385b","3ad6dfa7-734d-41b0-ad66-15f0cff630d0",
                "faa87e86-4db0-43d2-9e84-9d9002cd5668","a847a390-cb2a-4ff8-8d52-847f23b0a8fa",
                "ed4ddff8-2740-438f-9af1-2e8125a96239","8285d0a3-7455-4aaf-9a59-deb9da16be48",
                "39ca8dd8-c721-43dc-ae3f-63ab895f73a1","4af2f63c-065a-4ddf-a232-fce87b21468c",
                "5211856d-cae1-40f8-aa8e-281fca0cc74d","443bb9dd-d449-4591-beb4-b8bfb9f22346",
                "8c170f00-5dea-4ff0-b65d-9545b5d2978c","df5f8115-ff60-40b0-8ba6-eaa61db9bbdf",
                "7bbda645-1e29-4e72-be97-8eb60962e06c","bce4b21f-e461-48ec-bca6-1d05776eacf3",
                "4abf867f-f25c-435c-8165-4417e4b32a5a","955c1e64-2817-4ab5-9719-b2a4d56fc40f",
                "965fb832-1e42-4aec-905e-95800c24a81a","97267ed3-479f-4c8a-ba46-cbeb86f393bd",
                "84283c7a-5e7f-40b3-9a84-2bab1c2f7149","0e2e1ebd-a9a9-4b83-8039-23bade3bffa3",
                "e37d3fb1-73f0-48f1-ae40-e8800c17fc99","72332335-d7ac-4e15-971c-dd08d6d58a79",
                "b03e9fb6-950f-4090-9f2e-1c45ca6d8afa","fde15027-bbcc-4877-bef7-a659f7bbd6ba",
                "a89e9db8-3b32-44a7-a89d-d73cde32ad5e","0b243f09-a843-46a1-9b56-a5e3ca1cedc7",
                "d57d71c1-73e7-410b-8ea7-25727df28bc3","4f2586ce-516b-45a5-a0ed-b8adc09c2bea",
                "b3948fcd-4ff8-431f-8509-c12a16b5fbff","75c37bf7-5a24-450f-9bb9-331533673cf9",
                "29b9063c-a697-4325-a727-a4a25f3f3773","9a0e80a7-e999-4c32-90f7-d9e603407f00",
                "6c87b8c6-9ccb-4853-8502-ef2970c67040","8bead04f-b325-4bfd-958b-66e58c57a6ca",
                "ac72ed95-381b-47b9-a64f-83feda6cb816","77a72b5e-6af3-4b43-8883-3fb420a2d3dc",
                "2e2b7a9c-e1ea-4e5a-a028-60d5cc6c347f","010cc77d-6f6c-4372-93e9-ccc910a6462d",
                "893975ee-3ed7-4d18-8e8b-ed3786975693","e40599a5-d4bb-4699-b999-30e63b2702aa",
                "9eb7d521-f3de-4907-b10b-13a0d3b78656","176bae20-55cf-44e4-8417-57c7a7e71be3",
                "f5bd72ae-919f-481c-8203-f878ee1c651a","6fd6a6b6-dea3-453a-90a6-92c2354628b6",
                "a150fea5-8e2d-440a-95dd-1401ff21dd59","a5d7a298-59f1-40f0-8123-05959fa8063a",
                "11fb7cce-06ad-445b-b09f-9dc1050e4cd4","e30ed1b7-15e9-4c2c-bfc5-c1a9c2a50a49",
                "7987cc31-fc59-43c0-b2c1-9897fc15847b","9ef5bf69-5fec-4430-8bda-8f00925f1c67",
                "5c579448-5b9d-428c-9105-e075a07446a9","9c1656ee-ab71-4002-ac32-75f922c40777",
                "b105f825-06d2-4f24-82aa-03bb7b4dde14","88712125-67e5-4e47-a413-e55b522a641c",
                "ebd2b092-16eb-450a-9167-3b8a81ad02d6","9dadf368-5c48-470b-bbd2-7cb381b69c47",
                "4cd0a36e-27ae-49ed-bad7-4172aa4517da",
                "9e8375c3-f911-495c-831d-78cc2405e22b","ff27b77a-2a64-4d88-8db0-c80bfb22c1ae",
                "3752825d-20e1-4f44-aa3e-93c5088017f2","1284fcaa-0ba7-4108-ac20-e52f9ad970b4",
                "a3b42528-3a10-4bc9-9616-8ad216be8ab6","2c2c38cc-80b5-4a14-bbc7-d6a158bc0be9",
                "b0476c2c-89de-4ac7-88f1-28bbea729bfa","674e3c95-e630-4946-beda-41408cb199b1",
                "e3efa84e-a672-43bc-a6ab-d054c7615747","bfc72f40-792e-41c4-9716-042e9f4ea984",
                "855dcdd1-f683-4bf1-b6de-662b29ad1d93","b8ce1e48-7f1d-4e4f-879b-137e40edd9fc",
                "058e0b03-f61c-41ea-9087-2a81b9925562","ecd78191-f7b5-492d-90e2-76223af1615b",
                "d22c6d78-d7f4-48ed-8afd-32bf7d586db8","9f2f8ca2-811a-47ed-b78b-ff06034699e0",
                "67239a38-d8cc-4c79-9e93-69b18b979671","dc21c55f-32c2-486e-89d7-70dbc96a3505",
                "0c19f3db-170e-4b26-8dea-61f3133e242b","63414088-75b3-4c72-9549-02af0fad12a1",
                "e1d686d9-7f3d-4e51-8bb4-30313a256826","daae4914-43a1-4bf4-8561-46d34c2dd86f",
                "1982fa58-6935-4f71-bd5c-f0841240e63d","239751d8-5a64-4814-980e-f07ed1b4f46e",
                "ce80ab5d-01b6-458b-8eb2-3aa4f9c01555","db415f39-a8a0-4ae9-8800-1790c6851d0a",
                "8bc15a42-cd6d-4ae5-88a3-6be8549cb696","4fb16f83-a6ce-4366-a49e-e216838910ac",
                "069a624e-eea3-40f9-b196-4d6ec48dcdca","9874825b-a625-4b9c-8ec8-378bedab0e2f",
                "51bb1c21-403d-4ae1-be90-03635b87b767","5830adb7-cb32-4f97-8ef2-8232d1ca08e1",
                "e9c70271-fb4e-4104-a301-9ab1c5f2c63b","8af79eb6-c201-4ac5-b3af-2ed47edf14a0",
                "9ce2d713-bc6b-441a-b476-c0d68ee9bc79","a16bbbbc-9e7b-470e-91ca-8a525721589e",
                "89a44d9a-ac72-4642-8c5e-a5eed7fb6717")



query2 <- query
tmp <- query$results[[1]]
tmp <- tmp[which(!duplicated(tmp$cases)),]
query2$results[[1]] <- tmp

query2$results[[1]] <- query2$results[[1]] %>% filter(!(id %in% black_list ))

#query2$results[[1]] <- query2$results[[1]][1:100,]

expdat=GDCprepare(query = query2,
                  save = TRUE,
                  save.filename = "~/workspace/GDCdata/TARGET-AML.rda")

clin_gbm <- GDCquery_clinic("TARGET-AML", "clinical")





cases_to_use <- query2$results[[1]]$cases
samples.df <- cases_to_use %>% as.data.frame %>%  tidyr::separate(".",1:5 %>% as.character,)



query2$results[[1]] <- query2$results[[1]] %>% filter(!(id %in% black_list ))

expdat=GDCprepare(query = query2,
                  save = TRUE,
                  save.filename = "~/workspace/GDCdata/TARGET-AML.rda")
