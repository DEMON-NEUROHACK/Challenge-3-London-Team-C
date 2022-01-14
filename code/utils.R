message_parallel <- function(...) {
    system(sprintf('echo "%s"', paste0(..., collapse = "")))
}

vcfs2ids <- function(vcfs){
    stringr::str_split(basename(vcfs),"[.]",simplify = TRUE)[,1]
}

genotype_matrix_grouped <- function(vcfs,
                                    grouping_var = NULL,
                                    structural_variants = TRUE,
                                    save_prefix = NULL,
                                    na_fill=NULL,
                                    mc.cores = 10){
  if(!is.null(grouping_var)){
    vcf_tmp <- VariantAnnotation::readVcf(vcfs[1]) 
    groups <- unique(vcf_tmp@info[[grouping_var]])
  }  
  #### Read in all VCFs first #####
  ids <- vcfs2ids(vcfs = vcfs)
  vcf_obj_list <- parallel::mclapply(vcfs, function(vcf_file){
    message_parallel("Reading: ",vcf_file)
    VariantAnnotation::readVcf(vcf_file)
  }) %>% `names<-`(ids)
  ##### Iterate over groups #####
  genotypes_list <- lapply(groups, function(group){ 
    message("Group: ",group)
    genotype_matrix(vcfs = vcf_obj_list, 
                    structural_variants = structural_variants,
                    grouping_var = grouping_var,
                    group = group,
                    save_prefix=save_prefix, 
                    na_fill=na_fill,
                    mc.cores = mc.cores)
  }) %>% `names<-`(groups)
  return(genotypes_list)
}


genotype_matrix <- function(vcfs,
                            structural_variants = FALSE,
                            grouping_var = NULL,
                            group = NULL,
                            save_prefix=NULL,#"data/SV_VCFs/by_variant/SV_VCFs",
                            na_fill=NULL,
                            mc.cores = 10){
   
  if(methods::is(vcfs[[1]],"VCF")){ 
    ids <- names(vcfs)
  }else {
    ids <- vcfs2ids(vcfs = vcfs)  
    names(vcfs) <- ids
  }
  #### Iterate over VCFs #####
  mat <- parallel::mclapply(names(vcfs), function(nm){
    vcf_file <- vcfs[[nm]]
    if(methods::is(vcf_file,"VCF")){
      vcf <- vcf_file
    } else {
      message_parallel("Reading: ",vcf_file)
      vcf <- VariantAnnotation::readVcf(vcf_file)
    }
    message_parallel("Parsing: ",nm)
    #### Extract encoded genotypes ####
    if(structural_variants){
       vec <- encode_structural_variants(vcf = vcf,
                                         grouping_var = grouping_var,
                                         group = group)
       r <- t(as.matrix(vec))
    } else {
      snpmat <- VariantAnnotation::genotypeToSnpMatrix(x = vcf)
      r <- as(snpmat$genotype, "numeric")
    }
    return(r)
  }, mc.cores = mc.cores) %>% 
    do.call(what="rbind") %>%
    `rownames<-`(ids) %>%
    methods::as("sparseMatrix")
  #### Fill NAs ####
  if(!is.null(na_fill)) mat[is.na(mat)] <- na_fill
  #### Save matrix ####
  save_matrix(mat = mat,
              save_prefix = save_prefix, 
              suffix = "genotype_matrix",
              group = group)
  #### Report ####
  message(paste(dim(mat),collapse = " x "))
  return(mat)
}

save_matrix <- function(mat,
                        save_prefix,
                        suffix = "matrix",
                        group = NULL){
  #### Save matrix ####
  if(!is.null(save_prefix)){
    try({
      if(any(group=="all")) group <- NULL
      save_path <- paste0(save_prefix,
                          if(!is.null(group)) paste0(".",group) else NULL,
                          paste0(".",suffix,".tsv.gz"))
      save_path <- gsub(":","-",save_path)
      dir.create(dirname(save_path), 
                 showWarnings = FALSE, recursive = TRUE)
      message_parallel("Saving ==> ",save_path)
      data.table::fwrite(
        x = data.table::data.table(as.matrix(mat),
                                   keep.rownames = "id"),
        file = save_path, 
        sep="\t") 
    })
    
  } 
}



add_alleles <- function(vcf,
                        dat){
  snpmat <- suppressMessages(suppressWarnings(VariantAnnotation::genotypeToSnpMatrix(vcf)))
  dat$A1 <- as.character(snpmat$map$allele.1) 
  a2 <- snpmat$map$allele.2
  names(a2) <-  snpmat$map$snp.names 
  a2 <- as.character( unlist(a2)) 
  dat$A2 <- a2[rownames(dat)]
  return(dat)
}



vcf2dt <- function(vcf, 
                   add_alleles = TRUE,
                   as_granges = FALSE){
  # gr <- methods::as(methods::as(vcf,"VRanges"),"GRanges")
  ranges <- vcf@rowRanges
  inf <- VariantAnnotation::info(vcf)
  gen <- VariantAnnotation::geno(vcf) 
  nms <- names(vcf)
  dat <- cbind(CHR=GenomicRanges::seqnames(ranges),
               START=GenomicRanges::start(ranges),
               END=GenomicRanges::end(ranges),
               name=nms,
               inf,
               GT=as.vector(gen$GT)
               ) %>% `rownames<-`(nms) 
  if(add_alleles){
    dat <- add_alleles(vcf = vcf,
                       dat = dat)
  }
  if(as_granges){ 
    gr <- GenomicRanges::makeGRangesFromDataFrame(df = data.frame(dat),
                                                  seqnames.field = "CHR",
                                                  start.field = "START",
                                                  end.field = "END", 
                                                  keep.extra.columns = TRUE)
    return(gr)
  } else {
    return(dat)
  }
}

encode_structural_variants <- function(vcf,
                                       key = get_key(),
                                       genotype_col = "GT",
                                       grouping_var = NULL,
                                       group = NULL){
  # gr <-  methods::as(vcf,"VRanges") 
  #### Converting the file to df this way is far more efficient than methods::as() ####
  dat <- vcf2dt(vcf = vcf)
  #### Subset to just one type #####
  if(!is.null(grouping_var)){
    dat <- dat[dat[[grouping_var]]==group,]
  } 
  #### Make genotype numeric #### 
  dat$GT_int <- key[dat[[genotype_col]]]
  vec <- dat$GT_int %>% `names<-`(rownames(dat)) 
  return(vec)
}



gather_vcfs <- function(vcfs,
                        cols=NULL,
                        mc.cores = 10){
    ids <- vcfs2ids(vcfs = vcfs)
    grl <- parallel::mclapply(vcfs, function(vcf_file){
        message_parallel("Reading: ",vcf_file)
        vcf <- VariantAnnotation::readVcf(file = vcf_file)
        gr <- vcf2dt(vcf = vcf,
                      add_alleles = TRUE,
                      as_granges = TRUE) 
        GenomicRanges::mcols(gr)$id <- vcfs2ids(vcf_file)  
        if(!is.null(cols)){
            cols <- cols[cols %in% colnames(GenomicRanges::mcols(gr))]
            return(gr[,cols])
        } else {
            return(gr)
        }
    }, mc.cores = mc.cores) %>%
      `names<-`(ids)
    return(grl)
}

merge_vcfs <- function(grl,
                       save_prefix=NULL #"data/SNP_VCFs/by_variant/SNP_VCFs"
                       ){ 
    gr_dt <- grl %>%
        GenomicRanges::GRangesList() %>%
        unlist() %>%
        data.frame() %>%
        data.table::data.table()
    #### Save matrix ####
    if(!is.null(save_prefix)){
        save_path <- paste0(save_prefix,".merged.tsv.gz")
        dir.create(dirname(save_path), 
                   showWarnings = FALSE, recursive = TRUE)
        message("Saving ==> ",save_path)
        data.table::fwrite(
            x = data.table::data.table(gr_dt,
                                       keep.rownames = "id"),
            file = save_path, 
            sep="\t")
    } 
    return(gr_dt)
}

get_key <- function(){
  # Follows conventions used by VariantAnnotation package:
  # https://rdrr.io/bioc/VariantAnnotation/man/genotypeToSnpMatrix-methods.html
  # 0 = missing,
  # 1 = "0/0"
  # 2 = "0/1" or "1/0" 
  # 3 = "1/1"
  ## Assuming "0"=="." in notations
    key <- c(" "=0,
             "NULL"=0,
             "NA"=0,
             "./."=1,
             "0/0"=1,
             "1/."=2,
             "1/0"=2,
             "./1"=2,
             "0/1"=2, 
             "1/1"=3,
             # Additional values (guessed)
             "1"=2,
             "1/2"=4,
             "2/1"=4,
             "2/2"=5)
    names(key)[names(key)=="NA"] <- NA
    return(key)
    # data.table::fwrite(
    #     data.frame(
    #         key = names(key),
    #         value = unname(key)
    #     ),
    #     "data/genotype_encodings.csv"
    # )
}


 

get_gene_counts <- function(gr_dt, 
                            gene_col = "ensembl_gene_id", # "CSQ"
                            extract_string = NULL, # "ENSG"
                            grouping_var = NULL, # "SVTYPE"
                            key = get_key()
                            ){
    if(gene_col!="ensembl_gene_id"){
        gr_dt <- gr_dt %>% 
          dplyr::rename(ensembl_gene_id =  dplyr::all_of(gene_col)) 
    }
    # extract_ensg <- function(val){
    #     # val = gr_dt_long$ensembl_gene_id[1]
    #     # val <- paste(val,collapse = "|")
    #     split <- strsplit(gsub("\\|+","|",val),"[|]")[[1]]
    #     list(split[startsWith(split,"ENSG")])
    # }  
    # test <- gr_dt_long[1:10,] %>% 
    #     dplyr::mutate(ensembl_gene_id=extract_ensg(ensembl_gene_id))
    
    # gsub('^.*ENSG.*|*$', '',gr_dt_long$ensembl_gene_id[1])    
    #### Unnest lists in gene col ####
    gr_dt_long <- tidyr::unnest_longer(gr_dt, col = ensembl_gene_id)
    #### Only extract ensembl genes ####
    if(!is.null(extract_string)){
      gr_dt_long$ensembl_gene_id <- 
        stringr::str_extract(gr_dt_long$ensembl_gene_id,
                             paste0(extract_string,"(.*?)[|]")) %>% 
        trimws(whitespace = "[|]")
    } 
    summarised_vars <- c("QUAL","GT_int","SVLEN")
    summarised_vars <- summarised_vars[summarised_vars %in% colnames(gr_dt_long)]
    gene_counts <- gr_dt_long %>% 
        dplyr::mutate(GT_int = key[GT]) %>%
        # subset(length(ensembl_gene_id[[1]])>0) %>%
        dplyr::group_by(.dots=c("id", grouping_var,"ensembl_gene_id")) %>%
        dplyr::summarise_at(.vars = c("GT_int",summarised_vars), 
                            .funs = function(x)mean(x, na.rm=TRUE) ) %>%
        data.table::data.table()
    return(gene_counts)
}


scl <- function(x){
    scales::rescale(x, to = c(0.000000001, 1))
}

add_gene_length <- function(gene_counts,
                            grouping_var = NULL  # "SVLEN"
                            ){ 
    gene_len <- gene_length(gene_list = unique(gene_counts$ensembl_gene_id))
    gene_counts2 <- gene_counts %>% 
        data.table::merge.data.table(gene_len, 
                                     by.x="ensembl_gene_id", 
                                     by.y = "ENSEMBL", 
                                     all.x = TRUE) %>%
        subset(!is.na(TXLENGTH)) %>%
        dplyr::mutate_at(.vars = c(grouping_var,"TXLENGTH"),
                         .funs = abs) %>% 
        dplyr::group_by(.dots = c("id", grouping_var, "ensembl_gene_id")) %>% 
        dplyr::summarise_at(.vars = c("GT_int",grouping_var,"QUAL","TXLENGTH"),
                            .funs = mean, na.rm = TRUE) %>%
        dplyr::group_by(.dots = c("id", grouping_var, "ensembl_gene_id")) %>% 
        dplyr::summarise(GT_int = as.integer(GT_int),
                         # SVLEN = as.integer(SVLEN),
                         QUAL = as.numeric(QUAL),
                         TXLENGTH = as.integer(TXLENGTH),
                         count = dplyr::n()) %>%
        dplyr::ungroup() %>% 
        dplyr::mutate_at(.vars = c("count","GT_int",grouping_var,
                                   "QUAL","TXLENGTH"), 
                         .funs = list(z=scl)) %>%
        data.table::data.table()
    return(gene_counts2)
}

add_gene_score <- function(gene_counts,
                           formula_str = "GT_int ~ count + SVLEN + TXLENGTH + QUAL" 
                           ){ 
  if(!"GT_int" %in% colnames(gene_counts)){
    key <- get_key()
    gene_counts$GT_int <- key[gene_counts$GT]
  }
    model <- stats::lm(data = gene_counts,
                       formula = stats::as.formula(formula_str),
                       na.action = na.exclude)
    gene_counts$score <- as.numeric(
        scales::rescale(
            log10(
                scl(
                    stats::residuals(model)
                    )
            )
        )
    )  
    gene_counts$score[is.na(gene_counts$score)] <- 0 
    hist(gene_counts$score, breaks = 100)
    return(gene_counts)
}




create_gene_matrix <- function(gene_counts,
                               grouping_var = NULL,# "SVTYPE"
                               data_type = c("gene_genotypes",
                                             "gene_scores"),
                               formula_str = "id__type ~ ensembl_gene_id",
                               remove_nonvarying = TRUE,
                               save_prefix = NULL){
  rn_col <- trimws(strsplit(formula_str,"~")[[1]][1])
  data_type <- tolower(data_type[1])
  message("Filling matrix with: ",data_type)
  if(is.null(grouping_var)){
    groups <- "all"
  } else {
    groups <- unique(gene_counts[[grouping_var]])
  }
  gene_counts <- gene_counts %>% 
    dplyr::mutate(id__type=trimws(paste(id,grouping_var,sep="__"),
                                  whitespace = "__")) #%>%
  #### Iterate over groups ####
  mat_list <- lapply(groups, function(group){ 
    if(group!="all"){
      message("Group: ",group)
      gene_counts <- gene_counts[get(grouping_var)==group,]
    }  
    M <- data.table::dcast.data.table(
      data = gene_counts, 
      # formula = id ~ SVTYPE,
      formula = stats::as.formula(formula_str),
      fun.aggregate = "mean",
      fill = 0,
      drop = FALSE,
      value.var = if(data_type=="gene_scores")"score"else{"GT_int"}, 
      na.rm=TRUE) %>%
      tibble::column_to_rownames(rn_col) %>%
      as.matrix() %>%
      as("sparseMatrix")
    #### Remove features that don't vary across participants ####
    if(remove_nonvarying){
      M <- M[,Matrix::colSums(M, na.rm = TRUE)>0, drop=FALSE]
    }
    message(group," : ",paste(dim(M),collapse = " x "))
    #### Save matrix ####
    save_matrix(mat = M,
                save_prefix = save_prefix, 
                suffix = paste(data_type,"matrix",sep="_"),
                group = group)
    return(M)
  }) %>% `names<-`(groups) 
  #### Return ####
  if(length(mat_list)==1 && names(mat_list)[1]=="all"){
    return(mat_list[[1]])
  } else {
    return(mat_list)
  }
}


create_corr_heatmaps <- function(mat_list,
                                 data_type,
                                 remove_prefix = "LP6005878-DNA_",
                                 save_prefix=NULL#"figures/heatmaps/SV_"
                                 ){
  cor_mat_list <- lapply(names(mat_list), function(x){
    #### corr ####
    mat <- mat_list[[x]]
    cor_mat <- cor(as.matrix(Matrix::t(mat)))
    cor_mat[is.na(cor_mat)] <- 0
    diag(cor_mat) <- NA
    #### Remove prefixes in names ####
    cor_mat2 <- cor_mat
    colnames(cor_mat2) <- gsub(remove_prefix,"",colnames(cor_mat2))
    rownames(cor_mat2) <- gsub(remove_prefix,"",rownames(cor_mat2))
    #### plot ####
    if(!is.null(save_prefix)){
      path <- paste0(save_prefix,
                     gsub(":","-",x),".",data_type,".pdf")
      grDevices::pdf(path)
    } 
    stats::heatmap(cor_mat2, symm = TRUE)
    if(!is.null(save_prefix)){
      grDevices::dev.off()  
    }
    return(cor_mat)
  }) %>% `names<-`(names(mat_list))
  return(cor_mat_list)
}


createDT <- function (DF, caption = "", scrollY = 400) 
{
  data <- DT::datatable(DF, caption = caption, extensions = "Buttons", 
                        options = list(dom = "Bfrtip", buttons = c("copy", "csv", 
                                                                   "excel", "pdf", "print"), scrollY = scrollY, scrollX = T, 
                                       scrollCollapse = T, paging = F, columnDefs = list(list(className = "dt-center", 
                                                                                              targets = "_all"))))
  return(data)
}

drop_list_cols <- function(gr_dt){
  classes <- sapply(gr_dt,class)
  gr_dt[,c(names(..classes)[..classes!="AsIs"])]
}
