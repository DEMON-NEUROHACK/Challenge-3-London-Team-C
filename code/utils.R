message_parallel <- function(...) {
    system(sprintf('echo "%s"', paste0(..., collapse = "")))
}

vcfs2ids <- function(vcfs){
    stringr::str_split(basename(vcfs),"[.]",simplify = TRUE)[,1]
}


gather_vcfs <- function(vcfs,
                        cols=NULL,
                        mc.cores = 10){
    grl <- parallel::mclapply(vcfs, function(vcf){
        message_parallel(vcf)
        out <- VariantAnnotation::readVcf(file = vcf)
        out <- methods::as(methods::as(out,"VRanges"),"GRanges")
        id <-  strsplit(basename(vcf),"[.]")[[1]][1]
        # names(out) <- paste(id,names(out),sep=":")
        GenomicRanges::mcols(out)$id <- id
        # out_dt <- data.table::data.table(data.frame(out))
        # out_dt <- out_dt[,.(seqnames, start, end, width, strand, 
        #                     QUAL, SVTYPE, ensembl_gene_id)]
        if(!is.null(cols)){
            cols <- cols[cols %in% colnames(GenomicRanges::mcols(gr))]
            return(out[,cols])
        } else {
            return(out)
        }
    }, mc.cores = mc.cores)
    return(grl)
}

merge_vcfs <- function(grl){ 
    gr_dt <- grl %>%
        GenomicRanges::GRangesList() %>%
        unlist() %>%
        data.frame() %>%
        data.table::data.table()
    return(gr_dt)
}

get_key <- function(){
    key <- c("./."=1,
      "0/0"=1,
      "1/."=0,
      "1/0"=0,
      "0/1"=0,
      "1"=0,
      "1/1"=2,
      "1/2"=3,
      "2/2"=4)
    return(key)
    # data.table::fwrite(
    #     data.frame(
    #         key = names(key),
    #         value = unname(key)
    #     ),
    #     "data/genotype_encodings.csv"
    # )
}

genotype_matrices <- function(grl,
                              grouping_var = NULL,
                              genotype_col = "GT",
                              save_prefix="data/SV_VCFs/by_variant/SV_VCFs",
                              key = get_key()
                              ){ 
    if(is.null(grouping_var)){
        types <- "all"
    } else {
        types <- unique(GenomicRanges::mcols(grl[[1]])[[grouping_var]])
    }
    genotypes_list <- lapply(types, function(x){
        mat <- lapply(grl,function(gr){
            # gr <- grl[[1]]
            #### Subset to just one type #####
            if(is.null(grouping_var)){
                GenomicRanges::mcols(gr)$all <- "all"
            } else {
                gr <- gr[GenomicRanges::mcols(gr)[[grouping_var]]==x,]
            } 
            #### Make genotype numeric #### 
            GenomicRanges::mcols(gr)$GT_int <- 
                key[GenomicRanges::mcols(gr)[[genotype_col]]]
            vec <- gr$GT_int %>% `names<-`(names(gr)) 
            return(vec)
        }) %>%
            `names<-`(all_ids) %>% 
            do.call(what = "cbind") %>%
            as("sparseMatrix") %>%
            Matrix::t()
        #### Save matrix ####
        save_path <- paste0(save_prefix,".",gsub(":","-",x),".tsv.gz")
        dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
        message("Saving ==> ",save_path)
        data.table::fwrite(
            x = data.table::data.table(as.matrix(mat), keep.rownames = "id"),
            file = save_path, 
            sep="\t")
        message(x," : ",paste(dim(mat),collapse = " x "))
        return(mat)
    }) %>% `names<-`(types)
    return(genotypes_list)
}

get_gene_counts <- function(gr_dt,
                            key = get_key(),
                            gene_col = "ensembl_gene_id" # "CSQ"
                            ){
    if(gene_col!="ensembl_gene_id"){
        gr_dt <- gr_dt %>% dplyr::rename(ensembl_gene_id =  gene_col)
        # gr_dt$ensembl_gene_id <- gr_dt[[gene_col]]
    }
    
    ensg <- lapply(gr_dt_long$ensembl_gene_id)
    test <- gr_dt_long[1:10,] %>% 
        dplyr::mutate(ensembl_gene_id=extract_ensg(ensembl_gene_id))
    
    gr_dt_long <- tidyr::unnest_longer(gr_dt, col = ensembl_gene_id)
    
    extract_ensg <- function(val){
        # val = gr_dt_long$ensembl_gene_id[1]
        # val <- paste(val,collapse = "|")
        split <- strsplit(gsub("\\|+","|",val),"[|]")[[1]]
        list(split[startsWith(split,"ENSG")])
    }  
    
    gene_counts <- gr_dt_long %>% 
        dplyr::mutate(GT_int = key[GT]) %>%
        # subset(length(ensembl_gene_id[[1]])>0) %>%
        dplyr::group_by(id, SVTYPE, ensembl_gene_id) %>%
        dplyr::summarise(QUAL = mean(QUAL,na.rm=TRUE), 
                         GT_int = mean(GT_int, na.rm=TRUE), 
                         SVLEN = mean(SVLEN, na.rm=TRUE)) %>%
        data.table::data.table()
    return(gene_counts)
}


scl <- function(x){
    scales::rescale(x, to = c(0.000000001, 1))
}

add_gene_length <- function(gene_counts){
    gene_len <- gene_length(gene_list = unique(gene_counts$ensembl_gene_id))
    sv_gene_counts <- gene_counts %>% 
        data.table::merge.data.table(gene_len, 
                                     by.x="ensembl_gene_id", 
                                     by.y = "ENSEMBL", 
                                     all.x = TRUE) %>%
        subset(!is.na(TXLENGTH)) %>%
        dplyr::mutate_at(.vars = c("SVLEN","TXLENGTH"),
                         .funs = abs) %>% 
        dplyr::group_by(id, SVTYPE, ensembl_gene_id) %>% 
        dplyr::summarise_at(.vars = c("GT_int","SVLEN","QUAL","TXLENGTH"),
                            .funs = mean, na.rm = TRUE) %>%
        dplyr::group_by(id, SVTYPE, ensembl_gene_id) %>% 
        dplyr::summarise(GT_int = as.integer(GT_int),
                         SVLEN = as.integer(SVLEN),
                         QUAL = as.numeric(QUAL),
                         TXLENGTH = as.integer(TXLENGTH),
                         count = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(id__type=paste(id,SVTYPE,sep="__")) %>%
        dplyr::mutate_at(.vars = c("count","GT_int","SVLEN","QUAL","TXLENGTH"), 
                         .funs = list(z=scl)) %>%
        data.table::data.table()
    return(sv_gene_counts)
}

add_gene_score <- function(sv_gene_counts){ 
    model <- stats::lm(data = sv_gene_counts,
                       formula = GT_int ~ count + SVLEN + TXLENGTH + QUAL,
                       na.action = na.exclude)
    sv_gene_counts$score <- as.numeric(
        scales::rescale(
            log10(
                scl(
                    stats::residuals(model)
                    )
            )
        )
    )  
    sv_gene_counts$score[is.na(sv_gene_counts$score)] <- 0 
    hist(sv_gene_counts$score, breaks = 100)
    return(sv_gene_counts)
}
