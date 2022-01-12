gene_length <- function(gene_list, 
                        keytype = "ENSEMBL",
                        columns = c("ENSEMBL","SYMBOL",
                                    "TXNAME","TXSTART","TXEND")){ 
    db <- Homo.sapiens::Homo.sapiens
    # AnnotationDbi::columns(db)
    len <- AnnotationDbi::select(db,
                                 keys = gene_list,
                                 keytype = "ENSEMBL", 
                                 columns = columns)
    mean_len <- len %>% 
        dplyr::mutate(TXLENGTH=abs(TXSTART-TXEND)) %>%
        dplyr::group_by(.dots= c(keytype,"SYMBOL") ) %>%
        dplyr::summarise(TXLENGTH=mean(TXLENGTH,na.rm=T)) %>%
        data.table::data.table()
    return(mean_len)
}