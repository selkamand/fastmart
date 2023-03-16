#' Default Cache
#'
#' @return default path to cache
#' @export
#'
#' @examples
#' fastmart_default_cache()
fastmart_default_cache <- function(){
  "~/.fastmart"
}

fastmart_init <- function(cache_dir = fastmart_default_cache()){
  if(file.exists(cache_dir))
    cli::cli_abort("Cache already exists at {.path {cache_dir}}. Doing nothing")
  else{
    cli::cli_alert("Creating Cache at: {.path {cache_dir}}")
    dir.create(cache_dir)
  }
}


get_sqlite_path <- function(cache_dir, GRCh){
  paste0(cache_dir, "/fastmart.GRCh",GRCh,".sqlite")
}


#' Cache Tables
#'
#' @param ensembl_version ensembl version. Older versions may not be supported. Will be ignored if GRCh != 38
#' @param cache_dir directory to store cache in
#' @param overwrite overwrite existing cache (flag)
#' @param GRCh assembly to use. Can only be set to
#' @return invisible(TRUE). This function is run for its side effects
#' @export
#'
#' @examples
#' \dontrun{
#' fastmart_cache_tables()
#' }
fastmart_cache_tables <- function(ensembl_version = 109, overwrite = FALSE, GRCh = c("38", "37"),
                                  cache_dir = fastmart_default_cache()
                                  ){
  # cli::cli_h2("Setup")
  GRCh <- rlang::arg_match(GRCh)
  cli::cli_progress_step("Creating MART object")

  assertions::assert_directory_exists(cache_dir, msg = "No existing cache found at {.path {cache_dir}}. You can run {.code fastmart_init({.path {cache_dir}})} to create a cache at this location.")
  path_fastmart_db <- get_sqlite_path(cache_dir, GRCh = GRCh)

  if(overwrite == FALSE & file.exists(path_fastmart_db)) cli::cli_abort("Found existing fastmart database at {.path {path_fastmart_db}}. \fTo overwrite existing database, supply {.arg overwrite = TRUE}")

  cli::cli_h1("Setup")
  # Fetch biomart (slow step)
  ensembl <- biomaRt::useEnsembl(biomart = 'genes',
                                 dataset = 'hsapiens_gene_ensembl',
                                 version = if(GRCh != "38") NULL else ensembl_version,
                                 GRCh = if(GRCh == "38") NULL else GRCh
                                 )

  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path_fastmart_db)
  cli::cli_progress_done()

  # Convert HGNC -> ENSEMPL gene ID -------------------------------------------------------------------------
  cli::cli_h1("HGNC -> ENSEMBL GENE")
  # Table name hgnc_to_ensembl
  cli::cli_progress_step("[HGNC -> ENSEMBL GENE] [RETRIEVING] retrieving map from ENSEMBL gene mart")
  df_gene_id_map <- biomaRt::getBM(
    attributes=c("ensembl_gene_id", "ensembl_gene_id_version", "hgnc_symbol", "gene_biotype", "chromosome_name", "start_position", "end_position"),
    filters = c("with_hgnc"),
    values = list(
      with_hgnc=TRUE
    ),
    mart = ensembl
  )

  df_gene_id_map[["hgnc_unique_id"]] <- paste(df_gene_id_map[["hgnc_symbol"]], df_gene_id_map[["chromosome_name"]], df_gene_id_map[["start_position"]],df_gene_id_map[["end_position"]])

  assertions::assert(
    !anyDuplicated(df_gene_id_map[["hgnc_unique_id"]]),
    msg = "ENSEMBL query returned ambiguous mappings between ENSEMBL gene IDs and hgnc_symbol,start_position,end_position"
    )

  cli::cli_progress_step("[HGNC -> ENSEMBL GENE] [CACHING] saving result to {.strong hgnc_to_ensembl} table in the fastmart database ({.path {path_fastmart_db}})")
  RSQLite::dbWriteTable(conn, "hgnc_to_ensembl", df_gene_id_map, overwrite = TRUE)

  cli::cli_progress_step("[HGNC -> ENSEMBL GENE] [INDEXING] indexing on symbol.chrom.start.end & ensembl_gene_id")
  RSQLite::dbExecute(conn, "CREATE INDEX hgnc_unique_id ON hgnc_to_ensembl(hgnc_unique_id)")
  RSQLite::dbExecute(conn, "CREATE INDEX hgnc_symbol_chr_start_end ON hgnc_to_ensembl(hgnc_symbol, chromosome_name, start_position, end_position)")
  RSQLite::dbExecute(conn, "CREATE INDEX ensembl_id ON hgnc_to_ensembl(ensembl_gene_id)", overwrite = TRUE)
  cli::cli_progress_done()
  cli::cli_alert_success("Successfully created HGNC -> Ensemble ID map at {.path {path_fastmart_db}}")
  cli::cli_alert_info("To actually convert your gene symbols, run {.code fastmart_hgnc_to_ensembl()}")



  # Annotate Gene Biotype gene ID -------------------------------------------------------------------------
  cli::cli_h1("ENSEMBL GENE <-> BIOTYPE")
  # table_name: ensembl_to_biotype
  cli::cli_progress_step("[ENSEMBL GENE <-> BIOTYPE] [RETRIEVING]")

  df_gene_biotype <- biomaRt::getBM(
    attributes=c(
      "ensembl_gene_id",
      "gene_biotype"
    ),
    mart = ensembl
  )

  assertions::assert_no_duplicates(df_gene_biotype[["ensembl_gene_id"]])
  cli::cli_progress_step("[ENSEMBL GENE <-> BIOTYPE] [CACHING] result to {.strong hgnc_to_ensembl} table in the fastmart database ({.path {path_fastmart_db}})")
  RSQLite::dbWriteTable(conn, "ensembl_to_biotype", df_gene_biotype, overwrite = TRUE)

  cli::cli_progress_step("[ENSEMBL GENE <-> BIOTYPE] [INDEXING]")
  RSQLite::dbExecute(conn, "CREATE INDEX ensembl_gene_id ON ensembl_to_biotype(ensembl_gene_id)")
  RSQLite::dbExecute(conn, "CREATE INDEX gene_biotype ON ensembl_to_biotype(gene_biotype)")

  cli::cli_progress_done()
  cli::cli_alert_success("Successfully created ENSEMBL GENE ID <-> BIOTYPE map at {.path {path_fastmart_db}}")
  cli::cli_alert_info("To annotate ENSEMBL genes with biotype, run {.code fastmart_convert_hgnc_to_ensembl()}")

  # Disconnect from database
  RSQLite::dbDisconnect(conn)


  return(invisible(TRUE))
}


#' Convert HGNC -> ENSEMBL IDs
#'
#'
#' Because ENSEMBL IDs are more granular than HGNC / HUGO symbols we require the start & stop coordinates
#'
#' @param hgnc_symbols HGNC symbols
#' @param chrom Chromosome of Gene. Required to resolve ambiguous mappings (character)
#' @param start Start position of Gene. Required to resolve ambiguous mappings  (numeric)
#' @param end End position of Gene. Required to resolve ambiguous mappings  (numeric)
#' @param exact_match  Is an exact match between chrom-start-end required (TRUE), or should the omst similar interval be returned (FALSE)  (flag)
#' @inheritParams fastmart_cache_tables
#'
#' @return matched ENSEMBL IDs (character)
#' @export
#'
#' @examples
#' \dontrun{
#' fastmart_convert_hgnc_to_ensembl('TBCE', chrom='1', start='235530675', end='235612283', GRCh = "37")
#' # > ENSG00000116957.8
#' }
#'
fastmart_convert_hgnc_to_ensembl <- function(hgnc_symbols, chrom, start, end, GRCh = c("38", "37"), exact_match = TRUE,cache_dir = fastmart_default_cache()){

  # Setup
  GRCh <- rlang::arg_match(GRCh)
  path_fastmart_db = get_sqlite_path(cache_dir, GRCh)
  assert_cach_dir_exists(cache_dir)
  assert_cache_db_exists(cache_dir, GRCh)

  # Connect to DB
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path_fastmart_db)

  # Preprocess query terms
  observed_unique_hgnc_ids <- paste(hgnc_symbols, chrom, start, end)
  unique_hgnc_id_set <- unique(observed_unique_hgnc_ids)
  sql_unique_hgnc_id_set <- format_as_sql(unique_hgnc_id_set)

  # Query DB

  if(exact_match){
    df_ids <- RSQLite::dbGetQuery(conn = conn, glue::glue_safe("
                                 SELECT
                                  *
                                 FROM
                                  hgnc_to_ensembl
                                 WHERE
                                 hgnc_unique_id IN {sql_unique_hgnc_id_set}
                                 "))

  # Return property of interest
    ensemble_ids <- df_ids[["ensembl_gene_id_version"]][match(observed_unique_hgnc_ids, df_ids[["hgnc_unique_id"]])]
  }
  else{
    key <- paste(hgnc_symbols, chrom, start, end)
    is_first_occurance <- !duplicated(key)
    df_input_uniq <- data.frame(
      hgnc_symbols=hgnc_symbols,
      chrom=chrom,
      start=suppressWarnings(as.numeric(start)),
      end=suppressWarnings(as.numeric(end)),
      key = key)[is_first_occurance,]
    df_input_uniq <- na.omit(df_input_uniq)

    hgnc_symbols_uniq <- unique(hgnc_symbols)
    sql_hgnc_symbols <- format_as_sql(unique(hgnc_symbols))

    df_ids <- RSQLite::dbGetQuery(conn = conn, glue::glue_safe("
                                 SELECT
                                  *
                                 FROM
                                  hgnc_to_ensembl
                                 WHERE
                                 hgnc_symbol IN {sql_hgnc_symbols}
                                 "))

    df_input_uniq[['ensembl_gene_id_version']] <- vapply(seq_len(nrow(df_input_uniq)), FUN = function(i){
      symbol = df_input_uniq[["hgnc_symbols"]][i]
      chr = df_input_uniq[["chrom"]][i]

      df_possible <- df_ids[df_ids[["hgnc_symbol"]] == symbol & df_ids[["chromosome_name"]] == chr, ]
      if(nrow(df_possible) == 0) return(NA_character_)

      closest_index <- most_similar_index(
        start1 = df_input_uniq[["start"]][i],
        end1 = df_input_uniq[["end"]][i],
        start2 = df_possible[["start_position"]],
        end2 = df_possible[["end_position"]]
        )

      df_possible[["ensembl_gene_id_version"]][closest_index]
      }, FUN.VALUE = character(1))

    ensemble_ids <- df_input_uniq[["ensembl_gene_id_version"]][match(key, df_input_uniq[["key"]])]
  }


  # Disconnect
  RSQLite::dbDisconnect(conn)
  return(ensemble_ids)

}

most_similar_index <- function(start1, end1, start2, end2){
  which.max(abs(start1-start2) + abs(end1-end2))
}


#' Annotate ENSEMBL GENES with BIOTYPE
#'
#' @param ensemble_gene_id ensembl gene id (can be versioned / not versioned)
#' @inheritParams fastmart_cache_tables
#'
#' @return matched biotype for each sample (character)
#' @export
#'
#' @examples
#' \dontrun{
#' fastmart_annotate_biotype('ENSG00000282455.1')
#' #> IG_D_gene
#' }
fastmart_annotate_biotype <- function(ensemble_gene_id, GRCh = c("38", "37"), cache_dir = fastmart_default_cache()){

  # Setup
  GRCh <- rlang::arg_match(GRCh)
  path_fastmart_db = get_sqlite_path(cache_dir, GRCh)
  assert_cach_dir_exists(cache_dir)
  assert_cache_db_exists(cache_dir, GRCh)

  # Connect to DB
  conn <- RSQLite::dbConnect(RSQLite::SQLite(), path_fastmart_db)

  # Preprocess query terms
  ensemble_gene_id <- fastmart_ensembl_drop_version(ensemble_gene_id)
  sql_ensemble_gene_id_set <- format_as_sql(unique(ensemble_gene_id))

  # Query DB
  df_result <- RSQLite::dbGetQuery(conn = conn, glue::glue_safe("
                               SELECT
                                *
                               FROM
                                ensembl_to_biotype
                               WHERE
                               ensembl_gene_id IN {sql_ensemble_gene_id_set}
                               "))

  # Disconnect
  RSQLite::dbDisconnect(conn)

  # Return property of interest
  df_result[["gene_biotype"]][match(ensemble_gene_id, df_result[["ensembl_gene_id"]])]
}

fastmart_ensembl_drop_version <-  function(ensembl_id_version){
  sub(x=ensembl_id_version, pattern = "\\.[0-9]+", replacement = "")
}

assert_cach_dir_exists <- assertions::assert_create(
  func = function(x){ file.exists(x) },
  default_error_msg = "Could not find cache directory: {.path {arg_value}}.\fPlease run fastmart_init('{.path {arg_value}}')  to initialise the cache directory or specify the path of an existing cach"
)

assert_cache_db_exists <- assertions::assert_create(
  func = function(x, GRCh){ file.exists(get_sqlite_path(x, GRCh)) },
  default_error_msg = "Could not find cache database: {.path {arg_value}}.\fPlease run {.code fastmart_cache_tables(cache_dir = '{.path {dirname(get_sqlite_path(x, GRCh))}}', GRCh = '{GRCh}')} to populate the cache before running "
)


format_as_sql <- function(x){
  paste0('(',paste0('"',x, '"', collapse = ", "), ')' )
}
