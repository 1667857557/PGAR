suppressPackageStartupMessages({
  library(Rcpp)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(S4Vectors)
  library(BSgenome)
  library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
})
snp_dat <- MungeSumstats::load_snp_loc_data("GRCh37", dbSNP=155)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)

Rcpp::cppFunction('
#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <string>
#include <algorithm>
using namespace Rcpp;

static inline int acgt_code(char c){
  c = (char)std::toupper((unsigned char)c);
  if(c==\'A\') return 1;
  if(c==\'C\') return 2;
  if(c==\'G\') return 3;
  if(c==\'T\') return 4;
  return 0;
}
static inline char code_char(int code){
  if(code==1) return \'A\';
  if(code==2) return \'C\';
  if(code==3) return \'G\';
  if(code==4) return \'T\';
  return \'N\';
}
static inline int comp_code(int code){
  if(code==1) return 4;
  if(code==4) return 1;
  if(code==2) return 3;
  if(code==3) return 2;
  return 0;
}
static inline bool is_delim(char c){
  return (c==\'/\' || c==\'|\' || c==\',\' || c==\';\' || c==\' \' || c==\'\\t\' || c==\'\\n\' || c==\'\\r\');
}
static inline bool is_pal_pair(int a, int b){
  int smin = std::min(a,b), smax = std::max(a,b);
  return (smin==1 && smax==4) || (smin==2 && smax==3);
}
static inline bool is_digits_only(const std::string &s){
  if(s.empty()) return false;
  for(size_t i=0;i<s.size();++i){
    if(s[i] < \'0\' || s[i] > \'9\') return false;
  }
  return true;
}
static inline std::string normalize_rsid(const std::string &s0){
  if(s0.empty()) return std::string();
  std::string s = s0;
  while(!s.empty() && (s.back()==\' \' || s.back()==\'\\t\' || s.back()==\'\\n\' || s.back()==\'\\r\')) s.pop_back();
  size_t p=0;
  while(p<s.size() && (s[p]==\' \' || s[p]==\'\\t\' || s[p]==\'\\n\' || s[p]==\'\\r\')) p++;
  if(p>0) s = s.substr(p);
  if(s.empty()) return std::string();
  if(is_digits_only(s)) return std::string("rs") + s;
  if(s.size()>=2){
    char r0 = (char)std::tolower((unsigned char)s[0]);
    char r1 = (char)std::tolower((unsigned char)s[1]);
    if(r0==\'r\' && r1==\'s\') return std::string("rs") + s.substr(2);
  }
  return s;
}

static inline void parse_two_alleles(const std::string &s, int &a1, int &a2, bool &multiallelic, bool &indel){
  a1 = 0; a2 = 0;
  multiallelic = false;
  indel = false;
  std::string tok;
  auto flush_tok = [&](){
    if(tok.size()==0) return;
    if(tok.size()!=1){
      indel = true;
    }else{
      int c = acgt_code(tok[0]);
      if(c){
        if(a1==0) a1 = c;
        else if(a2==0 && c!=a1) a2 = c;
        else if(c!=a1 && c!=a2) multiallelic = true;
      }
    }
    tok.clear();
  };
  for(size_t i=0;i<s.size();++i){
    char ch = s[i];
    if(is_delim(ch)) flush_tok();
    else tok.push_back(ch);
  }
  flush_tok();
  if(a1==0 || a2==0){
    if(a1==0 && a2!=0) std::swap(a1,a2);
  }
}

static inline void best_update(int q,
                               int prio,
                               int op,
                               const std::string &status,
                               const std::string &rsid,
                               int ra1, int ra2,
                               IntegerVector &best_prio,
                               IntegerVector &best_op,
                               CharacterVector &best_status,
                               CharacterVector &best_rsid,
                               CharacterVector &best_ra1,
                               CharacterVector &best_ra2){
  if(prio < best_prio[q]){
    best_prio[q] = prio;
    best_op[q] = op;
    best_status[q] = status;
    best_rsid[q] = rsid;
    best_ra1[q] = std::string(1, code_char(ra1));
    best_ra2[q] = std::string(1, code_char(ra2));
  }
}

// [[Rcpp::export]]
DataFrame cpp_best_harmonize_pos(IntegerVector qh,
                                 CharacterVector db_ref,
                                 CharacterVector db_alt,
                                 CharacterVector db_alleles_raw,
                                 NumericVector dbmaf,
                                 CharacterVector rsid_raw,
                                 CharacterVector qa1,
                                 CharacterVector qa2,
                                 NumericVector qeaf,
                                 int n_query,
                                 bool use_eaf,
                                 double maf_tol,
                                 double eaf_tol){

  IntegerVector best_prio(n_query);
  IntegerVector best_op(n_query);
  LogicalVector hit(n_query);
  CharacterVector best_status(n_query);
  CharacterVector best_rsid_match(n_query);
  CharacterVector best_ra1(n_query);
  CharacterVector best_ra2(n_query);
  CharacterVector rsid_pos(n_query);

  IntegerVector qa1c(n_query);
  IntegerVector qa2c(n_query);
  LogicalVector q_ok(n_query);

  for(int q=0;q<n_query;++q){
    best_prio[q] = 999999;
    best_op[q] = -1;
    hit[q] = false;
    best_status[q] = NA_STRING;
    best_rsid_match[q] = NA_STRING;
    best_ra1[q] = NA_STRING;
    best_ra2[q] = NA_STRING;
    rsid_pos[q] = NA_STRING;

    std::string s1 = (qa1[q] == NA_STRING) ? std::string() : as<std::string>(qa1[q]);
    std::string s2 = (qa2[q] == NA_STRING) ? std::string() : as<std::string>(qa2[q]);
    int c1 = (s1.size()==1) ? acgt_code(s1[0]) : 0;
    int c2 = (s2.size()==1) ? acgt_code(s2[0]) : 0;
    qa1c[q] = c1;
    qa2c[q] = c2;
    q_ok[q] = (c1!=0 && c2!=0 && s1.size()==1 && s2.size()==1);
  }

  int n = qh.size();
  for(int i=0;i<n;++i){
    int q = qh[i] - 1;
    if(q < 0 || q >= n_query) continue;
    hit[q] = true;

    std::string rs = (rsid_raw[i] == NA_STRING) ? std::string() : as<std::string>(rsid_raw[i]);
    rs = normalize_rsid(rs);
    if(rsid_pos[q] == NA_STRING && !rs.empty()){
      rsid_pos[q] = rs;
    }

    if(!q_ok[q]){
      if(best_prio[q] > 800){
        best_update(q, 800, -1, "no_or_invalid_query_allele", rs, 0, 0, best_prio, best_op, best_status, best_rsid_match, best_ra1, best_ra2);
      }
      continue;
    }

    int ra1 = 0, ra2 = 0;
    bool multiallelic = false, indel = false;

    if(db_ref[i] != NA_STRING && db_alt[i] != NA_STRING){
      std::string r = as<std::string>(db_ref[i]);
      std::string a = as<std::string>(db_alt[i]);
      if(r.size()==1 && a.size()==1){
        ra1 = acgt_code(r[0]);
        ra2 = acgt_code(a[0]);
      }else{
        indel = true;
      }
    }else{
      std::string dbs = (db_alleles_raw[i] == NA_STRING) ? std::string() : as<std::string>(db_alleles_raw[i]);
      parse_two_alleles(dbs, ra1, ra2, multiallelic, indel);
    }

    if(indel || ra1==0 || ra2==0){
      best_update(q, 650, -1, "indel_or_missing_ref", rs, 0, 0, best_prio, best_op, best_status, best_rsid_match, best_ra1, best_ra2);
      continue;
    }
    if(multiallelic){
      best_update(q, 600, -1, "multiallelic_ref", rs, ra1, ra2, best_prio, best_op, best_status, best_rsid_match, best_ra1, best_ra2);
      continue;
    }

    bool pal = is_pal_pair(ra1, ra2);

    int op = -1;
    int prio = 500;
    std::string status = "mismatch";

    int q1 = qa1c[q];
    int q2 = qa2c[q];

    if(q1==ra1 && q2==ra2){
      op = 0;
      prio = pal ? 120 : 100;
      status = pal ? "aligned_palindromic" : "aligned";
    }else if(q1==ra2 && q2==ra1){
      op = 1;
      prio = pal ? 130 : 110;
      status = pal ? "swapped_palindromic" : "swapped";
    }else{
      int c1 = comp_code(q1);
      int c2 = comp_code(q2);
      if(c1==ra1 && c2==ra2){
        op = 2;
        prio = pal ? 220 : 200;
        status = pal ? "strand_flipped_palindromic_ambiguous" : "strand_flipped";
      }else if(c1==ra2 && c2==ra1){
        op = 3;
        prio = pal ? 230 : 210;
        status = pal ? "strand_flipped_swapped_palindromic_ambiguous" : "strand_flipped_swapped";
      }else{
        op = -1;
        prio = 900;
        status = "mismatch";
      }
    }

    if(pal && use_eaf){
      double qf = std::numeric_limits<double>::quiet_NaN();
      if(!NumericVector::is_na(qeaf[q])) qf = qeaf[q];
      double dbf = std::numeric_limits<double>::quiet_NaN();
      if(!NumericVector::is_na(dbmaf[i])) dbf = dbmaf[i];

      if(!std::isnan(qf)){
        if(std::fabs(qf - 0.5) <= eaf_tol){
          status = "palindromic_ambiguous_eaf_near_0.5";
          prio += 50;
        }else if(!std::isnan(dbf)){
          double d1 = std::fabs(qf - dbf);
          double d2 = std::fabs(qf - (1.0 - dbf));
          if(std::min(d1,d2) > maf_tol){
            status = "palindromic_ambiguous_freq_disagree";
            prio += 40;
          }
        }
      }
    }else if(pal && !use_eaf){
      if(op==2 || op==3){
        status = "palindromic_ambiguous_strand";
        prio += 60;
      }
    }

    best_update(q, prio, op, status, rs, ra1, ra2, best_prio, best_op, best_status, best_rsid_match, best_ra1, best_ra2);
  }

  IntegerVector qi(n_query);
  for(int q=0;q<n_query;++q) qi[q] = q + 1;

  return DataFrame::create(
    Named("qi") = qi,
    Named("hit") = hit,
    Named("rsid_pos") = rsid_pos,
    Named("rsid_match") = best_rsid_match,
    Named("ref_a1") = best_ra1,
    Named("ref_a2") = best_ra2,
    Named("op") = best_op,
    Named("status") = best_status,
    Named("prio") = best_prio
  );
}
')

.pick_first_col <- function(nms, preferred){
  hit <- intersect(preferred, nms)
  if(length(hit)) hit[1] else NA_character_
}

.extract_db_fields_fast <- function(snp_hits){
  md <- GenomicRanges::mcols(snp_hits)
  nms <- names(md)
  
  ref_col <- if(all(c("REF","ALT") %in% nms)) "REF" else NA_character_
  alt_col <- if(all(c("REF","ALT") %in% nms)) "ALT" else NA_character_
  
  allele_preferred <- c("alleles_as_ambig","alleles_as_string","alleles","Alleles","allele","ALLELES","ALLELE","REF_ALT","ref_alt","ref/alt")
  freq_preferred <- c("MAF","maf","AF","af","EAF","eaf","Freq","freq","aaf","AAF")
  rs_preferred <- c("RefSNP_id","RefSNP","rsid","RSID","rsID")
  
  allele_col <- .pick_first_col(nms, allele_preferred)
  freq_col <- .pick_first_col(nms, freq_preferred)
  rs_col <- .pick_first_col(nms, rs_preferred)
  
  n <- length(snp_hits)
  
  db_ref <- if(!is.na(ref_col)) as.character(md[[ref_col]]) else rep(NA_character_, n)
  db_alt <- if(!is.na(alt_col)) as.character(md[[alt_col]]) else rep(NA_character_, n)
  dballeles <- if(!is.na(allele_col)) as.character(md[[allele_col]]) else rep(NA_character_, n)
  dbmaf <- if(!is.na(freq_col)) as.numeric(md[[freq_col]]) else rep(NA_real_, n)
  rsid <- if(!is.na(rs_col)) as.character(md[[rs_col]]) else rep(NA_character_, n)
  
  list(db_ref=db_ref, db_alt=db_alt, dballeles=dballeles, dbmaf=dbmaf, rsid=rsid)
}

harmonize_gwas_with_snplocs <- function(A,
                                        chr_col="CHR", pos_col="POS",
                                        ea_col="EA", oa_col="NEA",
                                        beta_col="BETA", eaf_col="EAF",
                                        snp_dat=NULL,
                                        use_eaf = !is.null(A[[eaf_col]]),
                                        chunk_size=100000L,
                                        maf_tol=0.10, eaf_tol=0.10,
                                        show_progress=TRUE){
  
  if(is.null(snp_dat)) snp_dat <- SNPlocs.Hsapiens.dbSNP155.GRCh37
  
  n <- nrow(A)
  chr <- as.character(A[[chr_col]])
  pos <- as.integer(A[[pos_col]])
  ea <- toupper(as.character(A[[ea_col]]))
  oa <- toupper(as.character(A[[oa_col]]))
  beta <- as.numeric(A[[beta_col]])
  eaf <- if(!is.null(A[[eaf_col]])) as.numeric(A[[eaf_col]]) else rep(NA_real_, n)
  
  out_rsid <- rep(NA_character_, n)
  out_ref_a1 <- rep(NA_character_, n)
  out_ref_a2 <- rep(NA_character_, n)
  out_status <- rep("no_hit", n)
  out_op <- rep(NA_integer_, n)
  
  style_target <- GenomeInfoDb::seqlevelsStyle(snp_dat)[1]
  
  idx_starts <- seq.int(1L, n, by=as.integer(chunk_size))
  n_chunks <- length(idx_starts)
  pb <- NULL
  if(show_progress) pb <- txtProgressBar(min=0, max=n_chunks, style=3)
  
  for(k in seq_along(idx_starts)){
    i1 <- idx_starts[k]
    i2 <- min(n, i1 + as.integer(chunk_size) - 1L)
    m <- i2 - i1 + 1L
    
    gr <- GenomicRanges::GRanges(seqnames=chr[i1:i2], ranges=IRanges::IRanges(start=pos[i1:i2], end=pos[i1:i2]))
    GenomeInfoDb::seqlevelsStyle(gr) <- style_target
    
    hits <- snpsByOverlaps(snp_dat, ranges=gr)
    if(length(hits) > 0){
      ol <- GenomicRanges::findOverlaps(gr, hits)
      qh <- S4Vectors::queryHits(ol)
      sh <- S4Vectors::subjectHits(ol)
      
      fields <- .extract_db_fields_fast(hits)
      db_ref <- fields$db_ref[sh]
      db_alt <- fields$db_alt[sh]
      dballeles <- fields$dballeles[sh]
      dbmaf <- fields$dbmaf[sh]
      rsid_ol <- fields$rsid[sh]
      
      best <- cpp_best_harmonize_pos(
        qh,
        db_ref,
        db_alt,
        dballeles,
        dbmaf,
        rsid_ol,
        ea[i1:i2],
        oa[i1:i2],
        eaf[i1:i2],
        m,
        isTRUE(use_eaf),
        maf_tol,
        eaf_tol
      )
      
      hit <- best$hit
      if(any(hit)){
        loc_all <- which(hit)
        gid_all <- (i1 - 1L) + loc_all
        
        out_status[gid_all] <- as.character(best$status[loc_all])
        out_op[gid_all] <- as.integer(best$op[loc_all])
        
        rs_pos <- as.character(best$rsid_pos)
        rs_match <- as.character(best$rsid_match)
        rs_pos[rs_pos == ""] <- NA_character_
        rs_match[rs_match == ""] <- NA_character_
        
        rs_fill <- rs_pos
        ok_use_match <- hit & !is.na(best$op) & (best$op >= 0L) & !is.na(rs_match)
        rs_fill[ok_use_match] <- rs_match[ok_use_match]
        need_match_if_pos_missing <- hit & (is.na(rs_fill) | !nzchar(rs_fill)) & !is.na(rs_match)
        rs_fill[need_match_if_pos_missing] <- rs_match[need_match_if_pos_missing]
        
        out_rsid[gid_all] <- rs_fill[loc_all]
        
        ok_ref <- hit & !is.na(best$op) & (best$op >= 0L)
        if(any(ok_ref)){
          loc_ok <- which(ok_ref)
          gid_ok <- (i1 - 1L) + loc_ok
          out_ref_a1[gid_ok] <- as.character(best$ref_a1[loc_ok])
          out_ref_a2[gid_ok] <- as.character(best$ref_a2[loc_ok])
        }
      }
    }
    
    if(show_progress) setTxtProgressBar(pb, k)
  }
  
  if(show_progress) close(pb)
  
  flip_comp <- function(x){
    x <- toupper(x)
    x[x=="A"] <- "t"
    x[x=="T"] <- "a"
    x[x=="C"] <- "g"
    x[x=="G"] <- "c"
    toupper(x)
  }
  
  idx_swap <- which(out_op == 1L)
  idx_sflip <- which(out_op == 2L)
  idx_sflip_swap <- which(out_op == 3L)
  
  if(length(idx_sflip)){
    ea[idx_sflip] <- flip_comp(ea[idx_sflip])
    oa[idx_sflip] <- flip_comp(oa[idx_sflip])
  }
  if(length(idx_sflip_swap)){
    ea[idx_sflip_swap] <- flip_comp(ea[idx_sflip_swap])
    oa[idx_sflip_swap] <- flip_comp(oa[idx_sflip_swap])
  }
  if(length(idx_swap)){
    tmp <- ea[idx_swap]; ea[idx_swap] <- oa[idx_swap]; oa[idx_swap] <- tmp
    beta[idx_swap] <- -beta[idx_swap]
    if(!is.null(A[[eaf_col]])) eaf[idx_swap] <- 1 - eaf[idx_swap]
  }
  if(length(idx_sflip_swap)){
    tmp <- ea[idx_sflip_swap]; ea[idx_sflip_swap] <- oa[idx_sflip_swap]; oa[idx_sflip_swap] <- tmp
    beta[idx_sflip_swap] <- -beta[idx_sflip_swap]
    if(!is.null(A[[eaf_col]])) eaf[idx_sflip_swap] <- 1 - eaf[idx_sflip_swap]
  }
  
  A$rsid <- out_rsid
  A$ref_a1 <- out_ref_a1
  A$ref_a2 <- out_ref_a2
  A$harmonize_status <- out_status
  A$harmonize_op <- out_op
  A[[ea_col]] <- ea
  A[[oa_col]] <- oa
  A[[beta_col]] <- beta
  if(!is.null(A[[eaf_col]])) A[[eaf_col]] <- eaf
  A
}
pacman::p_load("vroom","data.table","readr","tidyr","dplyr","devtools","ggplot2","tidyverse","GWAS.utils")
setwd("G:/HCC/Marginal_byChr/")
txt_files <- list.files(pattern = "\\.tsv$", full.names = TRUE)
rm(combined_data)
gc()
combined_data <- rbindlist(lapply(txt_files, vroom))

combined_data<-select(combined_data,c("GeneSymbol","ENSG","CHR","POS","EA","NEA","Beta","SE","PVAL","EAF","N"))
gc()

combined_data <- harmonize_gwas_with_snplocs(combined_data,
  chr_col="CHR",
  pos_col="POS",
  ea_col="EA",
  oa_col="NEA",
  beta_col="Beta",
  eaf_col="EAF",
  snp_dat=SNPlocs.Hsapiens.dbSNP155.GRCh37,
  use_eaf=TRUE,
  chunk_size=100000L,
  maf_tol=0.10,
  eaf_tol=0.10,
  show_progress=TRUE)

combined_data<-select(combined_data,c("GeneSymbol","ENSG","rsid","CHR","POS","EA","NEA","Beta","SE","PVAL","EAF","N"))
colnames(combined_data)<-c("Symbol","ENSGID","SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF","N")
combined_data<-na.omit(combined_data)
gc()
write_tsv(combined_data,file = "liver_eqtl.txt.gz")
combined_data<-select(combined_data,c("ENSGID","SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF","N"))
combined_data<-vroom("liver_eqtl.txt.gz")
gc()
# 假设你的数据框名为combined_data，ENSGID列包含基因ID
# 按ENSGID分割数据框
split_data <- split(combined_data, combined_data$ENSGID)

output_dir <- "split_files"
dir.create(output_dir, showWarnings = FALSE)
library(future.apply)
plan(multisession, workers = max(1, parallel::detectCores() - 1))

future_lapply(names(split_data), function(ensgid) {
  file_name <- file.path(output_dir, paste0(ensgid, ".tsv"))
  A <- dplyr::select(split_data[[ensgid]], -ENSGID)
  readr::write_tsv(A, file = file_name)
  NULL
}, future.seed = TRUE)

plan(sequential)

combined_data<-vroom("liver_eqtl.txt.gz")

suppressPackageStartupMessages({
  library(data.table)
})

dt <- as.data.table(combined_data)

output_dir <- "split_files"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
combined_data <- as.data.table(combined_data)
cores <- max(1, parallel::detectCores() - 1)

setDTthreads(cores)

combined_data[, {
  out <- .SD
  out[, ENSGID := NULL]
  fwrite(out, file = file.path(output_dir, paste0(.BY$ENSGID, ".tsv")),
         sep = "\t", quote = FALSE)
  NULL
}, by = ENSGID]
