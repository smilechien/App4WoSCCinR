
# app.R - WoS stable pipeline (Routine 1/2/3 -> long16 -> metatable -> summary + PNG)
# - Keeps homepage layout simple; adds top tabs: ReadMe, Log, then step-by-step outputs.
# - No markdown/litedown dependency.

AA_V14_VERSION <- 'aa_v14_fixed_10domains_regexsafe_v1'

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readxl)
  library(ggplot2)
  library(htmltools)
  library(base64enc)
  library(igraph)
  library(visNetwork)
})
options(shiny.maxRequestSize = 30 * 1024^2)

# ---- Paths / modules ----
APP_DIR <- getwd()
source(file.path(APP_DIR, "ipmodule.R"), local = FALSE)
source(file.path(APP_DIR, "appAAC_module_user.R"),   local = FALSE)
# ---- Optional: load SSplot panel renderer (same as apppubmed.R) ----
if (file.exists("renderSSplot.R")) {
  try(source("renderSSplot.R", local = FALSE, encoding = "UTF-8"), silent = TRUE)
  cat("[BOOT] sourced: renderSSplot.R\n")
} else {
  cat("[BOOT] renderSSplot.R not found; SSplot will use fallback message.\n")
}


# ---- Unified cluster palette (shared by Network / Sankey / SS / Chord) ----
.cluster_base_palette <- c(
  "#96ab03", "#6495ED", "#FF7F50", "#CCCCFF", "#9FE2BF",
  "#40E0D0", "#FFBF00", "#f00000", "#7ABF1B", "#C6CCBD",
  "#8A2BE2", "#00CED1", "#FFD700", "#FF69B4", "#A52A2A",
  "#2E8B57", "#1E90FF", "#FF8C00", "#708090"
)

.cluster_palette_map <- function(cluster_levels){
  clv <- unique(as.character(cluster_levels))
  clv <- clv[!is.na(clv) & nzchar(clv)]
  if (!length(clv)) clv <- "C1"
  cols <- .cluster_base_palette[((seq_along(clv)-1) %% length(.cluster_base_palette)) + 1]
  setNames(cols, clv)
}

.safe_addlog <- function(msg){
  # never crash app due to logging
  tryCatch({
    if (exists("addlog", mode="function")) addlog(msg)
  }, error=function(e) NULL)
}




# ---- helpers: Year ↔ Domain bipartite network (Top-N items, weight>1) ----
.build_year_domain_bipartite <- function(mt, domain, top_n = 20, min_w = 2, last_n_years = 10) {
  if (is.null(mt) || !is.data.frame(mt) || !nrow(mt) || !("year" %in% names(mt))) {
    return(list(nodes = data.frame(), edges = data.frame(from=character(), to=character(), value=numeric(), stringsAsFactors=FALSE)))
  }

  y <- suppressWarnings(as.integer(as.character(mt$year)))
  y <- y[is.finite(y)]
  if (!length(y)) {
    return(list(nodes = data.frame(), edges = data.frame(from=character(), to=character(), value=numeric(), stringsAsFactors=FALSE)))
  }
  ymax <- max(y, na.rm = TRUE)
  ymin <- ymax - (last_n_years - 1)

  mt2 <- mt
  mt2$year <- suppressWarnings(as.integer(as.character(mt2$year)))
  mt2 <- mt2[is.finite(mt2$year) & mt2$year >= ymin & mt2$year <= ymax, , drop = FALSE]
  if (!nrow(mt2)) {
    return(list(nodes = data.frame(), edges = data.frame(from=character(), to=character(), value=numeric(), stringsAsFactors=FALSE)))
  }

  if (domain == "Year") {
    yrs <- sort(unique(mt2$year))
    nodes <- data.frame(id=as.character(yrs), label=as.character(yrs), group="Year", value=1, stringsAsFactors=FALSE)
    return(list(nodes=nodes, edges=data.frame(from=character(), to=character(), value=numeric(), stringsAsFactors=FALSE)))
  }

  # pick item column
  pick_col <- function(cands) { hit <- intersect(cands, names(mt2)); if (length(hit)) hit[1] else NA_character_ }

  item_col <- NA_character_
  if (domain == "Journal") item_col <- pick_col(c("journal","Journal","SO","Source.Title","Source Title"))
  if (domain == "DocumentType") item_col <- pick_col(c("DocumentType","DT","document_type"))

  if (is.na(item_col)) {
    return(list(nodes = data.frame(), edges = data.frame(from=character(), to=character(), value=numeric(), stringsAsFactors=FALSE)))
  }

  rows <- data.frame(
    year = mt2$year,
    item = as.character(mt2[[item_col]]),
    stringsAsFactors = FALSE
  )
  rows <- rows[is.finite(rows$year) & nzchar(rows$item), , drop = FALSE]
  if (!nrow(rows)) {
    return(list(nodes = data.frame(), edges = data.frame(from=character(), to=character(), value=numeric(), stringsAsFactors=FALSE)))
  }

  suppressWarnings({
    edges <- dplyr::as_tibble(rows) |>
      dplyr::group_by(.data$year, .data$item) |>
      dplyr::summarise(value = dplyr::n(), .groups = "drop")
  })

  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  edges <- edges[is.finite(edges$value) & edges$value >= min_w, , drop = FALSE]
  if (!nrow(edges)) {
    return(list(nodes = data.frame(), edges = data.frame(from=character(), to=character(), value=numeric(), stringsAsFactors=FALSE)))
  }

  # Top-N items by total weight
  tot <- aggregate(value ~ item, edges, sum)
  tot <- tot[order(-tot$value), , drop=FALSE]
  top_items <- head(tot$item, top_n)

  edges <- edges[edges$item %in% top_items, , drop=FALSE]
  edges <- edges[order(edges$year, -edges$value), , drop=FALSE]

  yrs <- sort(unique(edges$year))
  nodes <- rbind(
    data.frame(id=as.character(yrs), label=as.character(yrs), group="Year", value=1, stringsAsFactors=FALSE),
    data.frame(id=as.character(top_items), label=as.character(top_items), group=domain, value=1, stringsAsFactors=FALSE)
  )
  nodes <- nodes[!duplicated(nodes$id), , drop=FALSE]

  list(
    nodes = nodes,
    edges = data.frame(from=as.character(edges$year), to=as.character(edges$item), value=as.numeric(edges$value), stringsAsFactors=FALSE)
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

.safe_source <- function(path){
  if (file.exists(path)) try(source(path, local = FALSE), silent = TRUE)
}




.safe_source(file.path("R", "wos_vba_port.R"))
.safe_source(file.path("R", "make_final_from_woslong.R"))
.safe_source(file.path("R", "wos_to_pubmed_glue.R"))
# ---- Load FLCA module (standalone) ----
if (file.exists('flca_ms_sil_module.R')) {
  .safe_source('flca_ms_sil_module.R')
}
if (exists('run_flca_ms_sil_runner', mode='function')) run_flca_ms_sil <- run_flca_ms_sil_runner


# ---------- PubMed visuals integration (Phase 2: minimal, no optimization) ----------
PUBMED_DIR <- file.path(getwd(), "pubmed")



# Load PubMed helper scripts if present (from provided PubMed visuals zip)
if (dir.exists(PUBMED_DIR)) {
  .safe_source(file.path(PUBMED_DIR, "flca_ms_sil_module.R"))
  .safe_source(file.path(PUBMED_DIR, "helper_ss_patch.R"))
  .safe_source(file.path(PUBMED_DIR, "silhouette.R"))
  .safe_source(file.path(PUBMED_DIR, "renderSSplot.R"))
  .safe_source(file.path(PUBMED_DIR, "kano.R"))
  .safe_source(file.path(PUBMED_DIR, "sankey.R"))
  .safe_source(file.path(PUBMED_DIR, "pubmed_parse_biblio.R"))
}

# Keep the app-local SSplot renderer highest priority when present.
if (file.exists(file.path(APP_DIR, "renderSSplot.R"))) {
  try(source(file.path(APP_DIR, "renderSSplot.R"), local = FALSE, encoding = "UTF-8"), silent = TRUE)
  cat("[BOOT] re-sourced app-local renderSSplot.R after PubMed helpers\n")
}

# Use renderSSplot.R directly; do not override render_panel() here.

# Override local AAC formula to use value (document count) first.
.calc_cluster_aac <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(NA_real_)
  base_all <- if ("value" %in% names(df)) {
    suppressWarnings(as.numeric(df$value))
  } else if ("value2" %in% names(df)) {
    suppressWarnings(as.numeric(df$value2))
  } else if ("sil_width" %in% names(df)) {
    suppressWarnings(as.numeric(df$sil_width))
  } else {
    numeric(0)
  }
  v <- base_all[is.finite(base_all) & base_all >= 1]
  v <- sort(v, decreasing = TRUE)
  if (length(v) == 0) return(NA_real_)
  if (length(v) == 1) return(0.5)
  if (length(v) == 2) {
    r <- v[1] / v[2]
    return(if (is.finite(r) && r > 0) r / (1 + r) else NA_real_)
  }
  r <- (v[1] / v[2]) / (v[2] / v[3])
  if (is.finite(r) && r > 0) r / (1 + r) else NA_real_
}




# ---------- FALLBACK: ensure FLCA + Major Sampling + Silhouette always available ----------
# If the full PubMed FLCA module is not present, we provide a lightweight fallback so the app
# still produces Top20 + cluster colors + SS/a* columns for visNetwork.
if (!exists("run_flca_ms_sil", mode = "function")) {
  run_flca_ms_sil <- function(nodes_full, edges_full, cfg = list(), verbose = FALSE) {
    if (is.null(nodes_full) || !is.data.frame(nodes_full) || nrow(nodes_full) == 0 ||
        is.null(edges_full) || !is.data.frame(edges_full) || nrow(edges_full) == 0) {
      return(list(modes = nodes_full, data = edges_full))
    }

    # Normalize columns
    nodes <- nodes_full
    if (!("name" %in% names(nodes))) {
      nm <- intersect(c("Name","ID","id","label"), names(nodes))
      nodes$name <- if (length(nm)) as.character(nodes[[nm[1]]]) else as.character(seq_len(nrow(nodes)))
    }
    if (!("value" %in% names(nodes))) nodes$value <- 1
    if (!("value2" %in% names(nodes))) nodes$value2 <- nodes$value
    if (!("value2_strength" %in% names(nodes))) nodes$value2_strength <- nodes$value2
    if (!("n_pubmed" %in% names(nodes))) nodes$n_pubmed <- nodes$value
    if (!("single_author" %in% names(nodes))) nodes$single_author <- pmax(0, nodes$value - nodes$value2)
    # Ensure FLCA/node summary columns exist for DT tables
    if (!("carac" %in% names(nodes))) nodes$carac <- NA_character_
    if (!("ss" %in% names(nodes))) {
      if ("SS" %in% names(nodes)) nodes$ss <- nodes$SS else nodes$ss <- NA_real_
    }
    if (!("a*" %in% names(nodes))) {
      if ("a_star1" %in% names(nodes)) nodes[["a*"]] <- nodes$a_star1 else nodes[["a*"]] <- NA_real_
    }

    ed <- edges_full
    leader_col   <- intersect(c("Leader","leader","from","From"), names(ed))
    follower_col <- intersect(c("Follower","follower","to","To"), names(ed))
    wcd_col      <- intersect(c("WCD","wcd","weight","Weight","value"), names(ed))
    if (!length(leader_col) || !length(follower_col)) {
      return(list(modes = nodes_full, data = edges_full))
    }
    a <- as.character(ed[[leader_col[1]]])
    b <- as.character(ed[[follower_col[1]]])
    w <- if (length(wcd_col)) suppressWarnings(as.numeric(ed[[wcd_col[1]]])) else rep(1, length(a))
    w[!is.finite(w)] <- 1
    g <- igraph::graph_from_data_frame(data.frame(from=a, to=b, weight=w, stringsAsFactors=FALSE), directed=FALSE)

    # --- Major Sampling (Top20) fallback: keep highest "value" nodes
    target_n <- cfg$target_n %||% 20
    nodes_ord <- nodes[order(-suppressWarnings(as.numeric(nodes$value)), nodes$name), , drop=FALSE]
    keep <- head(nodes_ord$name, min(target_n, nrow(nodes_ord)))
    g2 <- igraph::induced_subgraph(g, vids = intersect(keep, igraph::V(g)$name))

    # --- Cluster fallback: Louvain on weighted undirected graph
    mem <- tryCatch({
      if (igraph::ecount(g2) > 0) igraph::membership(igraph::cluster_louvain(g2, weights=igraph::E(g2)$weight))
      else setNames(seq_along(igraph::V(g2)), igraph::V(g2)$name)
    }, error=function(e) setNames(seq_along(igraph::V(g2)), igraph::V(g2)$name))
    carac <- as.integer(factor(mem, levels=unique(mem)))
    carac <- setNames(carac, names(mem))

    modes <- nodes[nodes$name %in% names(mem), , drop=FALSE]
    modes$carac <- as.integer(carac[modes$name])

    # --- Silhouette fallback (rough): based on shortest-path distance with 1/weight
    modes$ssi <- 0
    if (igraph::vcount(g2) >= 3 && igraph::ecount(g2) > 0) {
      d <- igraph::distances(g2, weights = 1 / pmax(1e-9, igraph::E(g2)$weight))
      # silhouette per node
      cl <- modes$carac
      vnames <- modes$name
      for (i in seq_along(vnames)) {
        di <- d[vnames[i], ]
        # a(i): mean distance within cluster
        in_idx <- which(cl == cl[i] & vnames != vnames[i])
        a_i <- if (length(in_idx)) mean(di[vnames[in_idx]], na.rm=TRUE) else NA_real_
        # b(i): min mean distance to other clusters
        other_cls <- setdiff(unique(cl), cl[i])
        b_i <- NA_real_
        if (length(other_cls)) {
          b_i <- min(vapply(other_cls, function(cc){
            idx <- which(cl == cc)
            mean(di[vnames[idx]], na.rm=TRUE)
          }, FUN.VALUE = numeric(1)), na.rm=TRUE)
        }
        s_i <- if (is.finite(a_i) && is.finite(b_i) && max(a_i, b_i) > 0) (b_i - a_i) / max(a_i, b_i) else 0
        modes$ssi[i] <- s_i
      }
    }
    # a* fallback: scaled value2_strength within cluster
    modes$a_star1 <- ave(suppressWarnings(as.numeric(modes$value2_strength)), modes$carac,
                        FUN=function(z){ if (all(!is.finite(z))) return(rep(0, length(z)));
                          z[!is.finite(z)] <- 0; if (max(z) == 0) rep(0, length(z)) else z/max(z) })

    # Return edges restricted to kept nodes
    ed2 <- igraph::as_data_frame(g2, what="edges")
    data <- data.frame(Leader=ed2$from, Follower=ed2$to, WCD=ed2$weight, stringsAsFactors=FALSE)
    list(modes = modes, data = data)
  }
}

if (!exists("build_domain_flca_top20", mode = "function")) {
  build_domain_flca_top20 <- function(edges, term_list, seed = 123) {
    # edges: data.frame(from,to,weight) or Leader/Follower/WCD
    if (is.null(edges) || !is.data.frame(edges) || nrow(edges) == 0) {
      return(list(nodes = data.frame(), data = data.frame()))
    }
    # Convert to Leader/Follower/WCD for reuse
    ed <- edges
    if (!("Leader" %in% names(ed))) {
      if (all(c("from","to") %in% names(ed))) {
        ed$Leader <- as.character(ed$from)
        ed$Follower <- as.character(ed$to)
      }
    }
    if (!("WCD" %in% names(ed))) {
      if ("weight" %in% names(ed)) ed$WCD <- suppressWarnings(as.numeric(ed$weight)) else ed$WCD <- 1
    }
    ed <- ed[, c("Leader","Follower","WCD"), drop=FALSE]
    # Nodes: frequency from term_list, plus degree from edges
    all_terms <- unlist(lapply(term_list, unique), use.names = FALSE)
    freq <- sort(table(all_terms), decreasing=TRUE)
    node_names <- sort(unique(c(ed$Leader, ed$Follower)))
    nodes <- data.frame(name=node_names, stringsAsFactors=FALSE)
    nodes$value <- as.integer(freq[nodes$name]); nodes$value[is.na(nodes$value)] <- 0L
    deg <- table(c(ed$Leader, ed$Follower))
    nodes$value2 <- as.integer(deg[nodes$name]); nodes$value2[is.na(nodes$value2)] <- 0L
    nodes$value2_strength <- nodes$value2
    nodes$n_pubmed <- nodes$value
    nodes$single_author <- 0L

    cfg <- list(target_n = 20)
    res <- run_flca_ms_sil(nodes, ed, cfg = cfg, verbose = FALSE)
    list(nodes = res$modes, data = res$data)
  }
}


.load_pubmed_nodes_edges <- function(){
  nodes_path <- file.path(PUBMED_DIR, "nodes.csv")
  edges_path <- file.path(PUBMED_DIR, "data_edges.csv")
  if (!file.exists(nodes_path) || !file.exists(edges_path)) {
    return(list(nodes = NULL, edges = NULL))
  }
  nodes <- try(read.csv(nodes_path, stringsAsFactors = FALSE), silent = TRUE)
  edges <- try(read.csv(edges_path, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(nodes, "try-error")) nodes <- NULL
  if (inherits(edges, "try-error")) edges <- NULL
  list(nodes = nodes, edges = edges)
}

.norm_key <- function(x){
  x <- tolower(trimws(as.character(x)))
  x <- gsub("\\s+", " ", x)
  x
}

.pubmed_add_carac <- function(nodes){
  if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) == 0) return(nodes)
  if (!all(c("value","value2") %in% names(nodes))) return(nodes)
  nodes$value  <- suppressWarnings(as.numeric(nodes$value))
  nodes$value2 <- suppressWarnings(as.numeric(nodes$value2))
  mx <- mean(nodes$value2, na.rm=TRUE)
  my <- mean(nodes$value,  na.rm=TRUE)
  nodes$carac <- ifelse(nodes$value2 >= mx & nodes$value >= my, "A",
                 ifelse(nodes$value2 >= mx & nodes$value <  my, "B",
                 ifelse(nodes$value2 <  mx & nodes$value >= my, "C", "D")))
  nodes
}

# ---------- helpers ----------
.trim_na <- function(x){
  x <- ifelse(is.na(x), "", as.character(x))
  trimws(x)
}

# Standardize / clean country names for mapping and labels
.clean_country_label <- function(x){
  x <- .trim_na(x)
  if (!length(x)) return(x)
  x <- gsub("#[0-9]+\\b", "", x, perl = TRUE)           # remove cluster suffix like #2
  x <- gsub("\\([^)]*\\)", "", x, perl = TRUE)       # remove trailing (C2) etc.
  x <- gsub("[;|]+.*$", "", x, perl = TRUE)              # keep first country-like token block only
  x <- trimws(gsub("\\s+", " ", x, perl = TRUE))

  alias_map <- list(
    "United States" = c("united states of america","united states","u s a","u s","usa","us"),
    "United Kingdom" = c("united kingdom","u k","uk","great britain","britain","england","scotland","wales","northern ireland","north ireland"),
    "Germany" = c("germany"),
    "Italy" = c("italy"),
    "Norway" = c("norway"),
    "India" = c("india"),
    "Brazil" = c("brazil"),
    "Portugal" = c("portugal"),
    "Sweden" = c("sweden"),
    "Finland" = c("finland"),
    "France" = c("france"),
    "Spain" = c("spain"),
    "Netherlands" = c("netherlands"),
    "Belgium" = c("belgium"),
    "Switzerland" = c("switzerland"),
    "Austria" = c("austria"),
    "Denmark" = c("denmark"),
    "Ireland" = c("ireland","republic of ireland"),
    "China" = c("china","peoples republic of china","people s republic of china","pr china","p r china"),
    "Taiwan" = c("taiwan","roc","r o c","republic of china"),
    "Hong Kong" = c("hong kong"),
    "Macau" = c("macau","macao"),
    "Japan" = c("japan"),
    "South Korea" = c("south korea","korea","republic of korea"),
    "Singapore" = c("singapore"),
    "Australia" = c("australia"),
    "New Zealand" = c("new zealand"),
    "Canada" = c("canada"),
    "Mexico" = c("mexico")
  )

  norm <- tolower(x)
  norm <- gsub("[[:punct:]]", " ", norm)
  norm <- gsub("\\s+", " ", trimws(norm))

  out <- x
  for (i in seq_along(norm)) {
    ni <- norm[i]
    if (!nzchar(ni)) { out[i] <- ""; next }

    found <- NA_character_
    for (canon in names(alias_map)) {
      aliases <- alias_map[[canon]]
      hit <- vapply(aliases, function(a) {
        a2 <- gsub("\\s+", " ", trimws(tolower(a)))
        grepl(paste0("(^| )", a2, "$"), ni, perl = TRUE)
      }, logical(1))
      if (any(hit)) { found <- canon; break }
    }

    if (!is.na(found)) {
      out[i] <- found
    } else {
      # fallback: if comma-separated, keep last segment; then remove leading state/zip noise
      xi <- trimws(tail(strsplit(x[i], ",", fixed = TRUE)[[1]], 1))
      xi <- gsub("^[A-Z]{2}\\s+[0-9]{3,6}\\s+", "", xi, perl = TRUE)
      xi <- gsub("^[A-Z]{2}\\s+", "", xi, perl = TRUE)
      xi <- trimws(gsub("\\s+", " ", xi, perl = TRUE))
      out[i] <- xi
    }
  }
  out
}

.std_country <- function(x){
  x <- .clean_country_label(x)
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | !nzchar(x)] <- ""
  x
}

.us_state_map <- local({
  nm <- c(
    "AL"="Alabama","AK"="Alaska","AZ"="Arizona","AR"="Arkansas","CA"="California","CO"="Colorado",
    "CT"="Connecticut","DE"="Delaware","FL"="Florida","GA"="Georgia","HI"="Hawaii","ID"="Idaho",
    "IL"="Illinois","IN"="Indiana","IA"="Iowa","KS"="Kansas","KY"="Kentucky","LA"="Louisiana",
    "ME"="Maine","MD"="Maryland","MA"="Massachusetts","MI"="Michigan","MN"="Minnesota","MS"="Mississippi",
    "MO"="Missouri","MT"="Montana","NE"="Nebraska","NV"="Nevada","NH"="New Hampshire","NJ"="New Jersey",
    "NM"="New Mexico","NY"="New York","NC"="North Carolina","ND"="North Dakota","OH"="Ohio","OK"="Oklahoma",
    "OR"="Oregon","PA"="Pennsylvania","RI"="Rhode Island","SC"="South Carolina","SD"="South Dakota",
    "TN"="Tennessee","TX"="Texas","UT"="Utah","VT"="Vermont","VA"="Virginia","WA"="Washington",
    "WV"="West Virginia","WI"="Wisconsin","WY"="Wyoming","DC"="District Of Columbia"
  )
  as.list(nm)
})

.china_region_aliases <- list(
  "Beijing"=c("beijing","peking","北京市","北京"),
  "Shanghai"=c("shanghai","上海市","上海"),
  "Tianjin"=c("tianjin","天津市","天津"),
  "Chongqing"=c("chongqing","重庆市","重慶市","重庆","重慶"),
  "Guangdong"=c("guangdong","广东省","廣東省","guangzhou","shenzhen","dongguan","foshan","zhuhai","guangdong province","canton"),
  "Guangxi"=c("guangxi","广西","廣西","guangxi zhuang autonomous region","nanning"),
  "Fujian"=c("fujian","福建省","福建","xiamen","fuzhou","quanzhou"),
  "Zhejiang"=c("zhejiang","浙江省","浙江","hangzhou","ningbo","wenzhou"),
  "Jiangsu"=c("jiangsu","江苏省","江蘇省","江苏","江蘇","nanjing","suzhou","wuxi","changzhou"),
  "Shandong"=c("shandong","山东省","山東省","山东","山東","jinan","qingdao","yantai"),
  "Henan"=c("henan","河南省","河南","zhengzhou"),
  "Hebei"=c("hebei","河北省","河北","shijiazhuang"),
  "Hubei"=c("hubei","湖北省","湖北","wuhan"),
  "Hunan"=c("hunan","湖南省","湖南","changsha"),
  "Sichuan"=c("sichuan","四川省","四川","chengdu"),
  "Shaanxi"=c("shaanxi","陕西省","陝西省","陕西","陝西","xian","xi'an"),
  "Shanxi"=c("shanxi","山西省","山西","taiyuan"),
  "Liaoning"=c("liaoning","辽宁省","遼寧省","辽宁","遼寧","shenyang","dalian"),
  "Jilin"=c("jilin","吉林省","吉林","changchun"),
  "Heilongjiang"=c("heilongjiang","黑龙江省","黑龍江省","黑龙江","黑龍江","harbin"),
  "Anhui"=c("anhui","安徽省","安徽","hefei"),
  "Jiangxi"=c("jiangxi","江西省","江西","nanchang"),
  "Yunnan"=c("yunnan","云南省","雲南省","云南","雲南","kunming"),
  "Guizhou"=c("guizhou","贵州省","貴州省","贵州","貴州","guiyang"),
  "Gansu"=c("gansu","甘肃省","甘肅省","甘肃","甘肅","lanzhou"),
  "Qinghai"=c("qinghai","青海省","青海","xining"),
  "Inner Mongolia"=c("inner mongolia","nei mongol","内蒙古","內蒙古","hohhot"),
  "Xinjiang"=c("xinjiang","新疆","urumqi","ü rümqi","u rumqi"),
  "Tibet"=c("tibet","xizang","西藏","lhasa"),
  "Ningxia"=c("ningxia","宁夏","寧夏","yinchuan"),
  "Hainan"=c("hainan","海南省","海南","haikou"),
  "Hong Kong"=c("hong kong","香港","hongkong","hk"),
  "Macau"=c("macau","macao","澳门","澳門"),
  "Taiwan"=c("taiwan","taipei","new taipei","taichung","tainan","kaohsiung","hsinchu","keelung","taoyuan","chiayi","miaoli","changhua","nantou","yunlin","pingtung","yilan","hualien","taitung","penghu","台灣","臺灣","台北","臺北","新北","台中","臺中","台南","臺南","高雄","新竹","基隆","桃園","嘉義","苗栗","彰化","南投","雲林","屏東","宜蘭","花蓮","台東","臺東","澎湖")
)

.normalize_region_value <- function(country = "", region = "", institute = "", dept = "", extra = ""){
  country <- .std_country(country)
  txt <- paste(region %||% "", institute %||% "", dept %||% "", extra %||% "", sep = " | ")
  txt <- trimws(gsub("\\s+", " ", as.character(txt), perl = TRUE))
  txt_low <- tolower(txt)
  reg <- trimws(as.character(region %||% ""))
  reg_low <- tolower(reg)

  if (!nzchar(txt_low) && !nzchar(country)) return("")

  if (identical(country, "United States")) {
    if (nzchar(reg)) {
      rg <- gsub("[^A-Za-z ]", " ", reg)
      rg <- trimws(gsub("\\s+", " ", rg))
      rgu <- toupper(rg)
      if (nchar(rgu) == 2 && !is.null(.us_state_map[[rgu]])) return(.us_state_map[[rgu]])
      title_rg <- tools::toTitleCase(tolower(rg))
      if (title_rg %in% unlist(.us_state_map, use.names = FALSE)) return(title_rg)
    }
    for (abbr in names(.us_state_map)) {
      pat <- paste0("(^|[ ,;|()])", abbr, "([ ,;|()0-9]|$)")
      if (grepl(pat, txt, perl = TRUE)) return(.us_state_map[[abbr]])
    }
    for (nm in unlist(.us_state_map, use.names = FALSE)) {
      if (grepl(paste0("(^|[ ,;|()])", tolower(nm), "([ ,;|()]|$)"), txt_low, perl = TRUE)) return(nm)
    }
    return("")
  }

  if (country %in% c("China","Hong Kong","Macau","Taiwan")) {
    for (canon in names(.china_region_aliases)) {
      aliases <- .china_region_aliases[[canon]]
      hit <- any(vapply(aliases, function(a) {
        a_low <- tolower(trimws(as.character(a)))
        grepl(a_low, txt_low, fixed = TRUE)
      }, logical(1)))
      if (hit) return(canon)
    }
    if (nzchar(reg)) {
      rg2 <- trimws(gsub("\\s+", " ", reg, perl = TRUE))
      # if region itself already looks like a Chinese province/city, keep it for later mapping
      if (grepl("省|市|自治区|自治區|特别行政区|特別行政區", rg2)) return(rg2)
    }
    return("")
  }

  if (country %in% c("United Kingdom","Canada","Australia","Japan","South Korea","India","Brazil","Germany","France","Italy","Spain","Portugal","Sweden","Finland","Netherlands","Belgium","Switzerland","Austria","Denmark","Ireland","New Zealand","Mexico")) {
    if (nzchar(reg)) return(tools::toTitleCase(tolower(reg)))
    return("")
  }

  if (nzchar(reg)) return(tools::toTitleCase(tolower(reg)))
  ""
}

.apply_region_normalization <- function(df){
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(df)
  for (side in c("FA","CA")) {
    ccol <- paste0(side, "country")
    rcol <- paste0(side, "region")
    icol <- paste0(side, "institute")
    dcol <- if (side == "FA") "FADept" else "CAdept"
    if (!(ccol %in% names(df))) df[[ccol]] <- ""
    if (!(rcol %in% names(df))) df[[rcol]] <- ""
    if (!(icol %in% names(df))) df[[icol]] <- ""
    if (!(dcol %in% names(df))) df[[dcol]] <- ""
    extra <- if ("citation" %in% names(df)) df$citation else ""
    df[[rcol]] <- vapply(seq_len(nrow(df)), function(i) {
      .normalize_region_value(df[[ccol]][i], df[[rcol]][i], df[[icol]][i], df[[dcol]][i], extra[i])
    }, character(1))
  }
  df
}

.domain_doc_counts_from_metatable <- function(meta_df, dom){
  if (is.null(meta_df) || !is.data.frame(meta_df) || !nrow(meta_df)) return(setNames(numeric(), character()))
  dom <- as.character(dom %||% "")
  get_terms <- function(v) unique(trimws(as.character(v[!is.na(v) & nzchar(trimws(as.character(v)))])))
  counts <- character()
  if (dom == "Keyword") {
    cols <- grep("^A[0-9]+$", names(meta_df), value = TRUE)
    if (!length(cols)) return(setNames(numeric(), character()))
    per_row <- apply(meta_df[, cols, drop = FALSE], 1, get_terms)
  } else if (dom %in% c("Region", "State/Province")) {
    cols <- intersect(c("FAregion","CAregion"), names(meta_df))
    if (!length(cols)) return(setNames(numeric(), character()))
    per_row <- lapply(seq_len(nrow(meta_df)), function(i) get_terms(unlist(meta_df[i, cols, drop = TRUE], use.names = FALSE)))
  } else if (dom == "Country") {
    cols <- intersect(c("FAcountry","CAcountry"), names(meta_df))
    if (!length(cols)) return(setNames(numeric(), character()))
    per_row <- lapply(seq_len(nrow(meta_df)), function(i) get_terms(unlist(meta_df[i, cols, drop = TRUE], use.names = FALSE)))
  } else if (dom == "Department") {
    cols <- intersect(c("FADept","CAdept"), names(meta_df))
    if (!length(cols)) return(setNames(numeric(), character()))
    per_row <- lapply(seq_len(nrow(meta_df)), function(i) get_terms(unlist(meta_df[i, cols, drop = TRUE], use.names = FALSE)))
  } else if (dom == "Institute") {
    cols <- intersect(c("FAinstitute","CAinstitute"), names(meta_df))
    if (!length(cols)) return(setNames(numeric(), character()))
    per_row <- lapply(seq_len(nrow(meta_df)), function(i) get_terms(unlist(meta_df[i, cols, drop = TRUE], use.names = FALSE)))
  } else if (dom == "Author") {
    cols <- intersect(c("FAauthor","CAauthor"), names(meta_df))
    if (!length(cols)) return(setNames(numeric(), character()))
    per_row <- lapply(seq_len(nrow(meta_df)), function(i) get_terms(unlist(meta_df[i, cols, drop = TRUE], use.names = FALSE)))
  } else return(setNames(numeric(), character()))
  all_terms <- unlist(per_row, use.names = FALSE)
  if (!length(all_terms)) return(setNames(numeric(), character()))
  tb <- table(all_terms)
  setNames(as.numeric(tb), names(tb))
}

h_index_cites <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_integer_)
  x <- sort(x, decreasing = TRUE)
  sum(x >= seq_along(x))
}

split_kp <- function(x){
  x <- .trim_na(x)
  if (!nzchar(x)) return(character(0))
  x <- gsub("\\s*;+\\s*", ";", x)
  x <- gsub("^;+|;+$", "", x)
  toks <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  toks <- trimws(toks)
  toks[toks != ""]
}

kp_to_Acols <- function(kp_vec, max_k = 13){
  toks <- lapply(kp_vec, split_kp)
  mat <- t(vapply(toks, function(v){
    v <- v[seq_len(min(length(v), max_k))]
    if (length(v) < max_k) v <- c(v, rep("", max_k - length(v)))
    v
  }, FUN.VALUE = character(max_k)))
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(out) <- paste0("A", seq_len(max_k))
  out
}

ensure_col <- function(df, nm, val){
  if (!nm %in% names(df)) df[[nm]] <- val
  df
}

# h-index for a vector of counts
h_index <- function(counts){
  counts <- sort(as.integer(counts), decreasing = TRUE)
  if (!length(counts)) return(0L)
  h <- 0L
  for (i in seq_along(counts)){
    if (counts[i] >= i) h <- i else break
  }
  h
}


# ---------- Summary helpers (10 domains; VBA-style counting) ----------
# n: distinct article_id + element (so FA/CA same term in same paper counted once)
# h: h-index on per-article element occurrence counts (usually 1s, but kept for consistency)
h_index <- function(x){
  x <- sort(as.integer(x[!is.na(x)]), decreasing = TRUE)
  if (!length(x)) return(0L)
  sum(x >= seq_along(x))
}

compute_aac <- function(values){
  v <- suppressWarnings(as.numeric(values))
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) return(NA_real_)
  v <- sort(v, decreasing = TRUE)

  k <- length(v)
  if (k >= 3) {
    v1 <- v[1]; v2 <- v[2]; v3 <- v[3]
    if (v2 == 0 || v3 == 0) return(NA_real_)
    r <- (v1 / v2) / (v2 / v3)
    return(r / (1 + r))
  } else if (k == 2) {
    r <- 1
    return(r / (1 + r))
  } else { # k == 1
    return(0)
  }
}

.split_semis <- function(x){
  x <- as.character(x)
  x[is.na(x)] <- ""
  out <- strsplit(x, "\\s*;\\s*")
  out
}

.long_from_cols <- function(df, cols, domain, cite_col = "citation"){
  cols <- cols[cols %in% names(df)]
  if (!length(cols)) {
    return(tibble::tibble(Domain=character(), Element=character(), n=integer(), h=integer(), AAC=numeric()))
  }

  has_cite <- cite_col %in% names(df)

  tmp <- dplyr::bind_rows(lapply(cols, function(cl){
    tibble::tibble(
      article_id = as.character(df$article_id),
      Element    = as.character(df[[cl]]),
      cites      = if (has_cite) suppressWarnings(as.numeric(df[[cite_col]])) else NA_real_
    )
  })) %>%
    dplyr::mutate(Element = .trim_na(Element)) %>%
    dplyr::filter(!is.na(Element) & Element != "")

  if (!nrow(tmp)) {
    return(tibble::tibble(Domain=domain, Element=character(), n=integer(), h=integer(), AAC=numeric()))
  }

  # IMPORTANT: 同一篇文章同一個 element 只算一次，但要保留 cites
  tmp <- tmp %>% dplyr::distinct(article_id, Element, .keep_all = TRUE)

  n_tbl <- tmp %>% dplyr::count(Element, name="n") %>% dplyr::arrange(dplyr::desc(n), Element)

  if (has_cite) {
    h_tbl <- tmp %>%
      dplyr::group_by(Element) %>%
      dplyr::summarise(h = h_index_cites(cites), .groups="drop")
  } else {
    # fallback（沒 citation 才走這個，否則 h 幾乎都會是 1）
    h_tbl <- tmp %>%
      dplyr::group_by(Element, article_id) %>%
      dplyr::summarise(k = dplyr::n(), .groups="drop") %>%
      dplyr::group_by(Element) %>%
      dplyr::summarise(h = h_index(k), .groups="drop")
  }

  out <- n_tbl %>%
    dplyr::left_join(h_tbl, by="Element") %>%
    dplyr::mutate(Domain = domain) %>%
    dplyr::select(Domain, Element, n, h)

  top_vals <- out$n[seq_len(min(5, nrow(out)))]
  aac <- compute_aac(top_vals)

  out %>% dplyr::mutate(AAC = aac) %>% dplyr::slice_head(n=5)
}



.long_from_keywordsplus <- function(df, domain="Keyword"){
  # WoS keyword plus column name variants
  kp_candidates <- c("KeywordsPlus","Keywords Plus","Keywords.Plus","Keyword Plus","KeywordPlus","Keywords_Plus","KP")
  kp_candidates <- kp_candidates[kp_candidates %in% names(df)]
  if (!length(kp_candidates)) {
    return(tibble::tibble(Domain=domain, Element=character(), n=integer(), h=integer(), AAC=numeric()))
  }
  col <- kp_candidates[1]

  cite_col <- detect_cite_col(df)

  sp <- .split_semis(df[[col]])
  tmp <- dplyr::bind_rows(lapply(seq_along(sp), function(i){
    if (!length(sp[[i]])) return(NULL)
    tibble::tibble(
      article_id = as.character(df$article_id[i]),
      Element    = sp[[i]],
      cites      = if (!is.na(cite_col)) suppressWarnings(as.numeric(df[[cite_col]][i])) else NA_real_
    )
  })) %>%
    dplyr::mutate(Element = .trim_na(Element)) %>%
    dplyr::filter(!is.na(Element) & Element != "") %>%
    dplyr::distinct(article_id, Element, .keep_all = TRUE)

  if (!nrow(tmp)) {
    return(tibble::tibble(Domain=domain, Element=character(), n=integer(), h=integer(), AAC=numeric()))
  }

  n_tbl <- tmp %>% dplyr::count(Element, name="n") %>% dplyr::arrange(dplyr::desc(n), Element)

  if (!is.na(cite_col)) {
    h_tbl <- tmp %>% dplyr::group_by(Element) %>% dplyr::summarise(h = h_index_cites(cites), .groups="drop")
  } else {
    # if no citation column, fall back to occurrence-based h (will usually be 1) rather than NA
    h_tbl <- tmp %>%
      dplyr::group_by(Element, article_id) %>%
      dplyr::summarise(k = dplyr::n(), .groups="drop") %>%
      dplyr::group_by(Element) %>%
      dplyr::summarise(h = h_index(k), .groups="drop")
  }

  out <- n_tbl %>%
    dplyr::left_join(h_tbl, by="Element") %>%
    dplyr::mutate(Domain = domain) %>%
    dplyr::select(Domain, Element, n, h)

  top_vals <- out$n[seq_len(min(5, nrow(out)))]
  aac <- compute_aac(top_vals)
  out %>% dplyr::mutate(AAC = aac) %>% dplyr::slice_head(n=5)
}


build_summary10 <- function(wide32, long32){
  # Build a single "summary source" table using the CORRECT metadata from woslong32,
  # and KeywordPlus / citation info from woswide32.
  w <- as.data.frame(wide32)
  l <- as.data.frame(long32)

  # ensure article_id exists and aligned
  if (!("article_id" %in% names(w))) {
    if ("row_id" %in% names(w)) w$article_id <- w$row_id else w$article_id <- seq_len(nrow(w))
  }
  if (!("article_id" %in% names(l))) {
    if ("row_id" %in% names(l)) l$article_id <- l$row_id else l$article_id <- seq_len(nrow(l))
  }
  w$article_id <- as.character(w$article_id)
  l$article_id <- as.character(l$article_id)

  # Normalize CAauthor column name
  if ("Caauthor" %in% names(l) && !("CAauthor" %in% names(l))) names(l)[names(l)=="Caauthor"] <- "CAauthor"

  # helper to pick first existing column among candidates
  pick <- function(df, cands){
    nms <- names(df); nms_l <- tolower(nms)
    c_l <- tolower(cands)
    hit <- which(nms_l %in% c_l)
    if (length(hit)) return(nms[hit[1]])
    for (c in c_l){
      h <- which(grepl(c, nms_l, fixed = TRUE))
      if (length(h)) return(nms[h[1]])
    }
    NA_character_
  }

  # author fields (IMPORTANT: avoid the bogus 'CA' column that is actually year in woslong16)
 # author fields (IMPORTANT: avoid the bogus 'CA' column that is actually year in woslong16)
col_fa <- if ("FAauthor" %in% names(l)) "FAauthor"
          else if ("FA_author" %in% names(l)) "FA_author"
          else if ("FA Author" %in% names(l)) "FA Author"
          else NA_character_

col_ca <- if ("CAauthor" %in% names(l)) "CAauthor"
          else if ("Caauthor" %in% names(l)) "CAauthor"
          else if ("CA_author" %in% names(l)) "CA_author"
          else if ("CA Author" %in% names(l)) "CA Author"
          else NA_character_
  # year/journal/document/discipline typically come from wide32
  col_year    <- pick(w, c("Publication Year","Publication.Year","year"))
  col_journal <- pick(w, c("Journal","Source Title","SO","journal"))
  col_doc     <- pick(w, c("Document Type","Document.Type","DT","DocumentType"))
  col_d1      <- pick(w, c("Discipince1","Discipline1","Research Areas","WC"))
  col_d2      <- pick(w, c("Discipince2","Discipline2","WC2"))

  # keyword plus variants from wide32
  col_kp <- pick(w, c("Keywords Plus","Keywords.Plus","Keyword Plus","KeywordPlus","KP"))
  if (is.na(col_kp)) {
    # sometimes it was carried into long32
    col_kp <- pick(l, c("Keywords Plus","Keywords.Plus","Keyword Plus","KeywordPlus","KeywordsPlus"))
  }

  # citation / times cited column detection (include 'citation' from your exports)
  cite_col <- detect_cite_col(w)
  if (is.na(cite_col)) cite_col <- detect_cite_col(l)

  # merge wide + long on article_id
  src <- l %>%
    dplyr::select(
      article_id,
      dplyr::any_of(c("FAcountry","CAcountry","FAinstitute","CAinstitute","FADept","CAdept","FAregion","CAregion"))
    ) %>%
    dplyr::left_join(
      w %>% dplyr::select(
        article_id,
        citation = dplyr::any_of(cite_col),
        year = dplyr::any_of(col_year),
        journal = dplyr::any_of(col_journal),
        DocumentType = dplyr::any_of(col_doc),
        Discipince1 = dplyr::any_of(col_d1),
        Discipince2 = dplyr::any_of(col_d2),
        KeywordsPlus = dplyr::any_of(col_kp)
      ),
      by = "article_id"
    ) %>%
    dplyr::mutate(
      FAauthor = if (!is.na(col_fa) && col_fa %in% names(l)) as.character(l[[col_fa]]) else "",
      CAauthor = if (!is.na(col_ca) && col_ca %in% names(l)) as.character(l[[col_ca]]) else "",
      KeywordsPlus = ifelse(is.na(KeywordsPlus), "", as.character(KeywordsPlus)),
      year = ifelse(is.na(year), "", as.character(year)),
      journal = ifelse(is.na(journal), "", as.character(journal)),
      DocumentType = ifelse(is.na(DocumentType), "", as.character(DocumentType)),
      Discipince1 = ifelse(is.na(Discipince1), "", as.character(Discipince1)),
      Discipince2 = ifelse(is.na(Discipince2), "", as.character(Discipince2))
    )

  # Standardize country labels and normalize region labels
  if ("FAcountry" %in% names(src)) src$FAcountry <- .std_country(src$FAcountry)
  if ("CAcountry" %in% names(src)) src$CAcountry <- .std_country(src$CAcountry)
  src <- .apply_region_normalization(as.data.frame(src))

  # now build domain summaries (top5) from src
  dfm <- src %>% dplyr::mutate(article_id = as.character(article_id))

  s_country    <- .long_from_cols(dfm, c("FAcountry","CAcountry"), "Country")
  s_institute  <- .long_from_cols(dfm, c("FAinstitute","CAinstitute"), "Institute")
  s_department <- .long_from_cols(dfm, c("FADept","CAdept"), "Department")
  s_author     <- .long_from_cols(dfm, c("FAauthor","CAauthor"), "Author")
  s_journal    <- .long_from_cols(dfm, c("journal"), "Journal")
  s_year       <- .long_from_cols(dfm, c("year"), "Year")
  s_doctype    <- .long_from_cols(dfm, c("DocumentType"), "Document Type")
  s_discip     <- .long_from_cols(dfm, c("Discipince1","Discipince2"), "Discipince")
  s_keyword    <- .long_from_keywordsplus(dfm, domain="Keyword")
  s_region     <- .long_from_cols(dfm, c("FAregion","CAregion"), "Region")

  dplyr::bind_rows(s_country, s_institute, s_department, s_author, s_journal, s_year, s_doctype, s_discip, s_keyword, s_region)
}



plot_summary10 <- function(df_sum, total_n = NA_integer_, total_h = NA_integer_){
  domains <- c("Country","Institute","Department","Author","Journal","Year","Document Type","Discipline","Keyword","Region")

  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add=TRUE)
  par(mfrow=c(5,2))
  par(mar=c(1.2, 1.0, 2.4, 1.0))
  par(oma=c(0.5, 0.5, 2.2, 0.5))  # space for global title

  for (d in domains){
    dd <- df_sum %>% dplyr::filter(Domain %in% c(d, if (d=="Discipline") "Discipince" else d))

    plot.new()

    # domain title in red
    title(main=d, font.main=3, cex.main=1.55, line=0.5, col.main="red")
    # add items on the SAME line as the title (top margin)
    mtext("h", side = 3, line = 0.5, at = 0.82, adj = 1, cex = 1.2, font = 2)
    mtext("n", side = 3, line = 0.5, at = 0.95, adj = 1, cex = 1.2, font = 2)

    # subtitle headers (both columns)
    #text(0.05, 0.92, "", adj=0, cex=1.05, font=2)
    #text(0.82, 0.92, "h",       adj=1, cex=1.05, font=2)
    #text(0.95, 0.92, "n",       adj=1, cex=1.05, font=2)

    if (!nrow(dd)){
      text(0.5,0.5,"(No Data)", cex=1.2, font=2)
      next
    }

    dd <- dd %>% dplyr::arrange(dplyr::desc(n), Element) %>% dplyr::slice_head(n=5)

    y_positions <- seq(0.87, 0.12, length.out = nrow(dd))

    for (i in seq_len(nrow(dd))){
         term_mid <- substr(dd$Element[i], 1, 21)   # mid(term, 1, n)
     text(0.05, y_positions[i],  term_mid , adj=0, cex=1.2, font=2)
     text(0.82, y_positions[i], dd$h[i], adj=1, cex=1.2, font=2)
     text(0.95, y_positions[i], dd$n[i], adj=1, cex=1.2, font=2)

    }
  aac_y <- max(0.01, min(y_positions) - 0.10)
   text(0.8, aac_y, paste("AAC =", round(dd$AAC[1], 2)), cex = 1.1, font = 2, col = "red")
    #text(0.5, 0.001, paste("AAC =", round(dd$AAC[1], 2)), cex=1.1, font=2, col = "red")
  }

  # global title
  ttl <- if (!is.na(total_n) && !is.na(total_h)) {
    paste0("Summary Report (n=", total_n, ", h=", total_h, ")")
  } else if (!is.na(total_n)) {
    paste0("Summary Report (n=", total_n, ")")
  } else {
    "Summary Report"
  }
  mtext(ttl, side=3, outer=TRUE, line=0.3, font=2, cex=1.6)
}



# =========================================================
# Lotka helpers (author publication distribution)
# =========================================================
.build_lotka_input_from_long16 <- function(df_long16) {
  if (is.null(df_long16) || !is.data.frame(df_long16) || !nrow(df_long16)) {
    return(data.frame(papers = integer(), authors = integer(), stringsAsFactors = FALSE))
  }

  fa_col <- intersect(c("FAauthor", "FA", "first author"), names(df_long16))[1]
  ca_col <- intersect(c("CAauthor", "CA", "corresponding author"), names(df_long16))[1]
  if (is.na(fa_col) && is.na(ca_col)) {
    return(data.frame(papers = integer(), authors = integer(), stringsAsFactors = FALSE))
  }

  per_article <- lapply(seq_len(nrow(df_long16)), function(i) {
    vals <- character()
    if (!is.na(fa_col)) vals <- c(vals, as.character(df_long16[[fa_col]][i]))
    if (!is.na(ca_col)) vals <- c(vals, as.character(df_long16[[ca_col]][i]))
    vals <- trimws(vals)
    vals <- vals[!is.na(vals) & nzchar(vals)]
    unique(vals)
  })

  all_authors <- unlist(per_article, use.names = FALSE)
  all_authors <- trimws(as.character(all_authors))
  all_authors <- all_authors[!is.na(all_authors) & nzchar(all_authors)]
  if (!length(all_authors)) {
    return(data.frame(papers = integer(), authors = integer(), stringsAsFactors = FALSE))
  }

  author_pub <- sort(table(all_authors), decreasing = TRUE)
  tb <- as.data.frame(table(as.integer(author_pub)), stringsAsFactors = FALSE)
  names(tb) <- c("papers", "authors")
  tb$papers <- as.integer(as.character(tb$papers))
  tb$authors <- as.integer(tb$authors)
  tb <- tb[order(tb$papers), , drop = FALSE]
  rownames(tb) <- NULL
  tb
}

test_lotka <- function(df,
                       merge_tail   = TRUE,
                       min_expected = 5,
                       make_plot    = TRUE,
                       new_device   = FALSE,
                       device_width = 12,
                       device_height= 5,
                       quiet        = FALSE) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c("papers", "authors") %in% names(df)))

  df <- df[order(df$papers), c("papers", "authors")]
  df <- df[is.finite(df$papers) & is.finite(df$authors), ]
  df <- df[df$papers > 0 & df$authors > 0, ]

  if (nrow(df) < 3) {
    stop("Need at least 3 non-zero categories to test Lotka's law.")
  }

  fit <- lm(log(authors) ~ log(papers), data = df)
  a <- unname(coef(fit)[1])
  b <- unname(coef(fit)[2])
  c_est  <- -b
  A1_est <- exp(a)

  df$expected_raw <- A1_est / (df$papers ^ c_est)
  df$expected <- df$expected_raw / sum(df$expected_raw) * sum(df$authors)
  df$residual <- df$authors - df$expected

  df_test <- df[, c("papers", "authors", "expected")]
  if (merge_tail) {
    while (nrow(df_test) > 3 && any(df_test$expected < min_expected)) {
      i <- nrow(df_test)
      df_test$papers[i - 1]   <- paste0(df_test$papers[i - 1], "+")
      df_test$authors[i - 1]  <- df_test$authors[i - 1] + df_test$authors[i]
      df_test$expected[i - 1] <- df_test$expected[i - 1] + df_test$expected[i]
      df_test <- df_test[-i, ]
    }
  }

  k <- nrow(df_test)
  chisq_stat <- sum((df_test$authors - df_test$expected)^2 / df_test$expected)
  dfree <- k - 3
  p_value <- if (dfree > 0) pchisq(chisq_stat, df = dfree, lower.tail = FALSE) else NA_real_

  if (isTRUE(make_plot)) {
    if (isTRUE(new_device)) {
      try(dev.new(width = device_width, height = device_height), silent = TRUE)
    }
    oldpar <- par(c("mfrow", "mar", "oma", "mgp", "las", "xpd"))
    on.exit(par(oldpar), add = TRUE)
    layout(matrix(c(1, 2), nrow = 1, byrow = TRUE))
    par(mar = c(4.2, 4.2, 2.5, 1.2), oma = c(0, 0, 0.5, 0), mgp = c(2.3, 0.8, 0), las = 1, xpd = NA)

    plot(df$papers, df$authors,
         log  = "xy", pch = 16,
         xlab = "Number of papers (log scale)",
         ylab = "Number of authors (log scale)",
         main = "Lotka log-log plot")
    lines(df$papers, df$expected, lwd = 2, lty = 1)
    legend("topright", legend = c("Observed", "Fitted"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), bty = "n")

    ymax <- max(c(df$authors, df$expected)) * 1.08
    plot(df$papers, df$authors,
         type = "b", pch = 16, ylim = c(0, ymax),
         xlab = "Number of papers", ylab = "Number of authors",
         main = "Observed vs expected")
    lines(df$papers, df$expected, type = "b", pch = 1, lty = 2, lwd = 2)
    legend("topright", legend = c("Observed", "Expected"), pch = c(16, 1), lty = c(1, 2), lwd = c(1, 2), bty = "n")
  }

  if (!isTRUE(quiet)) {
    cat("\n=============================\n")
    cat("Lotka's law test results\n")
    cat("=============================\n")
    cat(sprintf("Estimated exponent c  = %.4f\n", c_est))
    cat(sprintf("Estimated A1          = %.4f\n", A1_est))
    cat(sprintf("R-squared             = %.4f\n", summary(fit)$r.squared))
    cat(sprintf("Chi-square statistic  = %.4f\n", chisq_stat))
    cat(sprintf("Degrees of freedom    = %s\n", ifelse(is.na(dfree), "NA", as.character(dfree))))
    cat(sprintf("P-value               = %s\n", ifelse(is.na(p_value), "NA", format(p_value, digits = 6))))
    if (!is.na(p_value)) {
      if (p_value > 0.05) {
        cat("Conclusion            = Data do not significantly differ from Lotka's law.\n")
      } else {
        cat("Conclusion            = Data significantly differ from Lotka's law.\n")
      }
    } else {
      cat("Conclusion            = Too few grouped categories for chi-square test.\n")
    }
  }

  invisible(list(
    model = fit,
    exponent_c = c_est,
    A1_est = A1_est,
    r_squared = summary(fit)$r.squared,
    chi_square = chisq_stat,
    df = dfree,
    p_value = p_value,
    observed_expected_table = df,
    test_table = df_test
  ))
}

.plot_lotka_result <- function(res) {
  if (is.null(res) || is.null(res$observed_expected_table)) {
    plot.new(); text(0.5, 0.5, "No Lotka result available")
    return(invisible(NULL))
  }
  df <- as.data.frame(res$observed_expected_table, stringsAsFactors = FALSE)
  oldpar <- par(c("mfrow", "mar", "oma", "mgp", "las", "xpd"))
  on.exit(par(oldpar), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1, byrow = TRUE))
  par(mar = c(4.2, 4.2, 2.5, 1.2), oma = c(0, 0, 1.4, 0), mgp = c(2.3, 0.8, 0), las = 1, xpd = NA)

  plot(df$papers, df$authors,
       log = "xy", pch = 16,
       xlab = "Number of papers (log scale)",
       ylab = "Number of authors (log scale)",
       main = "Lotka log-log plot")
  lines(df$papers, df$expected, lwd = 2)
  legend("topright", legend = c("Observed", "Fitted"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), bty = "n")

  ymax <- max(c(df$authors, df$expected), na.rm = TRUE) * 1.08
  plot(df$papers, df$authors,
       type = "b", pch = 16, ylim = c(0, ymax),
       xlab = "Number of papers", ylab = "Number of authors",
       main = "Observed vs expected")
  lines(df$papers, df$expected, type = "b", pch = 1, lty = 2, lwd = 2)
  legend("topright", legend = c("Observed", "Expected"), pch = c(16, 1), lty = c(1, 2), lwd = c(1, 2), bty = "n")

  chi_txt <- paste0(
    "Chi-square = ", format(round(res$chi_square, 4), nsmall = 4),
    "   |   df = ", ifelse(is.na(res$df), "NA", as.character(res$df)),
    "   |   p = ", ifelse(is.na(res$p_value), "NA", format(res$p_value, digits = 6)),
    "   |   c = ", format(round(res$exponent_c, 4), nsmall = 4)
  )
  mtext(chi_txt, outer = TRUE, cex = 0.95, font = 2)
  invisible(NULL)
}

.lotka_summary_table <- function(res) {
  if (is.null(res)) return(NULL)
  conclusion <- if (is.na(res$p_value)) {
    "Too few grouped categories for chi-square test"
  } else if (res$p_value > 0.05) {
    "Data do not significantly differ from Lotka's law"
  } else {
    "Data significantly differ from Lotka's law"
  }
  data.frame(
    Metric = c("Estimated exponent c", "Estimated A1", "R-squared", "Chi-square statistic", "Degrees of freedom", "P-value", "Conclusion"),
    Value = c(
      sprintf("%.4f", res$exponent_c),
      sprintf("%.4f", res$A1_est),
      sprintf("%.4f", res$r_squared),
      sprintf("%.4f", res$chi_square),
      ifelse(is.na(res$df), "NA", as.character(res$df)),
      ifelse(is.na(res$p_value), "NA", format(res$p_value, digits = 6)),
      conclusion
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


# ---------- UI (Unified: WoS base + PubMed-style homepage tabs) ----------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
    .small-note { font-size: 12px; color: #555; }
    .mono { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace; }
    .card { border: 1px solid #e5e5e5; border-radius: 10px; padding: 14px; margin-bottom: 12px; background: #fff; }
    .contact-fab{
      position: fixed;
      right: 18px;
      bottom: 18px;
      z-index: 9999;
      border-radius: 999px;
      padding: 10px 14px;
      background: #2d7ef7;
      color: #fff;
      border: 0;
      box-shadow: 0 6px 18px rgba(0,0,0,0.18);
    }
    .contact-fab:hover{ filter: brightness(0.95); }
    #cmc { background:#e6f4ff; border:2px solid #7bb6ff; }
  ")),
    tags$script(HTML("
      // ---- CMC: remember last value in localStorage (client-side only) ----
      (function(){
        const KEY = 'acsaas_cmc_last';
        function save(){
          const el = document.getElementById('cmc');
          if(!el) return;
          try{ localStorage.setItem(KEY, el.value || ''); }catch(e){}
        }
        function restore(){
          const el = document.getElementById('cmc');
          if(!el) return;
          let v = '';
          try{ v = localStorage.getItem(KEY) || ''; }catch(e){}
          if(v && !el.value){
            el.value = v;
            if(window.Shiny && Shiny.setInputValue){
              Shiny.setInputValue('cmc', v, {priority:'event'});
            }
          }
        }
        document.addEventListener('DOMContentLoaded', function(){
          restore();
          const el = document.getElementById('cmc');
          if(el){
            el.addEventListener('input', save);
            el.addEventListener('change', save);
          }
        });
      })();
    "))
  ),
  titlePanel("WoS stand alone (homepage tabs) — Summary first"),
  fluidRow(
    column(12,
      div(class="card",
        h4("Status / Log"),
        verbatimTextOutput("log_txt", placeholder = TRUE)
      )
    )
  ),
  sidebarLayout(
    sidebarPanel(
      passwordInput("cmc", "CMC / Trial code (optional)", value = "", placeholder = "Type test for once loading data chance, or enter 10-digit CMC"),

      fileInput(
        "wos_xls",
        "Upload WoS file (.xls/.xlsx, max 100MB)",
        accept = c(".xls", ".xlsx")
      ),

      actionButton("btn_run", "Run WoS (Upload)", class="btn-primary"),
      tags$details(
        tags$summary("Advanced (optional)"),
        textInput("base_letter", "Routine 1 base column", value = "AG"),
        tags$p(class="small-note", "If you don't know what this is, keep the default.")
      ),
      tags$hr(),
      textInput("wos_url", "URL (.xls/.xlsx/.csv/.txt)", value = "https://raw.githubusercontent.com/smilechien/woscc/main/LiangPo-Cheng.xls"),
      actionButton("btn_url_run", "Run URL", icon = icon("link")),
      actionButton("btn_demo", "Demo load (data/wosauthor.xls)", class="btn-success"),
      tags$hr(),
      downloadButton("dl_summary_png", "Download report PNG"),
      downloadButton("download_zip", "Download outputs (zip)"),
      tags$hr(),
      tags$p(class="small-note",
             "Order: CMC → Upload → Run WoS → URL → Run URL → Demo. URL run is blocked unless CMC is valid.")
    ),
mainPanel(
      tabsetPanel(id="main_tabs",
        tabPanel("ReadMe", verbatimTextOutput("readme_txt")),
        tabPanel("Summary",
          div(class="card",
              h4("Report (10 domains • Top 5 elements • AAC)"),
              tags$p(class="small-note", "Tip: use the download buttons on the left for PNG/ZIP."),
              plotOutput("summary_png", height = "820px")
          ),
          tags$details(
            tags$summary("Show report table"),
            DTOutput("tbl_summary")
          )
        ),
        tabPanel("AAC", aac_ui("aac")),
        tabPanel("woswide32", DTOutput("tbl_wide")),
        tabPanel("woslong32", DTOutput("tbl_long32")),
        tabPanel("woslong16", DTOutput("tbl_long16")),
        tabPanel("metatable", DTOutput("tbl_meta")),

        
        tabPanel("Network",
          fluidRow(
            column(3,
              wellPanel(
                selectInput("combo_domain", "Choose metadata domain", choices = c("Author","Journal/Year","Year/Articles","Term/Year","Country","Institute","Department","Keyword","State/Province","DocumentType"),
                  selected = "Author"
                ),
                helpText("Pipeline: edge building → FLCA/MS/SIL → Top20 → cluster color → interactive network. For single-valued domains (Year/Journal/DocumentType), isolated nodes are shown.")
              )
            ),
            column(9,
              visNetworkOutput("vn_combo", height = "650px"),
              tags$hr(),
              h4("Nodes (Top20 or isolated)"),
              DTOutput("combo_nodes"),
              tags$hr(),
              h4("Edges (Top20)"),
              DTOutput("combo_edges")
            )
          )
        ),
tabPanel("Chord",
  h4("Chord diagram (same Top-20 nodes/edges as Combo Network)"),
  fluidRow(
    column(4,
      selectInput("chord_domain", "Choose metadata domain (Chord)",
        choices = c("Author","Journal/Year","Year/Articles","Term/Year","Country","Institute","Department","Keyword","State/Province","DocumentType"),
        selected = "Author"
      )
    ),
    column(4,
      checkboxInput("chord_follow_combo", "Sync with Combo Network selection", value = TRUE)
    )
  ),
  tags$p(class="small-note",
         "Uses cluster colors from nodes$carac. If chorddiag is installed, an interactive chord is shown; otherwise a static circlize chord is used."),
  uiOutput("chord_ui"),
  tags$hr(),
  verbatimTextOutput("chord_debug")
),


tabPanel("Author",
          h4("Interactive network (Author, Top20 from WoS)"),
          visNetworkOutput("vn_author", height = "650px"),
          tags$hr(),
          h4("Top20 table"),
          DTOutput("tbl_dom_author"),
          tags$hr(),
          h4("Network nodes / edges"),
          fluidRow(
            column(6, DTOutput("tbl_author_nodes")),
            column(6, DTOutput("tbl_author_edges"))
          )
        ),
        tabPanel("Journal/Year",
          h4("Journal / Year network (from WoS)"),
          visNetworkOutput("vn_journal_year", height = "520px"),
          tags$hr(),
          DTOutput("tbl_dom_journal_year"),
          tags$hr(),
          DTOutput("tbl_journal_year")
        ),
        tabPanel("Year/Articles",
          plotOutput("plot_year_bar", height = "380px"),
          tags$hr(),
          visNetworkOutput("vn_year_articles", height="520px"), tags$hr(), DTOutput("tbl_dom_year_articles")
        ),
        tabPanel("Slope",
          fluidRow(
            column(4, selectInput("slope_domain", "Domain", choices = c("Country","Region","Keyword"), selected = "Country")),
            column(4, sliderInput("slope_recent_years", "Recent years", min = 5, max = 20, value = 10, step = 1)),
            column(4, numericInput("slope_top_n", "Top N", value = 10, min = 3, max = 20, step = 1))
          ),
          plotOutput("plot_year_slope", height = "520px"),
          tags$hr(),
          DTOutput("tbl_dom_year_count")
        ),
        tabPanel("Term/Year",
          plotOutput("plot_term_slope", height = "420px"),
          tags$hr(),
          visNetworkOutput("vn_term_year", height="520px"), tags$hr(), DTOutput("tbl_dom_term_year")
        ),
        tabPanel("Sankey",
         h4("Sankey (choose meta domain)"),
         fluidRow(
           column(4,
                  selectInput("sankey_domain", "Domain (Top-20 nodes)",
                              choices = c("Author","Journal/Year","Year/Articles","Term/Year","Country","State/Province","Institute","Department","Keyword","MeSH"),
                              selected = "Author")
           ),
           column(4,
                  numericInput("sankey_min_weight", "Min weight", value = 1, min = 0, step = 1)
           ),
           column(4,
                  tags$p("If networkD3 is installed, shows an in-app Sankey. Always provides SankeyMATIC text.")
           )
         ),
         uiOutput("sankey_ui"),
         tags$hr(),
         tags$p(tags$strong("SankeyMATIC text (copy/paste)")),
         verbatimTextOutput("sankeymatic_text")
),
tabPanel("SSplot",
                 h4("SSplot (Silhouette panel)"),
                 fluidRow(
                   column(4,
                          selectInput("ss_domain", "Domain (Top-20 nodes)",
                                      choices = c("Author","Journal/Year","Year/Articles","Term/Year","Country","State/Province","Institute","Department","Keyword","MeSH"),
                                      selected = "Author")
                   ),
                   column(4,
                          sliderInput("ss_font_scale", "Font scale (SSplot panel)",
                                      min = 0.6, max = 2.2, value = 1.3, step = 0.1)
                   ),
                   column(4, tags$div(style="padding-top:28px;", ""))
                 ),
                 plotOutput("ssplot_panel", height = "820px")
        ),
        tabPanel("SS kano",
         h4("SS Kano (SS vs a*)"),
         fluidRow(
           column(4,
                  selectInput("ss_kano_domain", "Domain (Top-20 nodes)",
                              choices = c("Author","Journal/Year","Year/Articles","Term/Year","Country","State/Province","Institute","Department","Keyword","MeSH"),
                              selected = "Author")
           ),
           column(4,
                  sliderInput("ss_kano_label_size", "Label size", min = 2, max = 10, value = 4, step = 0.5)
           ),
           column(4,
                  tags$p("Kano: x=a*, y=SS; bubble size=value; color=cluster.")
           )
         ),
         plotOutput("ss_kano_plot", height = "700px")
),
        tabPanel("TAAA",
          h4("TAAA: row-level representative topic"),
          tags$p(class = "small-note",
                 "Built from the WoS metatable keyword rows (A1..A13). Each article row is mapped to the dominant FLCA cluster among its terms in the Keyword domain. If a Profile column exists, Hungarian-mapped agreement tables are also shown."),
          tableOutput("taaa_table"),
          tags$hr(),
          h4("Frequency table of article themes"),
          tableOutput("taaa_theme_freq_table"),
          tags$hr(),
          h4("Profile vs article theme mapping"),
          tableOutput("taaa_profile_map_table"),
          tags$hr(),
          h4("Kappa summary after Hungarian mapping"),
          tableOutput("taaa_kappa_table"),
          tags$hr(),
          h4("K x K table after Hungarian mapping"),
          tableOutput("taaa_confusion_table")
        ),
        tabPanel("Lotka",
          h4("Lotka's law for author publication distribution"),
          tags$p(class = "small-note",
                 "Based on author publication counts from the current WoS run using FA/CA author occurrences per article. Left: log-log Lotka plot. Right: observed vs expected frequencies. The tables report the chi-square goodness-of-fit test."),
          plotOutput("lotka_plot", height = "460px"),
          tags$hr(),
          h4("Lotka chi-square summary"),
          tableOutput("lotka_test_table"),
          tags$hr(),
          h4("Observed vs expected table"),
          tableOutput("lotka_observed_expected_table"),
          tags$hr(),
          h4("Merged tail table used for chi-square test"),
          tableOutput("lotka_chisq_table")
        ),
        tabPanel("Download List", tags$p("Placeholder (not implemented yet).")),
        tabPanel("WebZIP", tags$p("Placeholder (not implemented yet).")),
        tabPanel("IP",
                 h4("IP pass status"),
                 textOutput("ip_pass_status"),
                 uiOutput("ip_access_mode_ui"),
                 verbatimTextOutput("ip_status"),
                 tags$div(class="small-note",
                          "Access is allowed if your IP is in IPlist.txt, or if CMC/test validation passes."),
                 DTOutput("ip_allowlist")
        ),
        tabPanel("Report",
  fluidRow(
    column(12,
      helpText("HTML report includes WoS summary + key domain visuals and PubMed Top20 visuals (if pipeline run)."),
      downloadButton("dl_report_html", "Download HTML report"),
      tags$hr(),
      uiOutput("report_preview")
    )
  )
),
        tabPanel("Most cited",
          fluidRow(
            column(12,
              helpText("Top cited articles from WoS upload. If DOI is available, links use doi.org; otherwise PubMed PMID link if present."),
              DTOutput("tbl_most_cited")
            )
          )
        ),
        
tabPanel("Kano",
         h4("Kano plots (WoS metatable domains)"),
         fluidRow(
           column(4,
                  selectInput("kano_domain", "Domain (Top-20 nodes)",
                              choices = c("Author","Journal/Year","Year/Articles","Term/Year",
                                          "Country","State/Province","Institute","Department","Keyword","MeSH"),
                              selected = "Author")
           ),
           column(4,
                  sliderInput("kano_label_size", "Label size",
                              min = 2, max = 10, value = 4, step = 0.5)
           ),
           column(4,
                  tags$div(style="padding-top:28px; font-weight:600;",
                           textOutput("kano_aac_header"))
           )
         ),
         tags$p(class="small-note",
                "Kano plots are drawn from the selected domain's Top-20 nodes after edge building → FLCA/MS/SIL (or fallback Top20)."),
         h5("Kano: value vs value2"),
         plotOutput("kano_vv2_plot", height = "650px"),
         tags$hr(),
         h5("Kano: SS vs a*"),
         plotOutput("kano_ss_astar_plot", height = "650px")
),

        tabPanel("Map",
                 h4("World map (Country counts; distribution only)"),
                 tags$p(class="small-note",
                        "Built from WoS metatable: FAcountry + CAcountry (pre-FLCA)."),
                 plotOutput("country_map", height = "650px")
        ),
        tabPanel("China",
                 h4("China (Provinces / Municipalities / SARs)"),
                 tags$p(class="small-note",
                        "Built from WoS metatable: FAregion/CAregion where FAcountry/CAcountry indicates China/Taiwan/Hong Kong/Macau."),
                 uiOutput("china_map"),
                 DTOutput("tbl_china_counts")
        ),
        tabPanel("USA",
                 h4("USA (Country)"),
                 tags$p(class="small-note",
                        "Built from WoS metatable: FAcountry + CAcountry where country indicates USA (pre-FLCA)."),
                 plotOutput("usa_country_map", height = "650px"),
                 DTOutput("tbl_usa_counts")
        ),

      )
    )
  )
)

# Floating author contact button (from PubMed homepage)
ui <- tagList(
  ui,
  actionButton("contact_btn", "Contact authors / Request CMC", class = "contact-fab")
)


.profile_col_index <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !ncol(df)) return(NA_integer_)
  nms <- trimws(tolower(as.character(names(df))))
  idx <- which(nms %in% c("profile", "profiles"))
  if (length(idx)) return(idx[1])
  NA_integer_
}

.solve_assignment_max <- function(tab) {
  m <- as.matrix(tab)
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  nr <- nrow(m); nc <- ncol(m)
  if (nr < 1L || nc < 1L) return(integer())
  n <- max(nr, nc)
  mm <- matrix(0, nrow = n, ncol = n)
  mm[seq_len(nr), seq_len(nc)] <- m
  if (requireNamespace("clue", quietly = TRUE)) {
    perm <- clue::solve_LSAP(mm, maximum = TRUE)
    return(as.integer(perm)[seq_len(nr)])
  }
  out <- integer(nr)
  used <- rep(FALSE, n)
  for (i in seq_len(nr)) {
    cand <- which(!used)
    if (!length(cand)) cand <- seq_len(n)
    j <- cand[which.max(mm[i, cand])][1]
    if (!length(j) || is.na(j)) j <- cand[1]
    used[j] <- TRUE
    out[i] <- j
  }
  out
}

.kappa_from_table <- function(tab) {
  m <- as.matrix(tab)
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  n <- sum(m)
  if (!is.finite(n) || n <= 0) return(NA_real_)
  po <- sum(diag(m)) / n
  pe <- sum(rowSums(m) * colSums(m)) / (n * n)
  if (!is.finite(pe) || abs(1 - pe) < .Machine$double.eps) return(NA_real_)
  as.numeric((po - pe) / (1 - pe))
}

.build_taaa_profile_agreement <- function(row_df) {
  if (is.null(row_df) || !is.data.frame(row_df) || !nrow(row_df)) return(NULL)
  if (!all(c("Profile", "Cluster") %in% names(row_df))) return(NULL)
  prof <- suppressWarnings(as.integer(as.character(row_df$Profile)))
  cl <- trimws(as.character(row_df$Cluster))
  keep <- is.finite(prof) & prof > 0L & nzchar(cl)
  if (!any(keep)) return(NULL)
  prof <- prof[keep]
  cl <- cl[keep]
  k <- suppressWarnings(as.integer(max(prof, na.rm = TRUE)))
  if (!is.finite(k) || k < 1L) return(NULL)
  clu_levels <- sort(unique(cl))
  if (length(clu_levels) != k) return(NULL)
  tab0 <- table(factor(prof, levels = seq_len(k)), factor(cl, levels = clu_levels))
  perm <- .solve_assignment_max(tab0)
  if (!length(perm)) return(NULL)
  map_df <- data.frame(
    Profile = seq_len(k),
    Cluster = clu_levels[perm[seq_len(min(k, length(perm)))]],
    stringsAsFactors = FALSE
  )
  map_vec <- stats::setNames(map_df$Profile, map_df$Cluster)
  mapped <- unname(map_vec[cl])
  tab1 <- table(factor(prof, levels = seq_len(k)), factor(mapped, levels = seq_len(k)))
  conf_df <- as.data.frame.matrix(tab1, stringsAsFactors = FALSE)
  conf_df <- cbind(Profile = rownames(conf_df), conf_df, stringsAsFactors = FALSE)
  rownames(conf_df) <- NULL
  kap <- .kappa_from_table(tab1)
  n <- sum(tab1)
  po <- if (n > 0) sum(diag(tab1)) / n else NA_real_
  pe <- if (n > 0) sum(rowSums(tab1) * colSums(tab1)) / (n * n) else NA_real_
  kappa_df <- data.frame(K = k, N = as.integer(n), Agreement = round(po, 4), Expected = round(pe, 4), Kappa = round(kap, 4), stringsAsFactors = FALSE, check.names = FALSE)
  list(mapping_table = map_df, kappa_table = kappa_df, confusion_table = conf_df)
}

.build_taaa_from_metatable <- function(meta_df, nodes_df) {
  if (is.null(meta_df) || !is.data.frame(meta_df) || !nrow(meta_df)) return(NULL)
  if (is.null(nodes_df) || !is.data.frame(nodes_df) || !nrow(nodes_df)) return(NULL)
  if (!all(c("name", "carac") %in% names(nodes_df))) return(NULL)
  term_cols <- grep("^A[0-9]+$", names(meta_df), value = TRUE)
  if (!length(term_cols)) return(NULL)

  nd <- as.data.frame(nodes_df, stringsAsFactors = FALSE)
  nd$name <- trimws(as.character(nd$name))
  nd$carac <- trimws(as.character(nd$carac))
  val_col <- if ("value" %in% names(nd)) "value" else names(nd)[1]
  nd$value <- suppressWarnings(as.numeric(nd[[val_col]]))
  nd$value[!is.finite(nd$value)] <- 0
  nd <- nd[nzchar(nd$name) & nzchar(nd$carac), , drop = FALSE]
  if (!nrow(nd)) return(NULL)

  rep_df <- do.call(rbind, lapply(split(nd, nd$carac), function(z) {
    z <- z[order(-z$value, z$name), , drop = FALSE]
    data.frame(carac = as.character(z$carac[1]), theme_name = as.character(z$name[1]), cluster_n = nrow(z), leader_value = as.numeric(z$value[1]), stringsAsFactors = FALSE)
  }))
  rep_df$cluster_n[!is.finite(rep_df$cluster_n)] <- 0
  rep_df$leader_value[!is.finite(rep_df$leader_value)] <- 0

  term_to_cluster <- stats::setNames(as.character(nd$carac), nd$name)
  pidx <- .profile_col_index(meta_df)
  profile_vec <- if (!is.na(pidx)) meta_df[[pidx]] else rep(NA, nrow(meta_df))

  out <- vector("list", nrow(meta_df))
  kk <- 0L
  for (i in seq_len(nrow(meta_df))) {
    vv <- unlist(meta_df[i, term_cols, drop = TRUE], use.names = FALSE)
    vv <- trimws(as.character(vv[!is.na(vv)]))
    vv <- unique(vv[nzchar(vv)])
    if (!length(vv)) next
    cl <- unname(term_to_cluster[vv])
    cl <- cl[!is.na(cl) & nzchar(cl)]
    if (!length(cl)) next
    tb <- as.data.frame(table(cl), stringsAsFactors = FALSE)
    names(tb) <- c("carac", "count")
    tb$count <- suppressWarnings(as.numeric(tb$count))
    tb$count[!is.finite(tb$count)] <- 0
    tb$theme_name <- rep_df$theme_name[match(tb$carac, rep_df$carac)]
    tb$cluster_n  <- rep_df$cluster_n[match(tb$carac, rep_df$carac)]
    tb$leader_value <- rep_df$leader_value[match(tb$carac, rep_df$carac)]
    tb$theme_name[is.na(tb$theme_name)] <- tb$carac[is.na(tb$theme_name)]
    tb$cluster_n[!is.finite(tb$cluster_n)] <- 0
    tb$leader_value[!is.finite(tb$leader_value)] <- 0
    carac_num <- suppressWarnings(as.numeric(tb$carac))
    carac_num[!is.finite(carac_num)] <- Inf
    tb <- tb[order(-tb$count, carac_num, -tb$leader_value, tb$theme_name), , drop = FALSE]
    chosen <- as.character(tb$carac[1])
    theme <- as.character(tb$theme_name[1])
    theme_count_str <- paste0(tb$theme_name, "(", tb$count, ")", collapse = ", ")
    kk <- kk + 1L
    out[[kk]] <- data.frame(Row = i, Profile = as.character(profile_vec[i] %||% ""), 主題名稱 = theme %||% "", 主題詞頻 = theme_count_str, Cluster = chosen %||% "", stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (kk == 0L) return(NULL)
  row_df <- do.call(rbind, out[seq_len(kk)])
  rownames(row_df) <- NULL
  freq_df <- as.data.frame(table(row_df$主題名稱), stringsAsFactors = FALSE)
  names(freq_df) <- c("主題名稱", "Frequency")
  freq_df$Frequency <- suppressWarnings(as.numeric(freq_df$Frequency))
  freq_df$Frequency[!is.finite(freq_df$Frequency)] <- 0
  freq_df$Cluster <- rep_df$carac[match(freq_df$主題名稱, rep_df$theme_name)]
  freq_df$群大小 <- rep_df$cluster_n[match(freq_df$主題名稱, rep_df$theme_name)]
  freq_df$群首值 <- rep_df$leader_value[match(freq_df$主題名稱, rep_df$theme_name)]
  freq_df$群大小[!is.finite(freq_df$群大小)] <- 0
  freq_df$群首值[!is.finite(freq_df$群首值)] <- 0
  freq_df <- freq_df[order(-freq_df$Frequency, -freq_df$群大小, -freq_df$群首值, freq_df$主題名稱), , drop = FALSE]
  rownames(freq_df) <- NULL
  agree <- .build_taaa_profile_agreement(row_df)
  list(row_table = row_df, freq_table = freq_df, profile_map_table = agree$mapping_table %||% NULL, kappa_table = agree$kappa_table %||% NULL, confusion_table = agree$confusion_table %||% NULL)
}


# ---------- Server ----------
server <- function(input, output, session){

  rv <- reactiveValues(
    log = "",
    woswide32 = NULL,
    woslong32 = NULL,
    woslong16 = NULL,
    metatable = NULL,
    combo_res = NULL,
    chord_res = NULL,
    chord_sig = "",
    summary = NULL,
    total_n = NA_integer_,
    total_h = NA_integer_,
    taaa_df = NULL,
    taaa_freq_df = NULL,
    taaa_profile_map_df = NULL,
    taaa_kappa_df = NULL,
    taaa_conf_df = NULL
  )


  lotka_dist_reactive <- reactive({
    if (is.null(rv$woslong16) || !is.data.frame(rv$woslong16) || !nrow(rv$woslong16)) return(NULL)
    .build_lotka_input_from_long16(rv$woslong16)
  })

  lotka_res_reactive <- reactive({
    df <- lotka_dist_reactive()
    if (is.null(df) || !is.data.frame(df) || nrow(df) < 3) return(NULL)
    tryCatch(test_lotka(df, make_plot = FALSE, quiet = TRUE), error = function(e) NULL)
  })

  # ---- AAC module ----
  aac_server("aac")

  # ---- IP status (ip.txt + IPlist.txt + CMC/Trial) ----
  output$ip_status <- renderPrint(list.files())

  rv$ip_access_type <- "checking"
  rv$ip_addr <- "..."

  observe({
    g <- tryCatch(
      ipm_gate_session(session,
                       cmc = input$cmc %||% "",
                       app_dir = APP_DIR,
                       inc_count_on_allow = FALSE),
      error = function(e) NULL
    )

    if (!is.null(g)) {
      rv$ip_access_type <- g$policy %||% "unknown"
      rv$ip_addr <- g$ip %||% "unknown"
    }
  })

  output$ip_pass_status <- renderText({
    paste0("Access mode: ",
           rv$ip_access_type,
           " | IP: ",
           rv$ip_addr)
  })

  output$ip_access_mode_ui <- renderUI({
    ip <- rv$ip_addr %||% "unknown"
    allow <- tryCatch(ipm_read_iplist(app_dir = APP_DIR), error = function(e) character())
    ipdf <- tryCatch(ipm_read_ip_txt(app_dir = APP_DIR), error = function(e) NULL)

    in_allow <- length(allow) && isTRUE(ip %in% allow)
    in_log <- !is.null(ipdf) && nrow(ipdf) && isTRUE(ip %in% ipdf$ip)

    mode <- if (in_allow) {
      "IP allowlisted (always allowed)"
    } else if (identical(rv$ip_access_type %||% "", "CMC pass")) {
      "CMC validated (allowed)"
    } else if (in_log) {
      "CMC required (trial already used)"
    } else {
      "Trial available (first upload free)"
    }

    col <- if (grepl("always allowed|validated", mode, ignore.case = TRUE)) "#1b7f3a" else if (grepl("Trial", mode, ignore.case = TRUE)) "#b36b00" else "#b00020"

    tags$div(style="margin:6px 0 10px 0;",
             tags$b("Access mode: "),
             tags$span(style=paste0("color:", col, "; font-weight:700;"), mode),
             tags$div(style="color:#666; font-size:12px; margin-top:2px;",
                      "Rule of thumb: allowlist > valid CMC > trial(first time) > blocked/CMC required."))
  })

  output$ip_allowlist <- DT::renderDT({
    allow <- tryCatch(ipm_read_iplist(app_dir = APP_DIR), error = function(e) character())
    ip <- rv$ip_addr %||% "unknown"
    df <- data.frame(
      allowlisted_ip = if (length(allow)) allow else character(),
      is_current_ip = if (length(allow)) (allow == ip) else logical(),
      stringsAsFactors = FALSE
    )
    DT::datatable(df, rownames = FALSE, options = list(pageLength = 25, dom = "tip"))
  })



  # ---- CMC + IP trial gate (ip.txt) ------------------------------------------
# Policy (as requested):
#   - Demo: always allowed (no CMC).
#   - Upload: if this IP has never used Upload before -> allow ONCE, then write ip.txt.
#             if IP already exists in ip.txt -> require 10-digit numeric CMC.
#   - URL autorun via link: cmc=test works ONLY when ip.txt has no record for this IP (one-time trial),
#                           otherwise must provide a valid 10-digit CMC in the URL.
#   - Manual "Run URL" button: 10-digit numeric CMC only (cmc=test is NOT allowed here).

.is_cmc_10 <- function(x){
  x <- as.character(x %||% "")
  x <- trimws(x)
  isTRUE(nzchar(x) && grepl("^[0-9]{10}$", x))
}

# Resolve a stable app directory (avoid getwd() drift under RStudio/runApp)
.app_dir <- local({
  d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (is.null(d) || !nzchar(d) || is.na(d)) d <- getwd()
  tryCatch(normalizePath(d, winslash = "/", mustWork = FALSE), error = function(e) d)
})

.ip_file <- function(){ file.path(.app_dir, "ip.txt") }

.read_ip_txt <- function(){
  f <- .ip_file()
  if (!file.exists(f)) {
    return(data.frame(ip=character(), recent_date=character(), total_count=integer(), cmc=character(),
                      stringsAsFactors = FALSE))
  }
  df <- tryCatch(utils::read.delim(f, sep="	", header=FALSE, stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(df) || ncol(df) < 4) {
    return(data.frame(ip=character(), recent_date=character(), total_count=integer(), cmc=character(),
                      stringsAsFactors = FALSE))
  }
  df <- df[,1:4, drop=FALSE]
  names(df) <- c("ip","recent_date","total_count","cmc")
  df$ip <- as.character(df$ip)
  df$recent_date <- as.character(df$recent_date)
  df$total_count <- suppressWarnings(as.integer(df$total_count))
  df$total_count[!is.finite(df$total_count)] <- 0L
  df$cmc <- as.character(df$cmc)
  df
}

.write_ip_txt <- function(df){
  f <- .ip_file()
  if (!is.data.frame(df) || !nrow(df)) return(invisible(FALSE))
  df <- df[, c("ip","recent_date","total_count","cmc"), drop=FALSE]
  tryCatch(utils::write.table(df, f, sep="	", row.names=FALSE, col.names=FALSE, quote=FALSE),
           error=function(e) NULL)
  invisible(TRUE)
}

.ip_seen <- function(ip){
  df <- .read_ip_txt()
  if (!nrow(df)) return(FALSE)
  ip <- as.character(ip %||% "unknown"); if (!nzchar(ip)) ip <- "unknown"
  any(df$ip == ip)
}

.upsert_ip <- function(ip, cmc, inc_upload_count = FALSE){
  ip <- as.character(ip %||% "unknown"); if (!nzchar(ip)) ip <- "unknown"
  cmc <- as.character(cmc %||% "")
  df <- .read_ip_txt()
  today <- as.character(Sys.Date())
  if (!nrow(df) || !(ip %in% df$ip)) {
    df <- rbind(df, data.frame(ip=ip, recent_date=today, total_count=0L, cmc=cmc, stringsAsFactors = FALSE))
  } else {
    i <- match(ip, df$ip)
    df$recent_date[i] <- today
    df$cmc[i] <- cmc
  }
  if (isTRUE(inc_upload_count)) {
    i <- match(ip, df$ip)
    df$total_count[i] <- as.integer(df$total_count[i] %||% 0L) + 1L
  }
  .write_ip_txt(df)
  df
}

.require_upload_gate <- function(ip){
  # first-time upload: allow once without CMC, but we will record in ip.txt after success
  if (!.ip_seen(ip)) return(TRUE)
  if (.is_cmc_10(input$cmc)) return(TRUE)

  showModal(modalDialog(
    title = "CMC required",
    "Upload has already been used from this IP. Please enter a 10-digit numeric CMC to upload again.",
    easyClose = TRUE,
    footer = modalButton("OK")
  ))
  FALSE
}

.require_manual_url_gate <- function(){
  if (.is_cmc_10(input$cmc)) return(TRUE)
  showModal(modalDialog(
    title = "CMC required",
    "Manual 'Run URL' requires a 10-digit numeric CMC. (cmc=test works only via URL link autorun.)",
    easyClose = TRUE,
    footer = modalButton("OK")
  ))
  FALSE
}

# ---- IP status (persist across UI redraws) ---- (persist across UI redraws) ----
  ip_msg <- reactiveVal(NULL)
  get_client_ip <- function() {
    ip <- tryCatch(session$request$REMOTE_ADDR, error = function(e) NA_character_)
    if (is.null(ip) || !nzchar(ip)) ip <- "unknown"
    ip
  }
  observeEvent(TRUE, {
    ip <- get_client_ip()
    ip_msg(paste0("Client IP: ", ip))
  }, once = TRUE)

  observeEvent(input$btn_ip_check, {
    ip <- get_client_ip()
    allow_raw <- input$ip_allowlist
    allowed <- TRUE
    reason <- "No allowlist set (allowed)."
    if (!is.null(allow_raw) && nzchar(trimws(allow_raw))) {
      allow <- unique(trimws(unlist(strsplit(allow_raw, "\n+"))))
      allow <- allow[nzchar(allow)]
      allowed <- ip %in% allow
      reason <- if (allowed) "Allowed (IP in allowlist)." else "Not allowed (IP not in allowlist)."
    }
    ip_msg(paste0("Client IP: ", ip, " — ", reason))
  })

  output$ip_status_ui <- renderUI({
    msg <- ip_msg()
    if (is.null(msg)) return(tags$div())
    tags$div(style="padding:8px 10px; border-left:4px solid #2b6cb0; background:#f7fafc; margin-bottom:10px;",
             tags$b("IP status: "), msg)
  })
  output$ip_debug <- renderText({
    paste0("REMOTE_ADDR=", tryCatch(session$request$REMOTE_ADDR, error=function(e) "NA"))
  })


  # PubMed cache (raw nodes/edges + computed Top20 with SS/a*)
 # (raw nodes/edges + computed Top20 with SS/a*)
  rv_pm <- reactiveValues(nodes = NULL, edges = NULL, modes = NULL, edges20 = NULL, last_run_msg = "")

  .pm_load <- function(){
    if (!dir.exists(PUBMED_DIR)) return(invisible(FALSE))
    ne <- .load_pubmed_nodes_edges()
    rv_pm$nodes <- ne$nodes
    rv_pm$edges <- ne$edges
    invisible(TRUE)
  }

  # initial load
  .pm_load()

  observeEvent(input$pm_reload, {
    .pm_load()
  })

  observeEvent(input$pm_run_pipeline, {
    if (is.null(rv_pm$nodes) || is.null(rv_pm$edges)) {
      showNotification("PubMed nodes/edges not found. Please reload PubMed data first.", type = "error")
      return(NULL)
    }
    if (!exists("run_flca_ms_sil_runner", mode = "function")) {
      showNotification("PubMed pipeline function run_flca_ms_sil_runner() not found (pubmed/flca_ms_sil_module.R).", type = "error", duration = NULL)
      return(NULL)
    }
    .safe_addlog("[PubMed] Running FLCA + MajorSamplingTop20 + Silhouette ...")
    out <- tryCatch(
      run_flca_ms_sil_runner(rv_pm$nodes, rv_pm$edges, cfg = list(), verbose = FALSE),
      error = function(e) e
    )
    if (inherits(out, "error")) {
      rv_pm$last_run_msg <- conditionMessage(out)
      .safe_addlog(paste0("[PubMed ERROR] ", rv_pm$last_run_msg))
      showNotification(rv_pm$last_run_msg, type = "error", duration = NULL)
      return(NULL)
    }
    rv_pm$modes   <- out$modes
    rv_pm$edges20 <- out$data
    rv_pm$last_run_msg <- sprintf("OK: Top20 nodes=%s edges=%s", nrow(out$modes), nrow(out$data))
    .safe_addlog(paste0("[PubMed] ", rv_pm$last_run_msg))
    showNotification(rv_pm$last_run_msg, type = "message")
  }, ignoreInit = TRUE)

  pm_nodes_current <- reactive({
    if (isTRUE(input$pm_use_top20) && !is.null(rv_pm$modes)) return(rv_pm$modes)
    rv_pm$nodes
  })
  pm_edges_current <- reactive({
    if (isTRUE(input$pm_use_top20) && !is.null(rv_pm$edges20)) return(rv_pm$edges20)
    rv_pm$edges
  })

  # Matched keys between WoS Author Top5 and PubMed nodes$name
  pm_match_tbl <- reactive({
    wos <- rv$summary
    if (is.null(wos) || !is.data.frame(wos)) return(data.frame())
    wos_a <- wos[wos$Domain == "Author", , drop=FALSE]
    if (!nrow(wos_a)) return(data.frame())
    pmn <- pm_nodes_current()
    if (is.null(pmn) || !is.data.frame(pmn) || !"name" %in% names(pmn)) return(data.frame())
    wkey <- .norm_key(wos_a$Element)
    pkey <- .norm_key(pmn$name)
    inter <- intersect(wkey, pkey)
    data.frame(
      key = sort(unique(inter)),
      stringsAsFactors = FALSE
    )
  })

  addlog <- function(x){
    rv$log <- paste0(rv$log, format(Sys.time(), "%H:%M:%S"), " ", x, "\n")
  }

  # safe logger (avoid calling global recursive addlog)
  .safe_addlog <- function(msg){
    tryCatch({ addlog(msg) }, error=function(e) NULL)
  }


  
  # ---- Helpers: build domain edges from WoS metatable and run FLCA + MS Top20 + SS/a* ----
  .trim0 <- function(x){
    x <- as.character(x)
    x[is.na(x)] <- ""
    trimws(x)
  }
  .nz <- function(x){ x <- .trim0(x); x[nzchar(x)] }

  .get_cols <- function(df, cands){
    nms <- names(df); nms_l <- tolower(nms)
    for (c in cands){
      h <- which(nms_l == tolower(c))
      if (length(h)) return(nms[h[1]])
    }
    for (c in cands){
      h <- which(grepl(tolower(c), nms_l, fixed = TRUE))
      if (length(h)) return(nms[h[1]])
    }
    NA_character_
  }

  .term_list_from_Acols <- function(df, a_prefix = "A", k = 13){
    cols <- paste0(a_prefix, seq_len(k))
    cols <- cols[cols %in% names(df)]
    lapply(seq_len(nrow(df)), function(i){
      if (!length(cols)) return(character())
      .nz(unlist(df[i, cols, drop=TRUE]))
    })
  }

  .cooc_edges <- function(term_list){
    # undirected co-occurrence within each article; returns Source/Target/WCD
    acc <- new.env(parent=emptyenv())
    for (items in term_list){
      items <- unique(.nz(items))
      if (length(items) < 2) next
      # all unordered pairs
      for (a_i in seq_len(length(items)-1)){
        for (b_i in (a_i+1):length(items)){
          a <- items[a_i]; b <- items[b_i]
          if (a > b) { tmp <- a; a <- b; b <- tmp }
          key <- paste(a,b,sep="||")
          acc[[key]] <- (acc[[key]] %||% 0L) + 1L
        }
      }
    }
    keys <- ls(acc)
    if (!length(keys)) return(data.frame(Source=character(), Target=character(), WCD=numeric(), stringsAsFactors = FALSE))
    sp <- strsplit(keys, "\\|\\|", fixed = FALSE)
    data.frame(
      Source = vapply(sp, `[[`, "", 1),
      Target = vapply(sp, `[[`, "", 2),
      WCD    = as.numeric(vapply(keys, function(k) acc[[k]], 0L)),
      stringsAsFactors = FALSE
    )
  }

  .bipartite_edges <- function(left_vec, right_vec){
    # edges between paired left/right for each row (article); counts duplicates
    l <- .trim0(left_vec); r <- .trim0(right_vec)
    ok <- nzchar(l) & nzchar(r)
    if (!any(ok)) return(data.frame(Source=character(), Target=character(), WCD=numeric(), stringsAsFactors = FALSE))
    df <- data.frame(L=l[ok], R=r[ok], stringsAsFactors = FALSE)
    agg <- df %>% dplyr::count(L, R, name="WCD")
    data.frame(Source=agg$L, Target=agg$R, WCD=as.numeric(agg$WCD), stringsAsFactors = FALSE)
  }

  .map_terms_to_pubmed <- function(wos_terms, pm_terms){
    # normalize keys and map to PubMed canonical terms
    norm <- function(x){
      x <- tolower(.trim0(x))
      x <- gsub("[^a-z0-9 ]+", " ", x)
      x <- gsub("\\s+", " ", x)
      trimws(x)
    }
    pm_terms <- unique(.nz(pm_terms))
    pm_key <- norm(pm_terms)
    idx <- match(norm(wos_terms), pm_key)
    out <- rep("", length(wos_terms))
    out[!is.na(idx)] <- pm_terms[idx[!is.na(idx)]]
    out
  }

  

  # ---- Normalizers for generic network building ----
  normalize_edges <- function(edges){
    if (is.null(edges) || !is.data.frame(edges) || !nrow(edges)) {
      return(data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE))
    }
    ed <- edges

    # map common column schemes to Leader/Follower/WCD
    nms <- names(ed)

    if (!("Leader" %in% nms)) {
      if ("Source" %in% nms) ed$Leader <- ed$Source
      else if ("from" %in% nms) ed$Leader <- ed$from
      else if ("node1" %in% nms) ed$Leader <- ed$node1
    }
    if (!("Follower" %in% nms)) {
      if ("Target" %in% nms) ed$Follower <- ed$Target
      else if ("to" %in% nms) ed$Follower <- ed$to
      else if ("node2" %in% nms) ed$Follower <- ed$node2
    }
    if (!("WCD" %in% nms)) {
      if ("weight" %in% nms) ed$WCD <- ed$weight
      else if ("value" %in% nms) ed$WCD <- ed$value
      else if ("value2" %in% nms) ed$WCD <- ed$value2
    }

    out <- data.frame(
      Leader   = as.character(ed$Leader %||% ""),
      Follower = as.character(ed$Follower %||% ""),
      WCD      = suppressWarnings(as.numeric(ed$WCD %||% 1)),
      stringsAsFactors = FALSE
    )
    out$Leader <- .trim0(out$Leader)
    out$Follower <- .trim0(out$Follower)
    out$WCD[!is.finite(out$WCD)] <- 1
    out <- out[nzchar(out$Leader) & nzchar(out$Follower), , drop=FALSE]
    out
  }

  normalize_nodes <- function(nodes){
    if (is.null(nodes) || !is.data.frame(nodes) || !nrow(nodes)) {
      return(data.frame(name=character(), stringsAsFactors = FALSE))
    }
    nd <- nodes
    if (!("name" %in% names(nd))) {
      if ("id" %in% names(nd)) nd$name <- nd$id
      else nd$name <- nd[[1]]
    }
    nd$name <- .trim0(as.character(nd$name))
    nd <- nd[nzchar(nd$name), , drop=FALSE]
    nd <- nd[!duplicated(nd$name), , drop=FALSE]
    nd
  }

.run_domain <- function(meta0, domain){
    has_runner <- exists("run_flca_ms_sil_runner", mode="function") || exists("run_flca_ms_sil", mode="function")
    if (!has_runner) {
      # Fallback: we can still build a basic Top20 network from raw edges (no FLCA/MS/SIL)
      # (return ok=TRUE so Network tab is not blank)
      # Note: edges0/nodes0 are built below; we only bypass FLCA later if needed.
    }
    df <- as.data.frame(meta0)

    # Build edges0 + nodes0 (name)
    edges0 <- NULL
    nodes0 <- NULL

    if (domain == "Author"){
      fa <- .get_cols(df, c("FA","FAauthor","first author"))
      ca <- .get_cols(df, c("CA","CAauthor","corresponding author"))
      term_list <- lapply(seq_len(nrow(df)), function(i) unique(.nz(c(df[[fa]][i], df[[ca]][i]))))
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    } else if (domain == "Country"){
      fa <- .get_cols(df, c("FAcountry","fa_country","country_fa"))
      ca <- .get_cols(df, c("CAcountry","ca_country","country_ca"))
      term_list <- lapply(seq_len(nrow(df)), function(i) unique(.nz(c(df[[fa]][i], df[[ca]][i]))))
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    } else if (domain == "Institute"){
      fa <- .get_cols(df, c("FAinstitute","fa_institute"))
      ca <- .get_cols(df, c("CAinstitute","ca_institute"))
      term_list <- lapply(seq_len(nrow(df)), function(i) unique(.nz(c(df[[fa]][i], df[[ca]][i]))))
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    } else if (domain == "Department"){
      fa <- .get_cols(df, c("FADept","FAdept","fa_dept"))
      ca <- .get_cols(df, c("CAdept","CADept","ca_dept"))
      term_list <- lapply(seq_len(nrow(df)), function(i) unique(.nz(c(df[[fa]][i], df[[ca]][i]))))
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    } else if (domain == "State/Province"){
      fa <- .get_cols(df, c("FAregion","fa_region","Region_FA"))
      ca <- .get_cols(df, c("CAregion","ca_region","Region_CA"))
      term_list <- lapply(seq_len(nrow(df)), function(i) unique(.nz(c(df[[fa]][i], df[[ca]][i]))))
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    
    } else if (domain == "Year") {
      cy <- .get_cols(df, c("year","Year","Publication Year","PY"))
      term_list <- lapply(.trim0(df[[cy]]), function(x) .nz(x))
      edges0 <- .cooc_edges(term_list)  # will be empty because each row has single term
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    } else if (domain == "Journal") {
      cj <- .get_cols(df, c("journal","Journal","Source Title","SO"))
      term_list <- lapply(.trim0(df[[cj]]), function(x) .nz(x))
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    } else if (domain == "DocumentType") {
      cd <- .get_cols(df, c("DocumentType","Document Type","DT"))
      term_list <- lapply(.trim0(df[[cd]]), function(x) .nz(x))
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    } else if (domain == "Discipline1") {
      c1 <- .get_cols(df, c("Discipince1","Discipline1","Research Areas","SC"))
      c2 <- .get_cols(df, c("Discipince2","Discipline2"))
      if (!is.null(c2)) {
        term_list <- lapply(seq_len(nrow(df)), function(i) unique(.nz(c(df[[c1]][i], df[[c2]][i]))))
      } else {
        term_list <- lapply(.trim0(df[[c1]]), function(x) .nz(x))
      }
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
} else if (domain == "Keyword"){
      term_list <- .term_list_from_Acols(df, "A", 13)
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    } else if (domain == "MeSH"){
      # Use matched PubMed terms in M1..M13 if present; else fall back to A1..A13
      if (any(grepl("^M\\d+$", names(df)))) {
        term_list <- .term_list_from_Acols(df, "M", 13)
      } else {
        term_list <- .term_list_from_Acols(df, "A", 13)
      }
      edges0 <- .cooc_edges(term_list)
      nodes0 <- data.frame(name=sort(unique(unlist(term_list))), stringsAsFactors = FALSE)
    } else if (domain == "Journal/Year"){
      cj <- .get_cols(df, c("journal","Journal","Source Title","SO"))
      cy <- .get_cols(df, c("year","Year","Publication Year","PY"))
      edges0 <- .bipartite_edges(df[[cj]], df[[cy]])
      nodes0 <- data.frame(name=sort(unique(c(edges0$Source, edges0$Target))), stringsAsFactors = FALSE)
    } else if (domain == "Year/Articles"){
      cy <- .get_cols(df, c("year","Year","Publication Year","PY"))
      ca <- .get_cols(df, c("article_id","PMID","UT","id"))
      edges0 <- .bipartite_edges(df[[cy]], df[[ca]])
      nodes0 <- data.frame(name=sort(unique(c(edges0$Source, edges0$Target))), stringsAsFactors = FALSE)
    } else if (domain == "Term/Year"){
      cy <- .get_cols(df, c("year","Year","Publication Year","PY"))
      term_list <- .term_list_from_Acols(df, "A", 13)
      years <- .trim0(df[[cy]])
      # build term-year edges
      acc <- new.env(parent=emptyenv())
      for (i in seq_len(nrow(df))){
        y <- years[i]
        if (!nzchar(y)) next
        terms <- unique(.nz(term_list[[i]]))
        if (!length(terms)) next
        for (t in terms){
          key <- paste(t, y, sep="||")
          acc[[key]] <- (acc[[key]] %||% 0L) + 1L
        }
      }
      keys <- ls(acc)
      if (!length(keys)) edges0 <- data.frame(Source=character(), Target=character(), WCD=numeric(), stringsAsFactors = FALSE)
      else {
        sp <- strsplit(keys, "\\|\\|", fixed = FALSE)
        edges0 <- data.frame(
          Source = vapply(sp, `[[`, "", 1),
          Target = vapply(sp, `[[`, "", 2),
          WCD    = as.numeric(vapply(keys, function(k) acc[[k]], 0L)),
          stringsAsFactors = FALSE
        )
      }
      nodes0 <- data.frame(name=sort(unique(c(edges0$Source, edges0$Target))), stringsAsFactors = FALSE)
    } else {
      return(list(ok=FALSE, domain=domain, err=paste0("Unknown domain: ", domain)))
    }

    edges0 <- normalize_edges(edges0)
    nodes0 <- normalize_nodes(nodes0)

    # ---- IMPORTANT: compute nodes$value/value2 BEFORE FLCA ----
    # User rule for Author domain:
    #   value  = document count (1st and/or corresponding author counts)
    #   value2 = sum(edge)
    # Do NOT use igraph strength as the primary value for Author.
    if (domain == "Author" && nrow(nodes0)) {
      nm <- trimws(as.character(nodes0$name))
      doc_map <- setNames(numeric(length(nm)), nm)
      if (exists("term_list") && length(term_list)) {
        for (terms in term_list) {
          terms <- unique(trimws(as.character(terms)))
          terms <- terms[nzchar(terms)]
          if (length(terms)) doc_map[terms] <- doc_map[terms] + 1
        }
      }
      edge_sum <- setNames(numeric(length(nm)), nm)
      if (nrow(edges0)) {
        fcol <- intersect(c("Leader","Source","from"), names(edges0))[1]
        tcol <- intersect(c("Follower","Target","to"), names(edges0))[1]
        wcol <- intersect(c("WCD","value","weight"), names(edges0))[1]
        if (!is.na(fcol) && !is.na(tcol) && !is.na(wcol)) {
          for (ii in seq_len(nrow(edges0))) {
            f <- trimws(as.character(edges0[[fcol]][ii])); t <- trimws(as.character(edges0[[tcol]][ii]))
            wv <- suppressWarnings(as.numeric(edges0[[wcol]][ii]))
            if (!is.finite(wv)) wv <- 0
            if (nzchar(f)) edge_sum[f] <- edge_sum[f] + wv
            if (nzchar(t)) edge_sum[t] <- edge_sum[t] + wv
          }
        }
      }
      nodes0$value <- as.numeric(doc_map[nm])
      nodes0$value2 <- as.numeric(edge_sum[nm])
      nodes0$value[!is.finite(nodes0$value) | nodes0$value < 1] <- 1
      nodes0$value2[!is.finite(nodes0$value2)] <- 0
    } else if (nrow(nodes0) && nrow(edges0)) {
      g0 <- tryCatch(
        igraph::graph_from_data_frame(edges0[, c("Leader","Follower","WCD")], directed = FALSE),
        error = function(e) NULL
      )
      if (!is.null(g0) && igraph::vcount(g0) > 0) {
        st <- tryCatch(igraph::strength(g0, weights = igraph::E(g0)$WCD), error = function(e) NULL)
        dg <- tryCatch(igraph::degree(g0), error = function(e) NULL)

        if (!is.null(st) && length(st)) {
          nodes0$value <- as.numeric(st[nodes0$name])
          nodes0$value[!is.finite(nodes0$value)] <- 0
        } else if (!("value" %in% names(nodes0))) {
          nodes0$value <- 1
        }

        if (!is.null(dg) && length(dg)) {
          nodes0$value2 <- as.numeric(dg[nodes0$name])
          nodes0$value2[!is.finite(nodes0$value2)] <- 0
        } else if (!("value2" %in% names(nodes0))) {
          nodes0$value2 <- nodes0$value
        }
      } else {
        if (!("value" %in% names(nodes0))) nodes0$value <- 1
        if (!("value2" %in% names(nodes0))) nodes0$value2 <- nodes0$value
      }
    } else if (nrow(nodes0)) {
      if (!("value" %in% names(nodes0))) nodes0$value <- 1
      if (!("value2" %in% names(nodes0))) nodes0$value2 <- nodes0$value
    }

    if (!nrow(nodes0)) {
        return(list(ok=FALSE, domain=domain, msg="No nodes"))
      }
      # If edges are empty (e.g., single-valued domains like Year/Journal/DocumentType),
      # still return isolated nodes (no FLCA).
      if (!nrow(edges0)) {
        # frequency
        freq <- sort(table(unlist(term_list)), decreasing = TRUE)
        modes <- data.frame(
          name = names(freq),
          value = as.integer(freq),
          value2 = 0L,
          carac = 1L,
          ssi = 0,
          a_star1 = 0,
          stringsAsFactors = FALSE
        )
        return(list(ok=TRUE, domain=domain, modes=modes, edges20=edges0[0,], edges_full=edges0, msg="Isolated nodes (no edges)"))
      }

out <- tryCatch(
  if (exists("run_flca_ms_sil_runner", mode="function")) {
    run_flca_ms_sil_runner(nodes0, edges0, cfg=list(), verbose=FALSE)
  } else {
    run_flca_ms_sil(nodes0, edges0, cfg=list(), verbose=FALSE)
  },
  error=function(e) e
)

if (inherits(out, "error")) {
  # Fallback: basic Top20 by weighted degree, cluster=1
  g0 <- tryCatch(
    igraph::graph_from_data_frame(edges0[, c("Leader","Follower","WCD")], directed=FALSE),
    error=function(e) NULL
  )

  if (!is.null(g0)) {
    strength0 <- tryCatch(
      igraph::strength(g0, weights = igraph::E(g0)$WCD),
      error=function(e) NULL
    )
  } else {
    strength0 <- NULL
  }

  if (!is.null(strength0) && length(strength0)) {
    ord <- order(as.numeric(strength0), decreasing=TRUE)
    top_names <- names(strength0)[ord[seq_len(min(20, length(ord)))]]

    if (domain == "Author") {
      # publication count should come from raw term_list, not nodes0$value
      freq <- sort(table(unlist(term_list)), decreasing = TRUE)
      pub_map <- as.numeric(freq)
      names(pub_map) <- names(freq)

      top_value  <- as.numeric(pub_map[top_names])         # publication count
      top_value2 <- as.numeric(strength0[top_names])       # weighted degree / strength
      top_value[is.na(top_value)] <- 0
      top_value2[is.na(top_value2)] <- 0
    } else {
      top_value  <- as.numeric(strength0[top_names])
      top_value2 <- as.numeric(strength0[top_names])
    }

    modes <- data.frame(
      name = top_names,
      value = top_value,
      value2 = top_value2,
      carac = 1L,
      ssi = 0,
      a_star1 = 0,
      stringsAsFactors = FALSE
    )

    edges20 <- edges0[
      edges0$Leader %in% top_names & edges0$Follower %in% top_names,
      , drop=FALSE
    ]
    edges20 <- edges20[order(edges20$WCD, decreasing=TRUE), , drop=FALSE]
    if (nrow(edges20) > 200) edges20 <- edges20[1:200, , drop=FALSE]

    return(list(
      ok=TRUE,
      domain=domain,
      modes=modes,
      edges20=edges20,
      edges_full=edges0,
      msg=paste0("Fallback (no FLCA): ", conditionMessage(out))
    ))
  }

  return(list(ok=FALSE, domain=domain, err=conditionMessage(out)))
}
    # Standardize columns for vis
    modes <- out$modes
    # Force FLCA output to inherit app-defined Author values:
    #   value  = document count
    #   value2 = sum(edge)
    if (domain == "Author" && is.data.frame(modes) && nrow(modes)) {
      n0 <- trimws(as.character(nodes0$name))
      nm <- trimws(as.character(modes$name))
      ixm <- match(nm, n0)
      if ("value" %in% names(nodes0))  modes$value  <- suppressWarnings(as.numeric(nodes0$value[ixm]))
      if ("value2" %in% names(nodes0)) modes$value2 <- suppressWarnings(as.numeric(nodes0$value2[ixm]))
      modes$value[!is.finite(modes$value) | modes$value < 1] <- 1
      modes$value2[!is.finite(modes$value2)] <- 0
    }
    edges20 <- out$data
    if ("follower" %in% names(edges20) && !("Follower" %in% names(edges20))) names(edges20)[names(edges20)=="follower"] <- "Follower"
    if (all(c("Source","Target") %in% names(edges20)) && !all(c("Leader","Follower") %in% names(edges20))) {
      names(edges20)[match(c("Source","Target"), names(edges20))] <- c("Leader","Follower")
    }
    if (!("WCD" %in% names(edges20))) {
      wcol <- intersect(c("weight","value"), names(edges20))
      if (length(wcol)) { edges20$WCD <- as.numeric(edges20[[wcol[1]]]) } else edges20$WCD <- 1
    }

    list(ok=TRUE, domain=domain, modes=modes, edges20=edges20, edges_full=edges0)
  }

  
  normalize_modes_for_table <- function(modes){
    if (is.null(modes) || !is.data.frame(modes) || !nrow(modes)) return(data.frame())
    m <- modes
    if (!"name" %in% names(m) && "id" %in% names(m)) m$name <- as.character(m$id)
    # standardize columns
    if ("SSi" %in% names(m) && !"ss" %in% names(m)) m$ss <- m$SSi
    if ("SS" %in% names(m) && !"ss" %in% names(m)) m$ss <- m$SS
    if ("a_star1" %in% names(m) && !"a*" %in% names(m)) m[["a*"]] <- m$a_star1
    if ("a_star" %in% names(m) && !"a*" %in% names(m)) m[["a*"]] <- m$a_star
    need <- c("name","value","value2","carac","a*","ss")
    for (k in need) if (!k %in% names(m)) m[[k]] <- NA
    m[, c(need, setdiff(names(m), need)), drop=FALSE]
  }

.normalize_modes_for_table <- normalize_modes_for_table


.as_vis_nodes_edges <- function(res){
    if (is.null(res) || !is.list(res) || !isTRUE(res$ok)) {
      return(list(nodes=data.frame(), edges=data.frame()))
    }
    modes <- res$modes
    edges20 <- res$edges20
    if (is.null(modes) || !is.data.frame(modes) || !nrow(modes)) return(list(nodes=data.frame(), edges=data.frame()))
    if (!("name" %in% names(modes))) {
      if ("id" %in% names(modes)) modes$name <- as.character(modes$id) else modes$name <- as.character(modes[[1]])
    }
    if (!("carac" %in% names(modes))) modes$carac <- "C1"
    # node size: prefer value then value2
    size <- rep(1, nrow(modes))
    if ("value" %in% names(modes)) size <- suppressWarnings(as.numeric(modes$value))
    else if ("value2" %in% names(modes)) size <- suppressWarnings(as.numeric(modes$value2))
    size[!is.finite(size)] <- 1
    nodes <- data.frame(
      id = as.character(modes$name),
      label = as.character(modes$name),
      group = as.character(modes$carac),
      value = size,
      stringsAsFactors = FALSE
    )
    # consistent colors by cluster (group)
    pal <- .cluster_palette_map(nodes$group)
    nodes$color.background <- unname(pal[nodes$group])
    nodes$color.border <- nodes$color.background
    nodes$font.color <- "#000000"
    
    if (is.null(edges20) || !is.data.frame(edges20) || !nrow(edges20)) {
      return(list(nodes=nodes, edges=data.frame()))
    }
    if (!("Leader" %in% names(edges20))) {
      if ("Source" %in% names(edges20)) edges20$Leader <- edges20$Source
    }
    if (!("Follower" %in% names(edges20))) {
      if ("Target" %in% names(edges20)) edges20$Follower <- edges20$Target
    }
    w <- if ("WCD" %in% names(edges20)) suppressWarnings(as.numeric(edges20$WCD)) else 1
    w[!is.finite(w)] <- 1
    edges <- data.frame(from=as.character(edges20$Leader), to=as.character(edges20$Follower), value=w, stringsAsFactors = FALSE)
    list(nodes=nodes, edges=edges, modes=modes, edges_full=res$edges_full)
  }


  output$readme_txt <- renderText({
    p <- "README.md"
    if (file.exists(p)) {
      paste(readLines(p, warn = FALSE), collapse = "
")
    } else {
      paste0(
        "README.md not found.

",
        "How to use:
",
        "1) Upload WoS file (wosauthor.xls/.xlsx)
",
        "2) Click Run
",
        "3) Check: Summary, woswide32, woslong32, woslong16, metatable

",
        "If tabs show nothing, open the Log (Status / Log) panel and check for errors."
      )
    }
  })
  output$log_txt <- renderText({ rv$log })

  # ---- Floating author contact modal (from PubMed homepage) ----
  observeEvent(input$contact_btn, {
    showModal(modalDialog(
      title = "Contact authors / Request CMC",
      easyClose = TRUE,
      footer = modalButton("Close"),
      tags$p("Please choose one of the following ways to contact the authors."),
      tags$h4("1) LINE account"),
      tags$p("LINE Official Account ID:"), tags$code("@onq5657t"),
      tags$p(tags$a("Open LINE add-friend page", href="https://line.me/R/ti/p/%40onq5657t",
                    target="_blank", rel="noopener noreferrer")),
      tags$hr(),
      tags$h4("2) Email"),
      tags$p(tags$a("raschonline.service@gmail.com", href="mailto:raschonline.service@gmail.com")),
      tags$hr(),
      tags$h4("3) ChatGPT group"),
      tags$p("If you need a CMC, please send a short message with your name and purpose.")
    ))
  }, ignoreInit = TRUE)

  
  # ---- unified WoS runner (used by Run + Demo) ----
  run_wos_file <- function(datapath, filename, base_letter = "AG") {

    req(datapath)
    rv$log <- ""
    rv$taaa_df <- NULL
    rv$taaa_freq_df <- NULL
    rv$taaa_profile_map_df <- NULL
    rv$taaa_kappa_df <- NULL
    rv$taaa_conf_df <- NULL
    .safe_addlog(paste0("[RUN] file=", filename, "  base=", base_letter))

    # Routine 1
    .safe_addlog("[R1] wosauthor -> woswide32")
    tryCatch({
      fmls <- names(formals(make_woswide32))
      if ("base_letter" %in% fmls){
        wide <- make_woswide32(datapath, base_letter = base_letter)
      } else {
        wide <- make_woswide32(datapath)
      }
      wide <- as.data.frame(wide)
      wide <- ensure_col(wide, "row_id", seq_len(nrow(wide)))
      wide <- ensure_col(wide, "article_id", wide$row_id)
      rv$woswide32 <- wide
      .safe_addlog(paste0("[R1] rows=", nrow(wide), " cols=", ncol(wide)))
    }, error = function(e){
      .safe_addlog(paste0("[ERROR] ", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      return(NULL)
    })

    req(rv$woswide32)

    # Routine 2/3
    .safe_addlog("[R2/3] woswide32 -> woslong32 (FA/CA metadata)")
    tryCatch({
      fmls <- names(formals(make_woslong32))
      if ("lookup_dir" %in% fmls){
        long32 <- make_woslong32(rv$woswide32, lookup_dir = "data")
      } else {
        long32 <- make_woslong32(rv$woswide32)
      }
      long32 <- as.data.frame(long32)
      long32 <- ensure_col(long32, "row_id", seq_len(nrow(long32)))
      long32 <- ensure_col(long32, "article_id", long32$row_id)
      rv$woslong32 <- long32
      # Normalize CAauthor column name (Caauthor -> CAauthor)
      if ("Caauthor" %in% names(rv$woslong32) && !"CAauthor" %in% names(rv$woslong32)) {
        names(rv$woslong32)[names(rv$woslong32)=="Caauthor"] <- "CAauthor"
      }
      .safe_addlog(paste0("[R2/3] rows=", nrow(long32), " cols=", ncol(long32)))
    }, error = function(e){
      .safe_addlog(paste0("[ERROR] ", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      return(NULL)
    })

    req(rv$woslong32)

    # woslong16 (derived)
    .safe_addlog("[F1] woslong32 -> woslong16 (derived)")
    tryCatch({
      df <- rv$woslong32
      # --- Override FA/CA institute/country/region/dept using make_final_from_woslong_df (addresses-based) ---
      if (exists("make_final_from_woslong_df")) {
        parsed <- tryCatch(make_final_from_woslong_df(df), error = function(e) NULL)
        if (is.data.frame(parsed) && nrow(parsed) == nrow(df)) {
          # map parsed columns back to df (keep existing if parsed missing)
          for (nm in c("FAcountry","CAcountry","FAregion","CAregion","FADept","CAdept","FAinstitute","CAinstitute")) {
            if (nm %in% names(parsed)) df[[nm]] <- parsed[[nm]]
          }
          # Optionally overwrite FAauthor/CAauthor if parsed provides cleaned versions
          if ("FAauthor" %in% names(parsed)) df[["FAauthor"]] <- parsed[["FAauthor"]]
          if ("CAauthor" %in% names(parsed)) df[["CAauthor"]] <- parsed[["CAauthor"]]
        }
      }

      # pick columns robustly
      pick <- function(cands){
        nms <- names(df); nms_l <- tolower(nms)
        c_l <- tolower(cands)
        hit <- which(nms_l %in% c_l)
        if (length(hit)) return(nms[hit[1]])
        for (c in c_l){
          # Avoid accidental matches for very short tokens like 'ca' inside 'publication'
          if (nchar(c) < 3) next
          h <- which(grepl(c, nms_l, fixed = TRUE))
          if (length(h)) return(nms[h[1]])
        }
        NA_character_
      }
      col_citation <- pick(c("citation","Citation"))
      col_year     <- pick(c("publication year","year","Publication.Year"))
      col_journal  <- pick(c("journal","source title","Source Title"))
      col_fa       <- pick(c("FAauthor","FA_author","FA Author","first author","first_author"))
      col_ca       <- pick(c("CAauthor","Caauthor","CA_author","CA Author","corresponding author","corresponding_author"))
      col_doctype  <- pick(c("document type","document.type","Document Type","Document Type"))
      col_disc1    <- pick(c("disci","discipince1","discipline1","research areas","Discipince1"))
      col_disc2    <- pick(c("discipince2","discipline2","wos categories","Discipince2"))

      need <- c("FAcountry","CAcountry","FAinstitute","CAinstitute","FADept","CAdept","FAregion","CAregion")
      for (nm in need) df <- ensure_col(df, nm, "")
      if ("FAcountry" %in% names(df)) df$FAcountry <- .std_country(df$FAcountry)
      if ("CAcountry" %in% names(df)) df$CAcountry <- .std_country(df$CAcountry)
      df <- .apply_region_normalization(as.data.frame(df))

      long16 <- tibble(
        article_id = df$article_id,
        citation = if (!is.na(col_citation)) df[[col_citation]] else "",
        year = if (!is.na(col_year)) df[[col_year]] else "",
        journal = if (!is.na(col_journal)) df[[col_journal]] else "",
        FAauthor = if (!is.na(col_fa)) df[[col_fa]] else "",
        CAauthor = if (!is.na(col_ca)) df[[col_ca]] else "",
        DocumentType = if (!is.na(col_doctype)) df[[col_doctype]] else "",
        Discipince1 = if (!is.na(col_disc1)) df[[col_disc1]] else "",
        Discipince2 = if (!is.na(col_disc2)) df[[col_disc2]] else "",
        FAcountry = df$FAcountry,
        CAcountry = df$CAcountry,
        FAinstitute = df$FAinstitute,
        CAinstitute = df$CAinstitute,
        FADept = df$FADept,
        CAdept = df$CAdept,
        FAregion = df$FAregion,
        CAregion = df$CAregion,
        KeywordsPlus = if ("Keywords.Plus" %in% names(df)) df[["Keywords.Plus"]] else if ("Keywords Plus" %in% names(df)) df[["Keywords Plus"]] else ""
      )
      rv$woslong16 <- as.data.frame(long16)
      # (kept stable like app260226.R) do NOT create FA/CA aliases here
      .safe_addlog(paste0("[F1] rows=", nrow(long16), " cols=", ncol(long16)))
    }, error = function(e){
      .safe_addlog(paste0("[ERROR] ", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      return(NULL)
    })

    req(rv$woslong16)

    # metatable: keywords plus -> A1..A13
    .safe_addlog("[F2] woslong16 -> metatable (A1..A13)")
    tryCatch({
      df16 <- rv$woslong16
      a <- kp_to_Acols(df16$KeywordsPlus, 13)
      meta <- bind_cols(df16 %>% select(-KeywordsPlus), a)
      rv$metatable <- as.data.frame(meta)
      .safe_addlog(paste0("[F2] rows=", nrow(meta), " cols=", ncol(meta)))

      # [F2b] TAAA from Keyword-domain FLCA clusters
      tryCatch({
        kw_res <- .run_domain(rv$metatable, "Keyword")
        kw_nodes <- kw_res$modes %||% kw_res$nodes20 %||% kw_res$nodes %||% NULL
        taaa_res <- .build_taaa_from_metatable(rv$metatable, kw_nodes)
        if (is.list(taaa_res)) {
          rv$taaa_df <- taaa_res$row_table %||% NULL
          rv$taaa_freq_df <- taaa_res$freq_table %||% NULL
          rv$taaa_profile_map_df <- taaa_res$profile_map_table %||% NULL
          rv$taaa_kappa_df <- taaa_res$kappa_table %||% NULL
          rv$taaa_conf_df <- taaa_res$confusion_table %||% NULL
          .safe_addlog(paste0("[F2b] TAAA rows=", nrow(rv$taaa_df %||% data.frame())))
        }
      }, error = function(e) {
        .safe_addlog(paste0("[WARN] TAAA build skipped: ", conditionMessage(e)))
      })
    }, error = function(e){
      .safe_addlog(paste0("[ERROR] ", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      return(NULL)
    })

    req(rv$metatable)

    # summary 10 domains top5
    .safe_addlog("[F3] summary (10 domains x top5; n,h,AAC)")
    tryCatch({
      all_sum <- build_summary10(rv$woswide32, rv$woslong32)
      rv$total_n <- dplyr::n_distinct(rv$woswide32$article_id)
      cc <- detect_cite_col(rv$woswide32)
      if (is.na(cc)) cc <- detect_cite_col(rv$woslong32)
      rv$total_h <- if (!is.na(cc)) h_index_cites(rv$woswide32[[cc]]) else NA_integer_
      rv$summary <- as.data.frame(all_sum)
      .safe_addlog(paste0("[F3] summary rows=", nrow(all_sum)))

      # [F4] Build domain visuals (WoS metatable -> FLCA/MS/SS) to match PubMed standalone tabs
      .safe_addlog("[F4] building domain visuals (Author, Journal/Year, Year/Articles, Term/Year, Country, Institute, Department, Keyword, Region)...")
      tryCatch({
        meta0 <- rv$metatable

        # ---- Map WoS KeywordPlus terms to PubMed terms (for MeSH-like domain) ----
        pm_terms <- character()
        if (!is.null(rv_pm$nodes) && is.data.frame(rv_pm$nodes) && "name" %in% names(rv_pm$nodes)) {
          pm_terms <- rv_pm$nodes$name
        }
        if (length(pm_terms)) {
          # build per-article mapped terms from A1..A13
          a_list <- .term_list_from_Acols(meta0, "A", 13)
          m_list <- lapply(a_list, function(v){
            vv <- .map_terms_to_pubmed(v, pm_terms)
            vv <- unique(.nz(vv))
            vv
          })
          # write M1..M13 columns
          M <- tibble::as_tibble(do.call(rbind, lapply(m_list, function(v){
            v <- v[seq_len(min(13, length(v)))]
            v <- c(v, rep("", 13-length(v)))
            t(v)
          })))
          names(M) <- paste0("M", seq_len(13))
          meta0 <- dplyr::bind_cols(meta0, M)
        }

        rv$dom_author       <- .run_domain(meta0, "Author")
        rv$dom_journal_year <- .run_domain(meta0, "Journal/Year")
        rv$dom_year_articles<- .run_domain(meta0, "Year/Articles")
        rv$dom_term_year    <- .run_domain(meta0, "Term/Year")
        rv$dom_country      <- .run_domain(meta0, "Country")
        rv$dom_institute    <- .run_domain(meta0, "Institute")
        rv$dom_department   <- .run_domain(meta0, "Department")
        rv$dom_keyword      <- .run_domain(meta0, "Keyword")
        rv$dom_mesh         <- .run_domain(meta0, "MeSH")
        rv$dom_region       <- .run_domain(meta0, "State/Province")
        .safe_addlog("[F4] domain visuals built.")
      }, error=function(e){ .safe_addlog(paste0("[ERROR][F4] ", conditionMessage(e))) })

      .safe_addlog("[DONE] All steps completed.")
    }, error = function(e){
      .safe_addlog(paste0("[ERROR] ", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      return(NULL)
    })

  
  }

  # ---- URL parameters: ?cmc=...&autorun=1&csv_url=... ----
  # Supports:
  # - csv_url ending with .xls/.xlsx/.csv -> WoS pipeline (same as Upload+Run)
  # - csv_url ending with .txt -> PubMed route (if pubmed_parse_biblio available; otherwise saves into pubmed/)
  # Robust to the user-provided pattern having extra '?' like "?cmc=test&?autorun=1&csv_url=..."
  .parse_url_params <- function(search){
    if (!is.character(search) || !nzchar(search)) return(list())
    s <- sub("^\\?", "", search)
    parts <- strsplit(s, "&", fixed = TRUE)[[1]]
    kv <- lapply(parts, function(p){
      p <- sub("^\\?", "", p)
      if (!nzchar(p)) return(NULL)
      sp <- strsplit(p, "=", fixed = TRUE)[[1]]
      k <- utils::URLdecode(sp[1])
      v <- if (length(sp) >= 2) utils::URLdecode(paste(sp[-1], collapse="=")) else ""
      k <- sub("^\\?", "", k)
      if (!nzchar(k)) return(NULL)
      setNames(list(v), k)
    })
    out <- do.call(c, kv)
    if (is.null(out)) list() else out
  }

  .download_to_temp <- function(url){
    if (!nzchar(url)) stop("Empty url")
    tf <- tempfile(fileext = paste0(".", tolower(tools::file_ext(url))))
    utils::download.file(url, tf, mode = "wb", quiet = TRUE)
    tf
  }

  .run_wos_csv <- function(csv_path, filename = "wosauthor.csv", base_letter = "AG"){
    rv$log <- ""
    .safe_addlog(paste0("[RUN] file=", filename, " (csv)  base=", base_letter))
    df <- tryCatch(utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) e)
    if (inherits(df, "error") || !is.data.frame(df) || !nrow(df)) {
      msg <- if (inherits(df, "error")) conditionMessage(df) else "CSV is empty"
      .safe_addlog(paste0("[ERROR] CSV read failed: ", msg))
      showNotification(paste0("CSV read failed: ", msg), type="error", duration = NULL)
      return(NULL)
    }
    df <- as.data.frame(df)
    # Heuristic: if it already looks like wide32, accept it as woswide32 and proceed to Routine 2/3
    df <- ensure_col(df, "row_id", seq_len(nrow(df)))
    df <- ensure_col(df, "article_id", df$row_id)
    rv$woswide32 <- df
    .safe_addlog(paste0("[R1] (CSV) rows=", nrow(df), " cols=", ncol(df)))
    # Continue with Routine 2/3+ using the same runner by calling downstream blocks via a tiny wrapper
    # We reuse run_wos_file's downstream sections by temporarily writing to rv and executing the remaining code path.
    # Easiest: call make_woslong32 and onward inline (same as in run_wos_file).
    .safe_addlog("[R2/3] woswide32 -> woslong32 (FA/CA metadata)")
    tryCatch({
      fmls <- names(formals(make_woslong32))
      if ("lookup_dir" %in% fmls){
        long32 <- make_woslong32(rv$woswide32, lookup_dir = "data")
      } else {
        long32 <- make_woslong32(rv$woswide32)
      }
      long32 <- as.data.frame(long32)
      long32 <- ensure_col(long32, "row_id", seq_len(nrow(long32)))
      long32 <- ensure_col(long32, "article_id", long32$row_id)
      rv$woslong32 <- long32
      if ("Caauthor" %in% names(rv$woslong32) && !"CAauthor" %in% names(rv$woslong32)) {
        names(rv$woslong32)[names(rv$woslong32)=="Caauthor"] <- "CAauthor"
      }
      .safe_addlog(paste0("[R2/3] rows=", nrow(long32), " cols=", ncol(long32)))
    }, error = function(e){
      .safe_addlog(paste0("[ERROR] ", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      return(NULL)
    })
    req(rv$woslong32)
    # Reuse the remainder by calling run_wos_file on a fake path is not safe; so we just call the same
    # woslong16/metatable/summary builder by delegating to the existing blocks: simplest is to call
    # run_wos_file on an xls, but we don't have that. Therefore: call the same code blocks by invoking
    # run_wos_file on csv_path only for bookkeeping is not appropriate. We'll just call the same blocks
    # by calling the helper already in memory: we can call the nested code by reusing run_wos_file's
    # logic via sourcing is not possible. So we call run_wos_file's remaining steps by triggering btn_run logic:
    # We'll do a minimal continuation by calling the same [F1..F4] code via a small local function.
    local_continue <- function(){
      # woslong16
      .safe_addlog("[F1] woslong32 -> woslong16 (derived)")
      tryCatch({
        df <- rv$woslong32
        if (exists("make_final_from_woslong_df")) {
          parsed <- tryCatch(make_final_from_woslong_df(df), error = function(e) NULL)
          if (is.data.frame(parsed) && nrow(parsed) == nrow(df)) {
            for (nm in c("FAcountry","CAcountry","FAregion","CAregion","FADept","CAdept","FAinstitute","CAinstitute")) {
              if (nm %in% names(parsed)) df[[nm]] <- parsed[[nm]]
            }
            if ("FAauthor" %in% names(parsed)) df[["FAauthor"]] <- parsed[["FAauthor"]]
            if ("CAauthor" %in% names(parsed)) df[["CAauthor"]] <- parsed[["CAauthor"]]
          }
        }
        pick <- function(cands){
          nms <- names(df); nms_l <- tolower(nms)
          c_l <- tolower(cands)
          hit <- which(nms_l %in% c_l)
          if (length(hit)) return(nms[hit[1]])
          for (c in c_l){
            if (nchar(c) < 3) next
            h <- which(grepl(c, nms_l, fixed = TRUE))
            if (length(h)) return(nms[h[1]])
          }
          NA_character_
        }
        col_citation <- pick(c("citation","Citation"))
        col_year     <- pick(c("publication year","year","Publication.Year","PY"))
        col_journal  <- pick(c("journal","source title","Source Title","SO"))
        col_fa       <- pick(c("FAauthor","FA_author","FA Author","first author","first_author"))
        col_ca       <- pick(c("CAauthor","Caauthor","CA_author","CA Author","corresponding author","corresponding_author"))
        col_doctype  <- pick(c("document type","document.type","Document Type","DT"))
        col_disc1    <- pick(c("disci","discipince1","discipline1","research areas","Discipince1"))
        col_disc2    <- pick(c("discipince2","discipline2","wos categories","Discipince2"))

        need <- c("FAcountry","CAcountry","FAinstitute","CAinstitute","FADept","CAdept","FAregion","CAregion")
        for (nm in need) df <- ensure_col(df, nm, "")
        if ("FAcountry" %in% names(df)) df$FAcountry <- .std_country(df$FAcountry)
        if ("CAcountry" %in% names(df)) df$CAcountry <- .std_country(df$CAcountry)
        df <- .apply_region_normalization(as.data.frame(df))

        long16 <- tibble::tibble(
          article_id = df$article_id,
          citation = if (!is.na(col_citation)) df[[col_citation]] else "",
          year = if (!is.na(col_year)) df[[col_year]] else "",
          journal = if (!is.na(col_journal)) df[[col_journal]] else "",
          FAauthor = if (!is.na(col_fa)) df[[col_fa]] else "",
          CAauthor = if (!is.na(col_ca)) df[[col_ca]] else "",
          DocumentType = if (!is.na(col_doctype)) df[[col_doctype]] else "",
          Discipince1 = if (!is.na(col_disc1)) df[[col_disc1]] else "",
          Discipince2 = if (!is.na(col_disc2)) df[[col_disc2]] else "",
          FAcountry = df$FAcountry,
          CAcountry = df$CAcountry,
          FAinstitute = df$FAinstitute,
          CAinstitute = df$CAinstitute,
          FADept = df$FADept,
          CAdept = df$CAdept,
          FAregion = df$FAregion,
          CAregion = df$CAregion,
          KeywordsPlus = if ("Keywords.Plus" %in% names(df)) df[["Keywords.Plus"]] else if ("Keywords Plus" %in% names(df)) df[["Keywords Plus"]] else ""
        )
        rv$woslong16 <- as.data.frame(long16)
        .safe_addlog(paste0("[F1] rows=", nrow(long16), " cols=", ncol(long16)))
      }, error=function(e){
        .safe_addlog(paste0("[ERROR] ", conditionMessage(e)))
        showNotification(conditionMessage(e), type="error", duration=NULL)
        return(NULL)
      })
      req(rv$woslong16)

      # metatable
      .safe_addlog("[F2] woslong16 -> metatable (A1..A13)")
      tryCatch({
        df16 <- rv$woslong16
        a <- kp_to_Acols(df16$KeywordsPlus, 13)
        meta <- dplyr::bind_cols(df16 %>% dplyr::select(-KeywordsPlus), a)
        meta <- as.data.frame(meta)
        if (!("Profile" %in% names(meta))) meta$Profile <- sample(1:3, nrow(meta), replace = TRUE)
        meta$Profile <- suppressWarnings(as.integer(as.character(meta$Profile)))
        bad_prof <- !is.finite(meta$Profile) | meta$Profile < 1L
        if (any(bad_prof)) meta$Profile[bad_prof] <- sample(1:3, sum(bad_prof), replace = TRUE)
        rv$metatable <- meta
        .safe_addlog(paste0("[F2] rows=", nrow(meta), " cols=", ncol(meta)))
      }, error=function(e){
        .safe_addlog(paste0("[ERROR] ", conditionMessage(e)))
        showNotification(conditionMessage(e), type="error", duration=NULL)
        return(NULL)
      })
      req(rv$metatable)

      # summary + domain visuals
      .safe_addlog("[F3] summary (10 domains x top5; n,h,AAC)")
      tryCatch({
        all_sum <- build_summary10(rv$woswide32, rv$woslong32)
        rv$total_n <- dplyr::n_distinct(rv$woswide32$article_id)
        cc <- detect_cite_col(rv$woswide32)
        if (is.na(cc)) cc <- detect_cite_col(rv$woslong32)
        rv$total_h <- if (!is.na(cc)) h_index_cites(rv$woswide32[[cc]]) else NA_integer_
        rv$summary <- as.data.frame(all_sum)
        .safe_addlog(paste0("[F3] summary rows=", nrow(all_sum)))

        .safe_addlog("[F4] building domain visuals ...")
        tryCatch({
          meta0 <- rv$metatable
          pm_terms <- character()
          if (!is.null(rv_pm$nodes) && is.data.frame(rv_pm$nodes) && "name" %in% names(rv_pm$nodes)) {
            pm_terms <- rv_pm$nodes$name
          }
          if (length(pm_terms)) {
            a_list <- .term_list_from_Acols(meta0, "A", 13)
            m_list <- lapply(a_list, function(v){
              vv <- .map_terms_to_pubmed(v, pm_terms)
              vv <- unique(.nz(vv))
              vv
            })
            M <- tibble::as_tibble(do.call(rbind, lapply(m_list, function(v){
              v <- v[seq_len(min(13, length(v)))]
              v <- c(v, rep("", 13-length(v)))
              t(v)
            })))
            names(M) <- paste0("M", seq_len(13))
            meta0 <- dplyr::bind_cols(meta0, M)
          }
          rv$dom_author       <- .run_domain(meta0, "Author")
          rv$dom_journal_year <- .run_domain(meta0, "Journal/Year")
          rv$dom_year_articles<- .run_domain(meta0, "Year/Articles")
          rv$dom_term_year    <- .run_domain(meta0, "Term/Year")
          rv$dom_country      <- .run_domain(meta0, "Country")
          rv$dom_institute    <- .run_domain(meta0, "Institute")
          rv$dom_department   <- .run_domain(meta0, "Department")
          rv$dom_keyword      <- .run_domain(meta0, "Keyword")
          rv$dom_mesh         <- .run_domain(meta0, "MeSH")
          rv$dom_region       <- .run_domain(meta0, "State/Province")
          .safe_addlog("[F4] domain visuals built.")
        }, error=function(e){ .safe_addlog(paste0("[ERROR][F4] ", conditionMessage(e))) })
        .safe_addlog("[DONE] All steps completed.")
      }, error=function(e){
        .safe_addlog(paste0("[ERROR] ", conditionMessage(e)))
        showNotification(conditionMessage(e), type="error", duration=NULL)
        return(NULL)
      })
    }
    local_continue()
  }

  .handle_pubmed_txt_url <- function(txt_path, url){
    .safe_addlog(paste0("[URL] PubMed txt detected: ", url))
    if (!dir.exists(PUBMED_DIR)) dir.create(PUBMED_DIR, recursive = TRUE, showWarnings = FALSE)
    target_txt <- file.path(PUBMED_DIR, "pubmed_input.txt")
    ok <- tryCatch({ file.copy(txt_path, target_txt, overwrite = TRUE); TRUE }, error=function(e) FALSE)
    # Try parse if parser exists
    parsed <- NULL
    fn <- NULL
    for (cand in c("pubmed_parse_biblio","parse_pubmed_biblio","parse_pubmed_txt")) {
      if (exists(cand, mode="function")) { fn <- get(cand, mode="function"); break }
    }
    if (!is.null(fn)) {
      .safe_addlog("[URL] Running PubMed parser ...")
      fmls <- names(formals(fn))
      parsed <- tryCatch({
        if (all(c("txt_path","out_dir") %in% fmls)) fn(txt_path = target_txt, out_dir = PUBMED_DIR)
        else if (all(c("path","out_dir") %in% fmls)) fn(path = target_txt, out_dir = PUBMED_DIR)
        else if ("txt_path" %in% fmls) fn(txt_path = target_txt)
        else if ("path" %in% fmls) fn(target_txt)
        else fn(target_txt)
      }, error=function(e) e)
      if (inherits(parsed, "error")) {
        .safe_addlog(paste0("[URL][PubMed ERROR] ", conditionMessage(parsed)))
        showNotification(conditionMessage(parsed), type="error", duration=NULL)
      } else {
        .safe_addlog("[URL] PubMed parse completed. Reload PubMed data.")
      }
    } else {
      .safe_addlog("[URL] No PubMed parser found; saved file to pubmed/pubmed_input.txt only.")
    }
    # Reload PubMed nodes/edges if csv files exist
    try(.pm_load(), silent = TRUE)
    try(updateTabsetPanel(session, "main_tabs", selected = "PubMed Visuals"), silent = TRUE)
  }

  # Run once on app load (and also if URL changes inside the same session)
  observeEvent(session$clientData$url_search, {
    params <- .parse_url_params(session$clientData$url_search %||% "")
    if (!length(params)) return(NULL)

    # cmc: sanitize first (strip spaces / trailing punctuation like ',' from copied URLs)
    cmc_val <- params$cmc %||% ""
    cmc_val <- trimws(as.character(cmc_val))
    cmc_val <- sub("[,;]+$", "", cmc_val)

    # set cmc input (passwordInput behaves like textInput for update)
    if (nzchar(cmc_val)) {
      try(updateTextInput(session, "cmc", value = cmc_val), silent = TRUE)
      .safe_addlog("[URL] cmc set via URL param.")
    }

    url <- params$csv_url %||% params$url %||% ""
    if (!nzchar(url)) return(NULL)
    # reflect csv_url into the URL box immediately
    try(updateTextInput(session, "wos_url", value = url), silent = TRUE)

    autorun <- params$autorun %||% "0"
    autorun <- isTRUE(autorun == "1" || tolower(autorun) == "true" || autorun == "yes")

    .safe_addlog(paste0("[URL] csv_url=", url, " autorun=", ifelse(autorun,"1","0")))

    # autorun=0: only fill URL box, do NOT download or run
    if (!isTRUE(autorun)) return(NULL)

    
# autorun=1: enforce link-based gate (cmc=test only for first-time IP)
ip <- get_client_ip()

cmc_val2 <- tolower(trimws(as.character(cmc_val %||% "")))
if (.is_cmc_10(cmc_val2)) {
  # valid CMC from link: allow (and record last CMC)
  .upsert_ip(ip, cmc_val2, inc_upload_count = FALSE)

} else if (cmc_val2 == "test") {
  # test from link: allow ONLY when this IP has never been recorded
  if (.ip_seen(ip)) {
    showModal(modalDialog(
      title = "Trial already used",
      "cmc=test can be used only once (first-time IP) via URL autorun. Please provide a 10-digit CMC in the URL.",
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
    return(NULL)
  }
  .upsert_ip(ip, "test", inc_upload_count = FALSE)

} else {
  showModal(modalDialog(
    title = "CMC invalid",
    "URL autorun requires cmc=test (first-time IP only) or a 10-digit numeric CMC.",
    easyClose = TRUE,
    footer = modalButton("OK")
  ))
  return(NULL)
}

tf <- tryCatch(.download_to_temp(url), error=function(e) e)
    if (inherits(tf, "error")) {
      .safe_addlog(paste0("[URL][ERROR] download failed: ", conditionMessage(tf)))
      showNotification(paste0("Download failed: ", conditionMessage(tf)), type="error", duration=NULL)
      return(NULL)
    }

    ext <- tolower(tools::file_ext(tf))
    # If url has no ext, try from the url
    if (!nzchar(ext)) ext <- tolower(tools::file_ext(url))

    if (ext %in% c("txt")) {
      if (autorun) .handle_pubmed_txt_url(tf, url)
      return(NULL)
    }

    if (ext %in% c("xls","xlsx")) {
      if (autorun) {
run_wos_file(tf, basename(url), input$base_letter)
        try(updateTabsetPanel(session, "main_tabs", selected = "Report"), silent = TRUE)
      }
      return(NULL)
    }

    if (ext %in% c("csv")) {
      if (autorun) {
.run_wos_csv(tf, basename(url), input$base_letter)
        try(updateTabsetPanel(session, "main_tabs", selected = "Report"), silent = TRUE)
      }
      return(NULL)
    }

    .safe_addlog(paste0("[URL] Unsupported file extension: .", ext))
    showNotification(paste0("Unsupported csv_url type: .", ext), type="error", duration=NULL)
  }, ignoreInit = FALSE)


  observeEvent(input$btn_url_run, {
    if (!.require_manual_url_gate()) return(NULL)
    url <- input$wos_url %||% ""
    url <- trimws(url)
    validate(need(nzchar(url), "Please paste a URL to a .xls/.xlsx/.csv/.txt file."))
    .safe_addlog(paste0("[URL] downloading: ", url))
    tf <- tryCatch(.download_to_temp(url), error=function(e) e)
    if (inherits(tf, "error")) {
      .safe_addlog(paste0("[URL][ERROR] download failed: ", conditionMessage(tf)))
      showNotification(paste0("Download failed: ", conditionMessage(tf)), type="error", duration=NULL)
      return(NULL)
    }
    ext <- tolower(tools::file_ext(tf))
    if (!nzchar(ext)) ext <- tolower(tools::file_ext(url))

    if (ext %in% c("xls","xlsx")) {
      withProgress(message = "Processing WoS file", detail = basename(url), value = 0.1, {
        incProgress(0.2, detail = "Reading workbook")
        run_wos_file(tf, basename(url), input$base_letter)
        incProgress(0.7, detail = "Building summary and visuals")
      })
      try(updateTabsetPanel(session, "main_tabs", selected = "Report"), silent = TRUE)
      return(NULL)
    }
    if (ext %in% c("csv")) {
      withProgress(message = "Processing WoS file", detail = basename(url), value = 0.1, {
        incProgress(0.2, detail = "Reading CSV")
        .run_wos_csv(tf, basename(url), input$base_letter)
        incProgress(0.7, detail = "Building summary and visuals")
      })
      try(updateTabsetPanel(session, "main_tabs", selected = "Report"), silent = TRUE)
      return(NULL)
    }
    if (ext %in% c("txt")) {
      .handle_pubmed_txt_url(tf, url)
      return(NULL)
    }

    .safe_addlog(paste0("[URL] Unsupported file extension: .", ext))
    showNotification(paste0("Unsupported URL type: .", ext), type="error", duration=NULL)
  }, ignoreInit = TRUE)

  observeEvent(input$btn_run, {
  ip <- get_client_ip()
  first_free <- !.ip_seen(ip) && !.is_cmc_10(input$cmc)

  if (!.require_upload_gate(ip)) return(NULL)

  req(input$wos_xls)
  withProgress(message = "Processing WoS upload", detail = input$wos_xls$name, value = 0.1, {
    incProgress(0.2, detail = "Reading workbook")
    run_wos_file(input$wos_xls$datapath, input$wos_xls$name, input$base_letter)
    incProgress(0.7, detail = "Building summary and visuals")
  })

  # record upload usage after SUCCESSFUL run
  if (isTRUE(first_free)) {
    .upsert_ip(ip, "free_upload", inc_upload_count = TRUE)
  } else {
    .upsert_ip(ip, as.character(input$cmc %||% ""), inc_upload_count = TRUE)
  }

  try(updateTabsetPanel(session, "main_tabs", selected = "Report"), silent = TRUE)
}, ignoreInit = TRUE)

  observeEvent(input$btn_demo, {
  # Demo is ALWAYS allowed (bypasses CMC and ip.txt gating)
  rv$log <- ""
  demo_xls  <- file.path(.app_dir, "data", "wosauthor.xls")
  demo_xlsx <- file.path(.app_dir, "data", "wosauthor.xlsx")
  datapath <- if (file.exists(demo_xls)) demo_xls else if (file.exists(demo_xlsx)) demo_xlsx else NA_character_

  if (!is.character(datapath) || !nzchar(datapath) || is.na(datapath)) {
    .safe_addlog(paste0("[DEMO][ERROR] Demo file not found. Expected: ", demo_xls, " (or .xlsx)"))
    showNotification("Demo file data/wosauthor.xls(.xlsx) not found", type = "error", duration = NULL)
    return(NULL)
  }

  .safe_addlog(paste0("[DEMO] Loading ", datapath))
  withProgress(message = "Processing demo WoS file", detail = basename(datapath), value = 0.1, {
    incProgress(0.2, detail = "Reading workbook")
    run_wos_file(datapath, basename(datapath), input$base_letter)
    incProgress(0.7, detail = "Building summary and visuals")
  })
  try(updateTabsetPanel(session, "main_tabs", selected = "Report"), silent = TRUE)
}, ignoreInit = TRUE)


  output$tbl_wide <- renderDT({
    req(rv$woswide32)
    datatable(rv$woswide32, options = list(pageLength = 10, scrollX = TRUE))
  })
  output$tbl_long32 <- renderDT({
    req(rv$woslong32)
    datatable(rv$woslong32, options = list(pageLength = 10, scrollX = TRUE))
  })
  output$tbl_long16 <- renderDT({
    req(rv$woslong16)
    datatable(rv$woslong16, options = list(pageLength = 10, scrollX = TRUE))
  })
  output$tbl_meta <- renderDT({
    req(rv$metatable)
    datatable(rv$metatable, options = list(pageLength = 10, scrollX = TRUE))
  })
 
  output$tbl_summary <- renderDT({
    req(rv$summary)
    df <- rv$summary

    dt <- DT::datatable(
      df,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE
    )

    # 只加粗 Domain / Element 欄
    if (all(c("Domain","Element") %in% names(df))) {
      dt <- DT::formatStyle(dt, c("Domain","Element"), fontWeight = "bold")
    }

    dt
  })
  
  # ---- Simple domain visuals (WoS) ----
  output$plot_year_bar <- renderPlot({
  req(rv$metatable)
  y <- suppressWarnings(as.integer(as.character(rv$metatable$year)))
  y <- y[is.finite(y)]
  req(length(y) > 0)

  ymax <- max(y, na.rm = TRUE)
  ymin <- ymax - 9L

  dd <- as.data.frame(table(y), stringsAsFactors = FALSE)
  names(dd) <- c("Year","Count")
  dd$Year <- suppressWarnings(as.integer(as.character(dd$Year)))
  dd$Count <- as.integer(dd$Count)

  dd <- dd[dd$Year >= ymin & dd$Year <= ymax, , drop = FALSE]
  dd <- dd[order(dd$Year), , drop = FALSE]

  ggplot(dd, aes(x = Year, y = Count)) +
    geom_col() +
    scale_x_continuous(breaks = dd$Year) +
    theme_minimal()
})

  .get_year_vector <- function() {
    if (!is.null(rv$metatable) && is.data.frame(rv$metatable) && "Year" %in% names(rv$metatable)) return(suppressWarnings(as.integer(as.character(rv$metatable$Year))))
    if (!is.null(rv$woswide32) && is.data.frame(rv$woswide32)) {
      ycol <- intersect(names(rv$woswide32), c("Publication.Year","Publication Year","Year","PY"))
      if (length(ycol)) return(suppressWarnings(as.integer(as.character(rv$woswide32[[ycol[1]]]))))
    }
    NULL
  }

  .domain_terms_by_row <- function(domain = "Country") {
    mt <- .mt()
    if (is.null(mt) || !nrow(mt)) return(vector("list", 0))
    dom <- tolower(domain)
    cols <- switch(dom,
      "country" = intersect(c("FAcountry","CAcountry","Country","country"), names(mt)),
      "region"  = intersect(c("FAregion","CAregion","Region","region"), names(mt)),
      "keyword" = grep("^A[0-9]+$", names(mt), value = TRUE),
      character(0)
    )
    if (!length(cols)) return(vector("list", nrow(mt)))
    out <- vector("list", nrow(mt))
    for (i in seq_len(nrow(mt))) {
      vv <- unlist(lapply(cols, function(cc) .split_terms_cell(mt[[cc]][i])), use.names = FALSE)
      vv <- trimws(as.character(vv))
      vv <- vv[nzchar(vv)]
      out[[i]] <- unique(vv)
    }
    out
  }

  tufte_sort_local <- function(df, x = "year", y = "value", group = "group", min.space = 0.05) {
    ids <- match(c(x, y, group), names(df))
    df <- df[, ids, drop = FALSE]
    names(df) <- c("x", "y", "group")
    tmp <- expand.grid(x = unique(df$x), group = unique(df$group), stringsAsFactors = FALSE)
    tmp <- merge(df, tmp, all.y = TRUE)
    tmp$y[is.na(tmp$y)] <- 0
    wide <- reshape2::dcast(tmp, group ~ x, value.var = "y")
    ord <- order(wide[, 2])
    wide <- wide[ord, , drop = FALSE]
    min.space <- min.space * diff(range(as.matrix(wide[, -1, drop = FALSE])))
    if (!is.finite(min.space)) min.space <- 1
    yshift <- numeric(nrow(wide))
    if (nrow(wide) >= 2) {
      for (i in 2:nrow(wide)) {
        mat <- as.matrix(wide[(i - 1):i, -1, drop = FALSE])
        d.min <- min(diff(mat))
        yshift[i] <- ifelse(d.min < min.space, min.space - d.min, 0)
      }
    }
    wide$yshift <- cumsum(yshift)
    long <- reshape2::melt(wide, id = c("group", "yshift"), variable.name = "x", value.name = "y")
    long$ypos <- long$y + long$yshift
    long
  }

  plot_slopegraph_local <- function(df, title = "Slopegraph") {
    if (!nrow(df)) return(NULL)
    left_x <- as.character(unique(df$x))[1]
    ylabs <- subset(df, x == left_x)$group
    yvals <- subset(df, x == left_x)$ypos
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = ypos, group = group)) +
      ggplot2::geom_line(color = "red") +
      ggplot2::geom_point(color = "white", size = 6) +
      ggplot2::geom_text(ggplot2::aes(label = y), size = 3) +
      ggplot2::scale_y_continuous(name = "", breaks = yvals, labels = ylabs) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
      ggplot2::labs(title = title)
  }

  .build_slope_counts <- reactive({
    yrs <- .get_year_vector()
    terms_by_row <- .domain_terms_by_row(input$slope_domain %||% "Country")
    if (is.null(yrs) || !length(yrs) || !length(terms_by_row)) return(data.frame())
    n <- min(length(yrs), length(terms_by_row))
    yrs <- yrs[seq_len(n)]
    terms_by_row <- terms_by_row[seq_len(n)]
    keep <- is.finite(yrs)
    yrs <- yrs[keep]
    terms_by_row <- terms_by_row[keep]
    if (!length(yrs)) return(data.frame())
    ymax <- max(yrs, na.rm = TRUE)
    ymin <- ymax - as.integer(input$slope_recent_years %||% 10) + 1L
    rows <- vector("list", 0)
    for (i in seq_along(yrs)) {
      if (!is.finite(yrs[i]) || yrs[i] < ymin || yrs[i] > ymax) next
      tm <- terms_by_row[[i]]
      tm <- tm[nzchar(tm)]
      if (!length(tm)) next
      tm <- unique(tm)
      rows[[length(rows)+1]] <- data.frame(group = tm, year = yrs[i], value = 1, stringsAsFactors = FALSE)
    }
    if (!length(rows)) return(data.frame())
    dd <- dplyr::bind_rows(rows) %>% dplyr::group_by(group, year) %>% dplyr::summarise(value = sum(value), .groups = "drop")
    topn <- as.integer(input$slope_top_n %||% 10)
    tops <- dd %>% dplyr::group_by(group) %>% dplyr::summarise(total = sum(value), .groups = "drop") %>% dplyr::arrange(dplyr::desc(total), group) %>% dplyr::slice_head(n = topn)
    dd <- dd[dd$group %in% tops$group, , drop = FALSE]
    all_years <- seq.int(ymin, ymax)
    dd <- tidyr::complete(dd, group, year = all_years, fill = list(value = 0))
    dd
  })

  output$plot_year_slope <- renderPlot({
    dd <- .build_slope_counts()
    validate(need(is.data.frame(dd) && nrow(dd) > 0, "No slopegraph data yet."))
    df <- tufte_sort_local(dd, x = "year", y = "value", group = "group", min.space = 0.05)
    df$x <- factor(as.character(df$x), levels = as.character(sort(unique(dd$year))))
    print(plot_slopegraph_local(df, title = paste0(input$slope_domain %||% "Country", " slopegraph")))
  })

  output$tbl_dom_year_count <- renderDT({
    dd <- .build_slope_counts()
    validate(need(is.data.frame(dd) && nrow(dd) > 0, "No year-count data."))
    wide <- reshape2::dcast(dd, group ~ year, value.var = "value", fill = 0)
    wide$total <- rowSums(wide[, -1, drop = FALSE])
    wide <- wide[order(-wide$total, wide$group), , drop = FALSE]
    DT::datatable(wide, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$plot_term_slope <- renderPlot({
    dd <- .build_slope_counts()
    validate(need(is.data.frame(dd) && nrow(dd) > 0, "No slopegraph data yet."))
    df <- tufte_sort_local(dd, x = "year", y = "value", group = "group", min.space = 0.05)
    df$x <- factor(as.character(df$x), levels = as.character(sort(unique(dd$year))))
    print(plot_slopegraph_local(df, title = paste0(input$slope_domain %||% "Country", " slopegraph")))
  })

  # ---- PubMed Top20 network plot ----
  output$pm_plot_network <- renderPlot({
    nodes <- pm_nodes_current()
    edges <- pm_edges_current()
    req(nodes, edges)
    # infer columns
    nms <- names(nodes)
    name_col <- if ("name" %in% nms) "name" else if ("Element" %in% nms) "Element" else nms[1]
    v <- as.character(nodes[[name_col]])
    v <- v[!is.na(v) & nzchar(v)]
    # edges: try Leader/Follower else from/to
    ecn <- names(edges)
    from_col <- intersect(ecn, c("Leader","from","node1","Source","source","V1"))[1]
    to_col   <- intersect(ecn, c("Follower","to","node2","Target","target","V2"))[1]
    validate(need(!is.na(from_col) && !is.na(to_col), "Edges need Leader/Follower (or from/to) columns."))
    ed <- edges[, c(from_col, to_col)]
    names(ed) <- c("from","to")
    ed$from <- as.character(ed$from); ed$to <- as.character(ed$to)
    ed <- ed[ed$from %in% v & ed$to %in% v, , drop=FALSE]
    g <- igraph::graph_from_data_frame(ed, directed = TRUE, vertices = data.frame(name=v))
    # size by value2 if available
    size_col <- intersect(names(nodes), c("value2","Value2","strength","size","n"))[1]
    vsz <- rep(6, igraph::vcount(g))
    if (!is.na(size_col)) {
      vv <- suppressWarnings(as.numeric(nodes[[size_col]]))
      vv <- vv[match(igraph::V(g)$name, as.character(nodes[[name_col]]))]
      vv[is.na(vv)] <- median(vv, na.rm=TRUE)
      vv <- sqrt(pmax(vv, 0))
      vsz <- 4 + 10*(vv/max(vv))
    }
    set.seed(1)
    lay <- igraph::layout_with_fr(g)
    plot(g, layout=lay, vertex.size=vsz, vertex.label.cex=0.7, vertex.label.dist=0.4,
         edge.arrow.size=0.3)
  })

  # ---- HTML report (download + preview) ----
  make_report_html <- reactive({
    req(rv$summary)
    # helper to capture plot to base64
    png_b64 <- function(plot_fun, width=2000, height=1500){
      tf <- tempfile(fileext = ".png")
      grDevices::png(tf, width=width, height=height, res=200)
      on.exit({ try(grDevices::dev.off(), silent=TRUE) }, add=TRUE)
      plot_fun()
      grDevices::dev.off()
      raw <- readBin(tf, "raw", n = file.info(tf)$size)
      base64enc::base64encode(raw)
    }
    # Summary (Top5 per domain) as HTML blocks (avoid PNG overlap)
    sum_tbl <- rv$summary
    sum_tbl <- sum_tbl[, intersect(names(sum_tbl), c("Domain","Element","n","h","AAC")), drop=FALSE]

    # Ensure numeric
    if ("n" %in% names(sum_tbl)) sum_tbl$n <- suppressWarnings(as.numeric(sum_tbl$n))
    if ("h" %in% names(sum_tbl)) sum_tbl$h <- suppressWarnings(as.numeric(sum_tbl$h))
    if ("AAC" %in% names(sum_tbl)) sum_tbl$AAC <- suppressWarnings(as.numeric(sum_tbl$AAC))

    # Domain order (keep your original order if present)
    dom_order <- c("Country","Institute","Department","Author","Journal","Year","Document Type","Discipline","Keyword","Region",
                   "State/Province","MeSH","Term/Year","Year/Articles","Journal/Year")
    doms <- unique(as.character(sum_tbl$Domain))
    doms <- c(intersect(dom_order, doms), setdiff(doms, dom_order))

    domain_block <- function(dom){
      dd <- sum_tbl[sum_tbl$Domain == dom, , drop=FALSE]
      if (!nrow(dd)) return(NULL)
      # keep top5 if more
      if ("n" %in% names(dd)) dd <- dd[order(-dd$n, dd$Element), , drop=FALSE]
      dd <- utils::head(dd, 5)
      aac <- if ("AAC" %in% names(dd) && any(is.finite(dd$AAC))) dd$AAC[which(is.finite(dd$AAC))[1]] else NA_real_

      htmltools::tags$div(
        style="border:1px solid #e6e6e6; border-radius:10px; padding:10px; background:#fff;",
        htmltools::tags$div(
          style="display:flex; align-items:flex-end; justify-content:space-between; margin-bottom:6px;",
          htmltools::tags$div(dom, style="color:#d32f2f; font-weight:700; font-size:15px;"),
          htmltools::tags$div(
            style="font-size:12px; color:#444; font-weight:600;",
            htmltools::HTML("h&nbsp;&nbsp;&nbsp;&nbsp;n")
          )
        ),
        htmltools::tags$table(
          style="border-collapse:collapse; width:100%; font-size:12px;",
          htmltools::tags$tbody(
            lapply(seq_len(nrow(dd)), function(i){
              htmltools::tags$tr(
                htmltools::tags$td(htmltools::htmlEscape(as.character(dd$Element[i])),
                                   style="padding:2px 4px; text-align:left; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; max-width:280px;"),
                htmltools::tags$td(ifelse(is.finite(dd$h[i]), sprintf("%d", as.integer(dd$h[i])), ""),
                                   style="padding:2px 4px; text-align:right; width:44px;"),
                htmltools::tags$td(ifelse(is.finite(dd$n[i]), sprintf("%d", as.integer(dd$n[i])), ""),
                                   style="padding:2px 4px; text-align:right; width:44px;")
              )
            })
          )
        ),
        htmltools::tags$div(
          style="margin-top:6px; font-size:12px; font-weight:700; color:#111;",
          sprintf("AAC = %s", ifelse(is.finite(aac), format(round(aac, 2), nsmall=2), "NA"))
        )
      )
    }

    summary_blocks <- htmltools::tags$div(
      style="display:grid; grid-template-columns: 1fr 1fr; gap:14px;",
      lapply(doms, domain_block)
    )

    # Year plots (build ggplot directly)
    # Year plots (build ggplot directly)
    year_bar_fun <- function(){
      df <- rv$woswide32; ycol <- intersect(names(df), c("Publication.Year","Publication Year","Year","PY"))[1]
      y <- suppressWarnings(as.integer(as.character(df[[ycol]]))); y <- y[!is.na(y)]
      dd <- as.data.frame(sort(table(y), decreasing = FALSE)); names(dd) <- c("Year","Count")
      print(ggplot(dd, aes(x=Year, y=Count)) + geom_col() + theme_minimal())
    }
    year_slope_fun <- function(){
      df <- rv$woswide32; ycol <- intersect(names(df), c("Publication.Year","Publication Year","Year","PY"))[1]
      y <- suppressWarnings(as.integer(as.character(df[[ycol]]))); y <- y[!is.na(y)]
      dd <- as.data.frame(sort(table(y), decreasing = FALSE)); names(dd) <- c("Year","Count")
      print(ggplot(dd, aes(x=as.integer(as.character(Year)), y=as.integer(Count))) +
              geom_line() + geom_point() + theme_minimal() + labs(x="Year", y="Count"))
    }
    term_slope_fun <- function(){
      kw <- rv$summary[rv$summary$Domain == "Keyword", , drop=FALSE]
      if (!nrow(kw)) { plot.new(); text(0.5,0.5,"No Keyword data"); return() }
      kw <- kw[order(-kw$n), , drop=FALSE]; kw$Rank <- seq_len(nrow(kw))
      print(ggplot(kw, aes(x=Rank, y=n)) + geom_line() + geom_point() + theme_minimal() +
              scale_x_continuous(breaks = kw$Rank, labels = kw$Element) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(x="Keyword (ranked)", y="Count"))
    }
    b64_year_bar   <- png_b64(year_bar_fun, 1200, 800)
    b64_year_slope <- png_b64(year_slope_fun, 1200, 800)
    b64_term_slope <- png_b64(term_slope_fun, 1200, 800)

    
    pm_net_fun <- function(){
      nodes <- pm_nodes_current()
      edges <- pm_edges_current()
      if (is.null(nodes) || is.null(edges)) { plot.new(); text(0.5,0.5,"No PubMed data"); return() }
      nms <- names(nodes)
      name_col <- if ("name" %in% nms) "name" else if ("Element" %in% nms) "Element" else nms[1]
      v <- as.character(nodes[[name_col]])
      v <- v[!is.na(v) & nzchar(v)]
      ecn <- names(edges)
      from_col <- intersect(ecn, c("Leader","from","node1","Source","source","V1"))[1]
      to_col   <- intersect(ecn, c("Follower","to","node2","Target","target","V2"))[1]
      if (is.na(from_col) || is.na(to_col)) { plot.new(); text(0.5,0.5,"Edges need Leader/Follower"); return() }
      ed <- edges[, c(from_col, to_col)]
      names(ed) <- c("from","to")
      ed$from <- as.character(ed$from); ed$to <- as.character(ed$to)
      ed <- ed[ed$from %in% v & ed$to %in% v, , drop=FALSE]
      g <- igraph::graph_from_data_frame(ed, directed = TRUE, vertices = data.frame(name=v))
      size_col <- intersect(names(nodes), c("value2","Value2","strength","size","n"))[1]
      vsz <- rep(6, igraph::vcount(g))
      if (!is.na(size_col)) {
        vv <- suppressWarnings(as.numeric(nodes[[size_col]]))
        vv <- vv[match(igraph::V(g)$name, as.character(nodes[[name_col]]))]
        vv[is.na(vv)] <- median(vv, na.rm=TRUE)
        vv <- sqrt(pmax(vv, 0))
        if (max(vv, na.rm=TRUE) > 0) vsz <- 4 + 10*(vv/max(vv, na.rm=TRUE))
      }
      set.seed(1)
      lay <- igraph::layout_with_fr(g)
      plot(g, layout=lay, vertex.size=vsz, vertex.label.cex=0.7, vertex.label.dist=0.4,
           edge.arrow.size=0.3)
    }

# PubMed plots if available
    b64_pm_net <- NULL
    if (!is.null(rv_pm$modes) && !is.null(rv_pm$edges20)) {
      b64_pm_net <- png_b64(pm_net_fun, 1400, 1000)
    }

    # Simple table html for summary
    sum_tbl <- rv$summary
    sum_tbl <- sum_tbl[, intersect(names(sum_tbl), c("Domain","Element","n","h","AAC")), drop=FALSE]
    tbl_html <- htmltools::tags$table(
      style="border-collapse:collapse; width:100%;",
      htmltools::tags$thead(
        htmltools::tags$tr(lapply(names(sum_tbl), function(h) htmltools::tags$th(h, style="border:1px solid #ccc; padding:6px; text-align:left;")))
      ),
      htmltools::tags$tbody(
        lapply(seq_len(nrow(sum_tbl)), function(i){
          htmltools::tags$tr(lapply(sum_tbl[i,], function(x) htmltools::tags$td(as.character(x), style="border:1px solid #ddd; padding:6px;")))
        })
      )
    )

    # Build HTML
    body <- htmltools::tagList(
      htmltools::tags$h2("WoS Standalone Report"),
      htmltools::tags$p(sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
      htmltools::tags$h3("WoS Summary (10 domains, Top5)"),
      summary_blocks,
      htmltools::tags$h3("Summary Table"),
      tbl_html,
      htmltools::tags$h3("Year/Articles"),
      htmltools::HTML(sprintf('<img style="max-width:100%%; height:auto;" src="data:image/png;base64,%s"/>', b64_year_bar)),
      htmltools::tags$h3("Slope (Year trend)"),
      htmltools::HTML(sprintf('<img style="max-width:100%%; height:auto;" src="data:image/png;base64,%s"/>', b64_year_slope)),
      htmltools::tags$h3("Term/Year (Keyword proxy, will bridge later)"),
      htmltools::HTML(sprintf('<img style="max-width:100%%; height:auto;" src="data:image/png;base64,%s"/>', b64_term_slope))
    )
    if (!is.null(b64_pm_net)) {
      body <- htmltools::tagAppendChildren(body,
        htmltools::tags$h3("PubMed Top20 Network"),
        htmltools::HTML(sprintf('<img style="max-width:100%%; height:auto;" src="data:image/png;base64,%s"/>', b64_pm_net))
      )
    } else {
      body <- htmltools::tagAppendChildren(body,
        htmltools::tags$h3("PubMed Top20 Network"),
        htmltools::tags$p("Not available yet. Run PubMed pipeline first.")
      )
    }
    htmltools::renderTags(htmltools::tags$html(htmltools::tags$head(
      htmltools::tags$meta(charset="utf-8"),
      htmltools::tags$title("Unified Report")
    ), htmltools::tags$body(style="font-family: Arial, sans-serif; margin: 16px;", body)))$html
  })

  output$report_preview <- renderUI({
    req(rv$summary)
    HTML(make_report_html())
  })

  output$dl_report_html <- downloadHandler(
    filename = function(){ sprintf("unified_report_%s.html", format(Sys.Date(), "%Y%m%d")) },
    content = function(file){
      writeLines(make_report_html(), con=file, useBytes=TRUE)
    }
  )

output$tbl_most_cited <- renderDT({
    req(rv$woswide32)
    w <- rv$woswide32
    cite_col <- detect_cite_col(w)
    if (is.na(cite_col)) {
      return(datatable(data.frame(Note = "No citation column found in WoS upload."), options = list(dom='t')))
    }

    # try common metadata columns
    pick_col <- function(cands){
      nms <- names(w); nms_l <- tolower(nms)
      c_l <- tolower(cands)
      hit <- which(nms_l %in% c_l)
      if (length(hit)) return(nms[hit[1]])
      for (c in c_l) {
        h <- which(grepl(c, nms_l, fixed = TRUE))
        if (length(h)) return(nms[h[1]])
      }
      NA_character_
    }

    col_title <- pick_col(c("ti","title"))
    col_year  <- pick_col(c("publication year","year"))
    col_jour  <- pick_col(c("journal","source title"))
    col_doi   <- pick_col(c("doi","di","doc","document object identifier"))
    col_pmid  <- pick_col(c("pmid"))
    col_cite  <- cite_col

    df <- w
    df[[col_cite]] <- suppressWarnings(as.numeric(df[[col_cite]]))

    # DOI fallback: extract from Citation text if no DOI column
    doi_from_text <- function(x){
      x <- as.character(x)
      m <- regexpr("10\\.[0-9]{4,9}/[^\\s;]+", x, perl = TRUE)
      out <- ifelse(m > 0, regmatches(x, m), NA_character_)
      out
    }
    if (is.na(col_doi)) {
      # try citation column names too
      col_citation_txt <- pick_col(c("citation","citation_raw"))
      if (!is.na(col_citation_txt)) {
        df$DOI_extracted <- doi_from_text(df[[col_citation_txt]])
        col_doi <- "DOI_extracted"
      }
    }

    # build link: prefer DOI, else PMID
    link_id <- rep(NA_character_, nrow(df))
    if (!is.na(col_doi) && col_doi %in% names(df)) link_id <- df[[col_doi]]
    link_id <- trimws(as.character(link_id))
    link_type <- ifelse(!is.na(link_id) & nzchar(link_id), "doi", "")
    if (!is.na(col_pmid) && col_pmid %in% names(df)) {
      pm <- trimws(as.character(df[[col_pmid]]))
      use_pmid <- (link_type == "") & !is.na(pm) & nzchar(pm)
      link_id[use_pmid] <- pm[use_pmid]
      link_type[use_pmid] <- "pmid"
    }

    make_href <- function(id, type){
      if (is.na(id) || !nzchar(id)) return("")
      if (type == "doi") return(paste0("https://doi.org/", id))
      if (type == "pmid") return(paste0("https://pubmed.ncbi.nlm.nih.gov/", id, "/"))
      ""
    }
    href <- mapply(make_href, link_id, link_type, USE.NAMES = FALSE)

    title_txt <- if (!is.na(col_title)) as.character(df[[col_title]]) else rep("", nrow(df))
    title_txt <- ifelse(is.na(title_txt) | !nzchar(title_txt), "(no title)", title_txt)
    title_link <- ifelse(nzchar(href), sprintf('<a href="%s" target="_blank" rel="noopener noreferrer">%s</a>', href, htmlEscape(title_txt)), htmlEscape(title_txt))

    out <- data.frame(
      Cited = df[[col_cite]],
      Year = if (!is.na(col_year)) df[[col_year]] else NA,
      Journal = if (!is.na(col_jour)) df[[col_jour]] else NA,
      Title = title_link,
      LinkType = ifelse(link_type=="", NA, link_type),
      LinkID = ifelse(link_type=="", NA, link_id),
      stringsAsFactors = FALSE
    )

    out <- out[order(out$Cited, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    out <- utils::head(out, 50)

    datatable(out, escape = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })
  output$summary_png <- renderPlot({
    req(rv$summary)
    plot_summary10(rv$summary, total_n = rv$total_n, total_h = rv$total_h)
  })

  
  output$dl_summary_png <- downloadHandler(
    filename = function(){ paste0("summary_top5_aac_", format(Sys.Date(), "%Y%m%d"), ".png") },
    content = function(file){
      req(rv$summary)
      png(file, width = 1400, height = 2300, res = 300)
      plot_summary10(rv$summary, total_n = rv$total_n, total_h = rv$total_h)
      dev.off()
    }
  )

output$download_zip <- downloadHandler(
    contentType = "application/zip",
    filename = function(){
      paste0("wos_allinone_outputs_", format(Sys.Date(), "%Y%m%d"), ".zip")
    },
    content = function(file){
      req(rv$woswide32, rv$woslong32, rv$woslong16, rv$metatable, rv$summary)

      tmpdir <- tempfile("wos_out_")
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

      write.csv(rv$woswide32, file.path(tmpdir, "woswide32.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      write.csv(rv$woslong32, file.path(tmpdir, "woslong32.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      write.csv(rv$woslong16, file.path(tmpdir, "woslong16.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      write.csv(rv$metatable, file.path(tmpdir, "wosmetatable.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      write.csv(rv$summary, file.path(tmpdir, "summary_report.csv"), row.names = FALSE, fileEncoding = "UTF-8")

        # save png (make it wider to avoid overlap)
      png_path <- file.path(tmpdir, "summary_report.png")

      w <- session$clientData$output_summary_png_width
      h <- session$clientData$output_summary_png_height
      if (is.null(w) || is.na(w)) w <- 1400
      if (is.null(h) || is.na(h)) h <- 900

      scale_w <- 2.0   # ← 1.6~2.2 自己調；2.0 通常就很明顯改善
      scale_h <- 1.9   # ← 可 1.0~1.4；不想變高就設 1.0

      png(png_path,
          width  = as.integer(w * scale_w),
          height = as.integer(h * scale_h),
          res = 300)

      plot_summary10(rv$summary, total_n = rv$total_n, total_h = rv$total_h)
      dev.off()

      # include app + scripts + data lookup
      dir.create(file.path(tmpdir, "R"), showWarnings = FALSE)
      file.copy(file.path("R","wos_vba_port.R"), file.path(tmpdir, "R", "wos_vba_port.R"), overwrite = TRUE)
      file.copy("app.R", file.path(tmpdir, "app.R"), overwrite = TRUE)
      if (file.exists("README.md")) file.copy("README.md", file.path(tmpdir, "README.md"), overwrite = TRUE)
      if (dir.exists("data")){
        dir.create(file.path(tmpdir, "data"), showWarnings = FALSE)
        file.copy(list.files("data", full.names = TRUE), file.path(tmpdir, "data"), overwrite = TRUE, recursive = TRUE)
      }

      oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE)
      setwd(tmpdir)
      zip::zipr(file, list.files(tmpdir, recursive = TRUE))
    }
  )
  output$tbl_pm_nodes <- renderDT({
    pmn <- pm_nodes_current()
    if (is.null(pmn) || !is.data.frame(pmn)) return(datatable(data.frame()))
    if (isTRUE(input$pm_show_matched_only)) {
      mk <- pm_match_tbl()
      if (nrow(mk)) pmn <- pmn[.norm_key(pmn$name) %in% mk$key, , drop=FALSE]
    }
    datatable(pmn, options=list(pageLength=10, scrollX=TRUE))
  })

  output$tbl_pm_edges <- renderDT({
    pme <- pm_edges_current()
    if (is.null(pme) || !is.data.frame(pme)) return(datatable(data.frame()))
    datatable(pme, options=list(pageLength=10, scrollX=TRUE))
  })

  output$tbl_match_keys <- renderDT({
    datatable(pm_match_tbl(), options=list(pageLength=10, scrollX=TRUE))
  })


  # ---- PubMed view selector (combo) ----
  output$pm_view_ui <- renderUI({
    v <- input$pm_view
    if (is.null(v) || !nzchar(v)) v <- "Network"
    if (v == "Sankey") return(uiOutput("pm_sankey_ui"))
    if (v == "SSplot") return(plotOutput("pm_ssplot_panel", height = 820))
    if (v == "SS Kano") return(plotOutput("pm_sskano_plot", height = 700))
    if (v == "Kano (value vs value2)") return(plotOutput("pm_kano_vv2_plot", height = 650))
    if (v == "Scatter (value vs value2)") return(plotOutput("pm_scatter_plot", height = 650))
    plotOutput("pm_network_plot", height = 650)
  })

  output$pm_network_plot <- renderPlot({
    nodes <- pm_nodes_current()
    edges <- pm_edges_current()
    if (is.null(nodes) || is.null(edges) || !is.data.frame(nodes) || !is.data.frame(edges)) {
      plot.new(); text(0.5,0.5,"No PubMed data. Click Reload PubMed data."); return()
    }
    e <- edges
    if ("Follower" %in% names(e) && !"follower" %in% names(e)) e$follower <- e$Follower
    if (!all(c("Leader","follower") %in% names(e))) {
      plot.new(); text(0.5,0.5,"Edges missing Leader/Follower."); return()
    }
    e$WCD <- if ("WCD" %in% names(e)) suppressWarnings(as.numeric(e$WCD)) else 1
    e$WCD[!is.finite(e$WCD)] <- 1

    nd <- nodes
    if (isTRUE(input$pm_show_matched_only)) {
      mk <- pm_match_tbl()
      if (nrow(mk)) nd <- nd[.norm_key(nd$name) %in% mk$key, , drop=FALSE]
      keep <- .norm_key(nd$name)
      e <- e[.norm_key(e$Leader) %in% keep & .norm_key(e$follower) %in% keep, , drop=FALSE]
    }

    nd$name <- as.character(nd$name)
    g <- igraph::graph_from_data_frame(e[,c("Leader","follower")], directed = TRUE, vertices = nd)
    igraph::E(g)$weight <- e$WCD
    vsz <- if ("value" %in% names(nd)) suppressWarnings(as.numeric(nd$value)) else rep(1, nrow(nd))
    vsz[!is.finite(vsz)] <- 1
    igraph::plot.igraph(g,
                        vertex.size = pmax(6, sqrt(pmax(1, vsz))) * 2,
                        vertex.label.cex = 0.8,
                        main = "PubMed Network (Top20)")
  })

  output$pm_sankey_ui <- renderUI({
    edges <- pm_edges_current()
    if (is.null(edges) || !is.data.frame(edges)) return(tags$div("No edges."))
    e <- edges
    if ("Follower" %in% names(e) && !"follower" %in% names(e)) e$follower <- e$Follower
    if (exists("render_author_sankey", mode="function")) {
      render_author_sankey(e[, intersect(c("Leader","follower","WCD"), names(e)), drop=FALSE])
    } else {
      tags$div("render_author_sankey() not found. Ensure pubmed/sankey.R is present.")
    }
  })

  output$pm_ssplot_panel <- renderPlot({
    d <- pm_nodes_current()
    if (is.null(d) || !is.data.frame(d) || !all(c("name","carac") %in% names(d))) {
      plot.new(); text(0.5,0.5,"Run PubMed pipeline first (Top20)."); return()
    }
    if (!any(c("ssi","SSi") %in% names(d))) {
      plot.new(); text(0.5,0.5,"No Silhouette (SS) yet. Click Run PubMed pipeline."); return()
    }
    if (exists("render_panel", mode="function") && exists("make_nodes0", mode="function") && exists("ensure_device", mode="function")) {
      sil_df <- d
      nodes0 <- make_nodes0(sil_df, nodes = d)
      created <- ensure_device(width=12, height=8, res=144)
      on.exit({ if (created) grDevices::dev.off() }, add = TRUE)
      render_panel(sil_df = sil_df, nodes0 = nodes0, nodes = d, font_scale = 1.3)
    } else {
      plot.new(); text(0.5,0.5,"SSplot renderer not found (pubmed/renderSSplot.R).")
    }
  })

  output$pm_sskano_plot <- renderPlot({
    d <- pm_nodes_current()
    if (is.null(d) || !is.data.frame(d) || !all(c("ssi","a_star1") %in% names(d))) {
      plot.new(); text(0.5,0.5,"No SS/a* yet. Click Run PubMed pipeline."); return()
    }
    ed <- pm_edges_current()
    if (exists("kano_plot_ss_astar", mode="function")) {
      p <- kano_plot_ss_astar(d, ed, label_size = input$pm_label_size, title_txt = "PubMed SS Kano (SS vs a*)")
      print(p)
    } else {
      plot(d$a_star1, d$ssi, xlab="a*", ylab="SS", main="PubMed SS Kano")
    }
  })

  output$pm_kano_vv2_plot <- renderPlot({
    d <- pm_nodes_current()
    if (is.null(d) || !is.data.frame(d) || !all(c("value","value2") %in% names(d))) {
      plot.new(); text(0.5,0.5,"No nodes."); return()
    }
    d <- .pubmed_add_carac(d)
    ed <- pm_edges_current()
    if (exists("kano_plot", mode="function")) {
      p <- kano_plot(d, ed, label_size = input$pm_label_size, title_txt = "PubMed Kano (value vs value2)")
      print(p)
    } else {
      plot(d$value2, d$value, xlab="value2", ylab="value", main="PubMed Kano")
    }
  })

  output$pm_scatter_plot <- renderPlot({
    d <- pm_nodes_current()
    if (is.null(d) || !is.data.frame(d) || !all(c("value","value2") %in% names(d))) {
      plot.new(); text(0.5,0.5,"No nodes."); return()
    }
    plot(d$value2, d$value, xlab="value2", ylab="value", main="PubMed Scatter (value2 vs value)")
    if ("name" %in% names(d)) text(d$value2, d$value, labels=d$name, pos=3, cex=0.7)
  })
  output$pm_plot_scatter <- renderPlot({
    pmn <- pm_nodes_current()
    if (is.null(pmn) || !is.data.frame(pmn) || !all(c("value","value2") %in% names(pmn))) {
      plot.new(); text(0.5,0.5,"No PubMed nodes loaded."); return()
    }
    d <- pmn
    if (isTRUE(input$pm_show_matched_only)) {
      mk <- pm_match_tbl()
      if (nrow(mk)) d <- d[.norm_key(d$name) %in% mk$key, , drop=FALSE]
    }
    d$value <- suppressWarnings(as.numeric(d$value))
    d$value2 <- suppressWarnings(as.numeric(d$value2))
    plot(d$value2, d$value, xlab="value2", ylab="value", main="PubMed: value vs value2")
    if ("name" %in% names(d)) text(d$value2, d$value, labels=d$name, pos=3, cex=0.7)
  })

  output$pm_plot_kano <- renderPlot({
    pmn <- pm_nodes_current()
    if (is.null(pmn) || !is.data.frame(pmn) || !all(c("value","value2") %in% names(pmn))) {
      plot.new(); text(0.5,0.5,"No PubMed nodes loaded."); return()
    }
    d <- pmn
    if (isTRUE(input$pm_show_matched_only)) {
      mk <- pm_match_tbl()
      if (nrow(mk)) d <- d[.norm_key(d$name) %in% mk$key, , drop=FALSE]
    }
    d <- .pubmed_add_carac(d)
    ed <- pm_edges_current()
    if (exists("kano_plot", mode="function")) {
      kano_plot(d, ed, title_txt = ifelse(isTRUE(input$pm_use_top20),
                                         "PubMed Kano (Top20, computed)",
                                         "PubMed Kano (raw nodes)"))
    } else {
      plot(d$value2, d$value, xlab="value2", ylab="value", main="PubMed Kano (fallback)")
      abline(v=mean(d$value2, na.rm=TRUE), h=mean(d$value, na.rm=TRUE), lty=2)
      text(d$value2, d$value, labels=d$name, pos=3, cex=0.7)
    }
  })

  output$pm_plot_ss_astar <- renderPlot({
    d <- pm_nodes_current()
    if (is.null(d) || !is.data.frame(d) || !all(c("ssi","a_star1") %in% names(d))) {
      plot.new(); text(0.5,0.5,"No SS/a* yet. Click 'Run PubMed pipeline'."); return()
    }
    ed <- pm_edges_current()
    if (exists("kano_plot_ss_astar", mode="function")) {
      kano_plot_ss_astar(d, ed, title_txt = "PubMed Kano (SS vs a*)")
    } else {
      plot(d$a_star1, d$ssi, xlab="a*", ylab="SS", main="PubMed: SS vs a* (fallback)")
      text(d$a_star1, d$ssi, labels=d$name, pos=3, cex=0.7)
      abline(v=median(d$a_star1, na.rm=TRUE), h=median(d$ssi, na.rm=TRUE), lty=2)
    }
  })

  output$pm_plot_ssplot <- renderPlot({
    d <- pm_nodes_current()
    if (is.null(d) || !is.data.frame(d) || !all(c("value","value2") %in% names(d))) {
      plot.new(); text(0.5,0.5,"No PubMed nodes loaded."); return()
    }
    plot(d$value2, d$value, xlab="value2", ylab="value",
         main = ifelse(isTRUE(input$pm_use_top20), "SSplot (Top20)", "SSplot (raw)"))
    text(d$value2, d$value, labels=d$name, pos=3, cex=0.7)
  })

  output$tbl_wos_author_top5 <- renderDT({
    s <- rv$summary
    if (is.null(s) || !is.data.frame(s)) return(datatable(data.frame()))
    a <- s[s$Domain == "Author", , drop=FALSE]
    datatable(a, options=list(dom='t', pageLength=5, scrollX=TRUE))
  })

  output$tbl_pm_matched <- renderDT({
    pmn <- pm_nodes_current()
    if (is.null(pmn) || !is.data.frame(pmn)) return(datatable(data.frame()))
    mk <- pm_match_tbl()
    if (!nrow(mk)) return(datatable(data.frame()))
    d <- pmn[.norm_key(pmn$name) %in% mk$key, , drop=FALSE]
    datatable(d, options=list(pageLength=10, scrollX=TRUE))
  })



  # ---- Domain tables (WoS summary) ----
  .dom_tbl <- function(dom, include_domain = FALSE){
    req(rv$summary)
    df <- rv$summary
    df <- df[df$Domain %in% dom, , drop=FALSE]
    if (!nrow(df)) return(data.frame(Note = "(No Data)"))
    if (!include_domain) {
      df <- df[, intersect(c("Element","n","h","AAC"), names(df)), drop=FALSE]
    } else {
      df <- df[, intersect(c("Domain","Element","n","h","AAC"), names(df)), drop=FALSE]
    }
    df
  }

  output$tbl_dom_author <- renderDT({
    datatable(.dom_tbl("Author"), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })

  output$tbl_dom_journal_year <- renderDT({
    datatable(.dom_tbl(c("Journal","Year"), include_domain = TRUE), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })

  output$tbl_dom_year_articles <- renderDT({
    datatable(.dom_tbl("Year"), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })

  output$tbl_dom_year_count <- renderDT({
    # show yearly counts from woswide32
    req(rv$woswide32)
    w <- rv$woswide32
    pick_col <- function(cands){
      nms <- names(w); nms_l <- tolower(nms); c_l <- tolower(cands)
      hit <- which(nms_l %in% c_l); if (length(hit)) return(nms[hit[1]])
      for (c in c_l){ h <- which(grepl(c, nms_l, fixed = TRUE)); if (length(h)) return(nms[h[1]]) }
      NA_character_
    }
    col_year <- pick_col(c("publication year","publication.year","year"))
    if (is.na(col_year)) return(datatable(data.frame(Note="No Year column."), options=list(dom='t')))
    yy <- suppressWarnings(as.integer(as.character(w[[col_year]])))
    tab <- as.data.frame(table(yy), stringsAsFactors = FALSE)
    names(tab) <- c("Year","Count")
    tab <- tab[order(suppressWarnings(as.integer(tab$Year))), , drop=FALSE]
    datatable(tab, options=list(pageLength=10, scrollX=TRUE))
  })

  output$tbl_dom_term_year <- renderDT({
    datatable(.dom_tbl("Term/Year"), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })

  output$tbl_dom_country <- renderDT({
    datatable(.dom_tbl("Country"), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })
  output$tbl_dom_institute <- renderDT({
    datatable(.dom_tbl("Institute"), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })
  output$tbl_dom_department <- renderDT({
    datatable(.dom_tbl("Department"), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })
  
  output$tbl_dom_mesh <- renderDT({
    datatable(.dom_tbl("MeSH"), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })

output$tbl_dom_state <- renderDT({
    datatable(.dom_tbl("State/Province"), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })
  output$tbl_dom_keyword <- renderDT({
    datatable(.dom_tbl("Keyword"), options = list(dom='t', pageLength=10, scrollX=TRUE))
  })

  # --------------------------
  # WoS visuals (match PubMed-style outputs)
  # --------------------------
  .pick_wos_col <- function(df, cands){
    nms <- names(df); nms_l <- tolower(nms); c_l <- tolower(cands)
    hit <- which(nms_l %in% c_l)
    if (length(hit)) return(nms[hit[1]])
    for (c in c_l){
      h <- which(grepl(c, nms_l, fixed = TRUE))
      if (length(h)) return(nms[h[1]])
    }
    NA_character_
  }

  .split_semis <- function(x){
    if (is.null(x) || !length(x)) return(character(0))
    x <- as.character(x)
    x <- trimws(x)
    if (!nzchar(x)) return(character(0))
    parts <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    parts[nzchar(parts)]
  }

  .build_author_network <- reactive({
    # Prefer woslong32 (has FA/CA parsed) when available; otherwise fall back to woswide32
    w <- NULL
    if (!is.null(rv$woslong32) && is.data.frame(rv$woslong32) && nrow(rv$woslong32) > 0) {
      w <- rv$woslong32
    } else {
      req(rv$woswide32)
      w <- rv$woswide32
    }

    # Prefer FAauthor / Caauthor; fallback to Authors first/last
    c_fa   <- .pick_wos_col(w, c("FAauthor","FA author","first author"))
    c_ca   <- .pick_wos_col(w, c("CAauthor","Caauthor","CA author","corresponding author"))
    c_auth <- .pick_wos_col(w, c("Authors","Author","AU"))

    fa <- if (!is.na(c_fa)) as.character(w[[c_fa]]) else rep(NA_character_, nrow(w))
    ca <- if (!is.na(c_ca)) as.character(w[[c_ca]]) else rep(NA_character_, nrow(w))

    if ((is.na(c_fa) || is.na(c_ca)) && !is.na(c_auth)) {
      au <- as.character(w[[c_auth]])
      fa2 <- vapply(au, function(z){ p <- .split_semis(z); if (length(p)) p[1] else NA_character_ }, FUN.VALUE = character(1))
      ca2 <- vapply(au, function(z){ p <- .split_semis(z); if (length(p)) p[length(p)] else NA_character_ }, FUN.VALUE = character(1))
      if (is.na(c_fa)) fa <- fa2
      if (is.na(c_ca)) ca <- ca2
    }

    df <- data.frame(from = fa, to = ca, stringsAsFactors = FALSE)
    df <- df[!is.na(df$from) & !is.na(df$to) & nzchar(df$from) & nzchar(df$to), , drop = FALSE]
    if (!nrow(df)) return(list(nodes = data.frame(), edges = data.frame(), modes = NULL, edges_full = data.frame()))

    # Keep full FA/CA pairs for document-count (self-loops included there)
    df_all <- df

    # Network edges exclude self-loops
    df <- df[df$from != df$to, , drop = FALSE]
    df$w <- 1L
    edges <- df %>%
      dplyr::group_by(from, to) %>%
      dplyr::summarise(value = sum(w), .groups = "drop")

    # --- FLCA + Major Sampling + Silhouette (Top20) ---
    # Build nodes input for FLCA runner
    all_names <- sort(unique(c(df_all$from, df_all$to)))

    # document count: prefer WoS summary table (Author domain) so values match the Top20 table exactly.
    # Fallback to FA/CA paper counts if summary is not ready.
    doc_map <- setNames(numeric(length(all_names)), all_names)
    sum_a <- tryCatch({
      if (!is.null(rv$summary) && is.data.frame(rv$summary)) rv$summary[rv$summary$Domain == "Author", , drop = FALSE] else NULL
    }, error = function(e) NULL)
    if (is.data.frame(sum_a) && nrow(sum_a) > 0 && all(c("Element", "n") %in% names(sum_a))) {
      nm <- trimws(as.character(sum_a$Element))
      nv <- suppressWarnings(as.numeric(sum_a$n))
      ok <- nzchar(nm) & is.finite(nv)
      mp <- setNames(nv[ok], nm[ok])
      hit <- intersect(all_names, names(mp))
      doc_map[hit] <- mp[hit]
    } else {
      for (ii in seq_len(nrow(df_all))) {
        ppl <- unique(c(as.character(df_all$from[ii]), as.character(df_all$to[ii])))
        ppl <- ppl[!is.na(ppl) & nzchar(ppl)]
        if (length(ppl)) doc_map[ppl] <- doc_map[ppl] + 1
      }
    }

    # value2 = sum(edge) excluding self-loops, not igraph strength label
    edge_sum <- setNames(numeric(length(all_names)), all_names)
    if (nrow(edges)) {
      for (ii in seq_len(nrow(edges))) {
        f <- as.character(edges$from[ii]); t <- as.character(edges$to[ii]); wv <- suppressWarnings(as.numeric(edges$value[ii]))
        if (!is.finite(wv)) wv <- 0
        if (nzchar(f)) edge_sum[f] <- edge_sum[f] + wv
        if (nzchar(t)) edge_sum[t] <- edge_sum[t] + wv
      }
    }

    nodes_in <- data.frame(
      name = all_names,
      value = as.numeric(doc_map[all_names]),
      value2 = as.numeric(edge_sum[all_names]),
      stringsAsFactors = FALSE
    )
    nodes_in$value[!is.finite(nodes_in$value)] <- 0
    nodes_in$value2[!is.finite(nodes_in$value2)] <- 0

    # edges input: Leader/Follower/WCD
    edges_in <- data.frame(
      Leader = as.character(edges$from),
      Follower = as.character(edges$to),
      WCD = as.numeric(edges$value),
      stringsAsFactors = FALSE
    )

    fl <- NULL
    # Prefer official runner if present; otherwise use fallback run_flca_ms_sil() defined in this app
    if (exists("run_flca_ms_sil_runner", mode = "function")) {
      fl <- tryCatch(
        run_flca_ms_sil_runner(nodes_in, edges_in, cfg = list(target_n = 20), verbose = FALSE),
        error = function(e) e
      )
    } else if (exists("run_flca_ms_sil", mode = "function")) {
      fl <- tryCatch(
        run_flca_ms_sil(nodes_in, edges_in, cfg = list(target_n = 20), verbose = FALSE),
        error = function(e) e
      )
    } else {
      fl <- simpleError("FLCA runner not available (run_flca_ms_sil_runner / run_flca_ms_sil missing).")
    }

    if (inherits(fl, "error")) {
      # If FLCA fails, still show a usable (unclustered) network rather than blank
      nodes <- data.frame(
        id = all_names,
        label = all_names,
        value = as.numeric(doc_map[all_names]),
        value2 = as.numeric(edge_sum[all_names]),
        group = "NA",
        title = paste0(
          all_names,
          "<br>n=", as.integer(round(as.numeric(doc_map[all_names]))),
          "<br>sum(edge)=", round(as.numeric(edge_sum[all_names]), 3)
        ),
        stringsAsFactors = FALSE
      )
      # keep top 20 by document count, then sum(edge)
      ord <- order(nodes$value, nodes$value2, decreasing = TRUE)
      nodes <- nodes[ord[seq_len(min(20, nrow(nodes)))], , drop = FALSE]
      edges2 <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, , drop = FALSE]
      ed <- data.frame(from = edges2$from, to = edges2$to, value = edges2$value, stringsAsFactors = FALSE)
      return(list(nodes = nodes, edges = ed, modes = NULL, edges_full = edges_in))
    }

    modes <- fl$modes
    edges20 <- fl$data  # Leader/Follower/WCD on Top20

    # Sync FLCA modes back to app-defined values: value=document count, value2=sum(edge)
    if (is.data.frame(modes) && nrow(modes)) {
      ixm <- match(trimws(as.character(modes$name)), trimws(as.character(nodes_in$name)))
      modes$value <- as.numeric(nodes_in$value[ixm])
      modes$value2 <- as.numeric(nodes_in$value2[ixm])
    }

    # nodes for visNetwork (Top20)
    modes <- modes[modes$name %in% unique(c(edges20$Leader, edges20$Follower)), , drop = FALSE]
    nodes <- data.frame(
      id = modes$name,
      label = modes$name,
      value = as.numeric(modes$value),
      group = as.character(modes$carac),
      title = paste0(
        modes$name,
        "<br>cluster=C", modes$carac,
        "<br>n=", round(as.numeric(modes$value), 3),
        if ("value2" %in% names(modes)) paste0("<br>sum(edge)=", round(as.numeric(modes$value2), 3)) else "",
        if ("SSi" %in% names(modes)) paste0("<br>SS=", round(as.numeric(modes$SSi), 3)) else "",
        if ("a_star" %in% names(modes)) paste0("<br>a*=", round(as.numeric(modes$a_star), 3)) else ""
      ),
      stringsAsFactors = FALSE
    )

    ed <- data.frame(
      from = as.character(edges20$Leader),
      to   = as.character(edges20$Follower),
      value = as.numeric(edges20$WCD),
      stringsAsFactors = FALSE
    )

    list(nodes = nodes, edges = ed, modes = modes, edges_full = edges_in)
  })


.build_journal_year_network <- reactive({
  # Prefer FLCA/MS/SS result from metatable (matches PubMed visuals)
  if (!is.null(rv$dom_journal_year) && isTRUE(rv$dom_journal_year$ok)) {
    return(.as_vis_nodes_edges(rv$dom_journal_year))
  }

  # Fallback: simple bipartite from woswide32
  req(rv$woswide32)
  w <- rv$woswide32

  c_j <- .pick_wos_col(w, c("Journal", "journal", "Source Title", "SO"))
  c_y <- .pick_wos_col(w, c("Publication Year", "Publication.Year", "Year", "year", "PY"))
  if (is.na(c_j) || is.na(c_y)) return(list(nodes = data.frame(), edges = data.frame()))

  j <- as.character(w[[c_j]])
  y <- as.character(w[[c_y]])
  df <- data.frame(journal = j, year = y, stringsAsFactors = FALSE)
  df <- df[nzchar(df$journal) & nzchar(df$year), , drop = FALSE]
  if (!nrow(df)) return(list(nodes = data.frame(), edges = data.frame()))

  df$w <- 1L
  edges <- df %>%
    dplyr::group_by(journal, year) %>%
    dplyr::summarise(value = sum(w), .groups = "drop")

  topj <- edges %>%
    dplyr::group_by(journal) %>%
    dplyr::summarise(n = sum(value), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::slice_head(n = 20)

  edges2 <- edges %>% dplyr::semi_join(topj, by = "journal")

  nodes <- dplyr::bind_rows(
    data.frame(id = edges2$journal, label = edges2$journal, group = "Journal", stringsAsFactors = FALSE),
    data.frame(id = edges2$year,    label = edges2$year,    group = "Year",    stringsAsFactors = FALSE)
  ) %>% dplyr::distinct(id, .keep_all = TRUE)

  ed <- data.frame(from = edges2$journal, to = edges2$year, value = edges2$value, stringsAsFactors = FALSE)
  list(nodes = nodes, edges = ed)
})


.build_dom_network <- function(dom_key){
  req(rv$metatable)
  res <- NULL
  if (dom_key == "journal_year") res <- rv$dom_journal_year
  if (dom_key == "year_articles") res <- rv$dom_year_articles
  if (dom_key == "term_year") res <- rv$dom_term_year
  if (dom_key == "country") res <- rv$dom_country
  if (dom_key == "institute") res <- rv$dom_institute
  if (dom_key == "department") res <- rv$dom_department
  if (dom_key == "keyword") res <- rv$dom_keyword
  if (dom_key == "mesh") res <- rv$dom_mesh
  if (dom_key == "region") res <- rv$dom_region
  .as_vis_nodes_edges(res)
}

  
  # ---- Combo network (any domain) ----
  observeEvent(list(rv$metatable, input$combo_domain), {
    req(rv$metatable, input$combo_domain)
    rv$combo_res <- .run_domain(rv$metatable, input$combo_domain)
  }, ignoreInit = FALSE)

  # Sync chord domain selector with combo network selector
  observeEvent(input$chord_domain, {
    req(input$chord_domain)
    if (isTRUE(input$chord_follow_combo)) {
      updateSelectInput(session, "combo_domain", selected = input$chord_domain)
    } else {
      # run chord domain without touching combo
      req(rv$metatable)
      rv$chord_res <- .run_domain(rv$metatable, input$chord_domain)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$combo_domain, {
    # keep chord selector in sync (when user changes combo)
    updateSelectInput(session, "chord_domain", selected = input$combo_domain)
  }, ignoreInit = FALSE)

  output$vn_combo <- renderVisNetwork({
    req(rv$combo_res)
    res <- rv$combo_res
    if (is.null(res) || isFALSE(res$ok)) {
      return(NULL)
    }
    g <- .as_vis_nodes_edges(res)
    if (!nrow(g$nodes)) return(NULL)

    visNetwork(g$nodes, g$edges, height = "650px") %>%
      visEdges(smooth = TRUE) %>%
      visNodes(font = list(size = 24, face = "arial", strokeWidth = 2, strokeColor = "white"), scaling = list(min = 10, max = 60)) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = FALSE) %>%
      visInteraction(hover = TRUE, tooltipDelay = 80)
  })

  output$combo_nodes <- renderDT({
    req(rv$combo_res)
    res <- rv$combo_res
    if (is.null(res) || isFALSE(res$ok)) return(datatable(data.frame()))
    nodes <- .normalize_modes_for_table(res$modes)
    datatable(nodes, options = list(pageLength = 20, scrollX = TRUE))
  })

  output$combo_edges <- renderDT({
    req(rv$combo_res)
    res <- rv$combo_res
    if (is.null(res) || isFALSE(res$ok) || !nrow(res$edges20)) return(datatable(data.frame()))
    datatable(res$edges20, options = list(pageLength = 20, scrollX = TRUE))
  })



# ---- Chord diagram (uses SAME Top-20 nodes/edges as Combo Network) ----
.build_chord_matrix <- function(nodes, edges20) {
  if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) < 2) return(NULL)
  if (is.null(edges20) || !is.data.frame(edges20) || nrow(edges20) < 1) return(NULL)

  # clean labels to avoid mismatch (e.g., "name<br>cluster=...")
  nodes$name <- gsub("<br>.*$", "", as.character(nodes$name))
  nodes <- nodes[!is.na(nodes$name) & nzchar(nodes$name), , drop = FALSE]
  if (nrow(nodes) < 2) return(NULL)

  # ensure unique node names for circlize sectors
  nodes <- nodes[!duplicated(nodes$name), , drop = FALSE]
  if (nrow(nodes) < 2) return(NULL)

  # normalize edges columns
  nms <- names(edges20)
  from_col <- intersect(nms, c("Leader","from","From","Source","source","node1","Node1"))[1]
  to_col   <- intersect(nms, c("Follower","to","To","Target","target","node2","Node2"))[1]
  w_col    <- intersect(nms, c("WCD","weight","Weight","value","Value"))[1]
  if (is.na(from_col) || is.na(to_col)) return(NULL)
  if (is.na(w_col)) { edges20$..w <- 1; w_col <- "..w" }

  from <- as.character(edges20[[from_col]])
  to   <- as.character(edges20[[to_col]])
  w    <- suppressWarnings(as.numeric(edges20[[w_col]]))
  w[!is.finite(w)] <- 0

  # filter to nodes set + positive weights
  keep <- from %in% nodes$name & to %in% nodes$name & w > 0
  from <- from[keep]; to <- to[keep]; w <- w[keep]
  if (!length(w)) return(NULL)

  # build symmetric matrix (undirected)
  mat <- matrix(0, nrow = nrow(nodes), ncol = nrow(nodes),
                dimnames = list(nodes$name, nodes$name))
  for (i in seq_along(w)) {
    a <- match(from[i], nodes$name)
    b <- match(to[i], nodes$name)
    if (!is.na(a) && !is.na(b) && a != b) {
      mat[a, b] <- mat[a, b] + w[i]
      mat[b, a] <- mat[b, a] + w[i]
    }
  }
  mat
}

# Safe wrapper for chorddiag across versions (prevents "argument 2 matches multiple formal arguments")
.safe_chorddiag_widget <- function(mat, group, groupColors, groupnamePadding = 20) {
  if (!requireNamespace("chorddiag", quietly = TRUE)) return(NULL)
  fn <- chorddiag::chorddiag
  fml <- names(formals(fn))

  args <- list()

  # first argument name varies across versions; prefer 'x' if present
  if ("x" %in% fml) {
    args$x <- mat
  } else if ("mat" %in% fml) {
    args$mat <- mat
  } else {
    args[[1]] <- mat
  }

  # group / colors (exact names only; try a few variants)
  if ("group" %in% fml) args$group <- group

  if ("groupColors" %in% fml) {
    args$groupColors <- groupColors
  } else if ("groupColours" %in% fml) {
    args$groupColours <- groupColors
  } else if ("col" %in% fml) {
    args$col <- groupColors
  }

  if ("groupnamePadding" %in% fml) {
    args$groupnamePadding <- groupnamePadding
  } else if ("groupNamePadding" %in% fml) {
    args$groupNamePadding <- groupnamePadding
  }

  do.call(fn, args)
}


.cluster_colors_from_carac <- function(nodes) {
  if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) < 1) return(NULL)

  # accept either carac or cluster column name
  if (!("carac" %in% names(nodes)) && ("cluster" %in% names(nodes))) nodes$carac <- nodes$cluster
  if (!("carac" %in% names(nodes))) nodes$carac <- "C1"

  # node names
  if (!("name" %in% names(nodes))) {
    if ("label" %in% names(nodes)) nodes$name <- as.character(nodes$label)
    else if ("id" %in% names(nodes)) nodes$name <- as.character(nodes$id)
    else nodes$name <- as.character(nodes[[1]])
  }
  nodes$name <- gsub("<br>.*$", "", as.character(nodes$name))
  nodes$carac <- as.character(nodes$carac)
  nodes$carac[is.na(nodes$carac) | !nzchar(nodes$carac)] <- "C1"

  # ensure unique node names; keep first occurrence
  keep <- !is.na(nodes$name) & nzchar(nodes$name) & !duplicated(nodes$name)
  nodes <- nodes[keep, , drop = FALSE]

  pal <- .cluster_palette_map(nodes$carac)
  cols <- unname(pal[nodes$carac])
  names(cols) <- nodes$name
  cols
}



# ---- Chord payload (matrix + cluster colors aligned to sectors) ----
chord_payload <- reactive({
  req(if (isTRUE(input$chord_follow_combo)) rv$combo_res else rv$chord_res)
  res <- if (isTRUE(input$chord_follow_combo)) rv$combo_res else rv$chord_res
  if (is.null(res) || isFALSE(res$ok)) return(list(ok=FALSE, reason="no data"))
  nodes <- res$modes
  edges20 <- res$edges20
  mat <- .build_chord_matrix(nodes, edges20)
  if (is.null(mat) || nrow(mat) < 2 || sum(mat) <= 0) return(list(ok=FALSE, reason="not enough"))

  # colors from nodes$carac aligned to rownames(mat)
  cols0 <- .cluster_colors_from_carac(nodes)
  rn <- rownames(mat)
  cols <- NULL
  if (!is.null(cols0) && length(cols0)) {
    cols <- cols0[match(rn, names(cols0))]
    names(cols) <- rn
    cols[is.na(cols)] <- "#B0B0B0"
  }

  # map node -> cluster -> color (for debug/log)
  # derive cluster per node name
  nodes2 <- nodes
  if (!is.null(nodes2) && is.data.frame(nodes2) && "name" %in% names(nodes2)) {
    nodes2$name <- gsub("<br>.*$", "", as.character(nodes2$name))
  }
  if (!("carac" %in% names(nodes2)) && ("cluster" %in% names(nodes2))) nodes2$carac <- nodes2$cluster
  if (!("carac" %in% names(nodes2))) nodes2$carac <- "C1"
  carac_by_name <- setNames(as.character(nodes2$carac), as.character(nodes2$name))
  df_map <- data.frame(
    node = rn,
    carac = unname(carac_by_name[rn]),
    color = if (is.null(cols)) NA_character_ else unname(cols),
    stringsAsFactors = FALSE
  )
  df_map$carac[is.na(df_map$carac) | !nzchar(df_map$carac)] <- "C1"

  list(ok=TRUE, mat=mat, cols=cols, map=df_map, edges_cols=names(edges20), n_nodes=nrow(nodes), n_edges=nrow(edges20))
})

# log chord cluster mapping once per domain+matrix signature
observeEvent(chord_payload(), {
  payload <- chord_payload()
  if (is.null(payload) || isFALSE(payload$ok)) return()
  dom <- if (isTRUE(input$chord_follow_combo)) input$combo_domain else input$chord_domain
  sig <- paste0(dom, "|", nrow(payload$mat), "|", sum(payload$mat))
  if (!identical(rv$chord_sig, sig)) {
    rv$chord_sig <- sig
    .safe_addlog(paste0("[CHORD] domain=", dom, " sectors=", nrow(payload$mat), " sumW=", sum(payload$mat)))
    # log top sectors with cluster/color
    mm <- head(payload$map, 20)
    for (i in seq_len(nrow(mm))) {
      .safe_addlog(paste0("[CHORD] ", mm$node[i], " | carac=", mm$carac[i], " | color=", mm$color[i]))
    }
  }
}, ignoreInit = TRUE)


output$chord_debug <- renderText({
  payload <- chord_payload()
  if (is.null(payload) || isFALSE(payload$ok)) {
    res <- if (isTRUE(input$chord_follow_combo)) rv$combo_res else rv$chord_res
    nodes <- if (is.null(res)) NULL else res$modes
    edges20 <- if (is.null(res)) NULL else res$edges20
    return(paste0(
      "Chord: not enough data. ",
      "nodes=", if (is.null(nodes)) "NULL" else nrow(nodes),
      ", edges20=", if (is.null(edges20)) "NULL" else nrow(edges20),
      ", edges cols=", if (is.null(edges20)) "NULL" else paste(names(edges20), collapse=",")
    ))
  }
  dom <- if (isTRUE(input$chord_follow_combo)) input$combo_domain else input$chord_domain
  paste0(
    "Chord OK | domain=", dom,
    " | sectors=", nrow(payload$mat),
    " | sumW=", sum(payload$mat),
    " | colored=", if (is.null(payload$cols)) 0 else sum(!is.na(payload$cols)),
    " | edges cols=", paste(payload$edges_cols, collapse=",")
  )
})

output$chord_ui <- renderUI({
  req(if (isTRUE(input$chord_follow_combo)) rv$combo_res else rv$chord_res)
  res <- if (isTRUE(input$chord_follow_combo)) rv$combo_res else rv$chord_res
  if (is.null(res) || isFALSE(res$ok)) {
    return(tags$div("Chord: no data yet (choose a domain)."))
  }

  nodes <- res$modes
  edges20 <- res$edges20

  mat <- .build_chord_matrix(nodes, edges20)
  if (is.null(mat) || nrow(mat) < 2 || sum(mat) <= 0) {
    return(tags$div("Chord: not enough data (need Top-20 nodes + edges). If Combo Network is visible, this usually means name mismatch or all edges were filtered out."))
  }

  # choose rendering path
  if (requireNamespace("chorddiag", quietly = TRUE)) {
    # interactive chorddiag (cluster-aware via group + groupColors)
    payload <- chord_payload()
    validate(need(!is.null(payload) && isTRUE(payload$ok), "Chord: not enough data."))

    rn <- rownames(payload$mat)

    # group per sector (cluster label) aligned to rn
    grp <- payload$map$carac
    names(grp) <- payload$map$node
    grp <- as.character(grp[rn])
    grp[is.na(grp) | !nzchar(grp)] <- "C1"
    grp_levels <- unique(grp)

    # sector colors aligned to rn (fallback gray)
    sector_cols <- payload$cols
    if (is.null(sector_cols) || !length(sector_cols)) {
      sector_cols <- rep("#B0B0B0", length(rn))
      names(sector_cols) <- rn
    } else {
      sector_cols <- sector_cols[rn]
      sector_cols[is.na(sector_cols) | !nzchar(sector_cols)] <- "#B0B0B0"
      names(sector_cols) <- rn
    }

    # groupColors: one color per group (take first sector color in that group)
    group_cols <- vapply(grp_levels, function(g) {
      idx <- which(grp == g)[1]
      unname(sector_cols[idx])
    }, FUN.VALUE = character(1))

    out <- tryCatch({
      .safe_chorddiag_widget(
        mat = payload$mat,
        group = factor(grp, levels = grp_levels),
        groupColors = unname(group_cols),
        groupnamePadding = 20
      )
    }, error = function(e) {
      .safe_addlog(paste0("[CHORD][ERROR] chorddiag failed: ", conditionMessage(e)))
      NULL
    })
    validate(need(!is.null(out), "Chord: chorddiag failed (see Log tab)."))

    htmlwidgets::prependContent(out, htmltools::tags$style("svg text { font-family: Arial; }"))
  } else if (requireNamespace("circlize", quietly = TRUE)) {
    # static fallback
    plotOutput("chord_plot", height = "700px")
  } else {
    tags$div(
      "Chord: need package 'chorddiag' (interactive) or 'circlize' (static).",
      tags$br(),
      tags$code("install.packages('chorddiag')"),
      "  or  ",
      tags$code("install.packages('circlize')")
    )
  }
})


output$chord_plot <- renderPlot({
  payload <- chord_payload()
  validate(need(!is.null(payload) && isTRUE(payload$ok), "Chord: not enough data."))
  validate(need(requireNamespace("circlize", quietly = TRUE), "Install 'circlize' for static chord."))

  circlize::circos.clear()
  cols <- payload$cols
  rn <- rownames(payload$mat)
  if (!is.null(cols) && length(cols)) {
    cols <- cols[rn]
    cols[is.na(cols) | !nzchar(cols)] <- "#B0B0B0"
    names(cols) <- rn
  }
  circlize::chordDiagram(payload$mat, grid.col = cols, annotationTrack = "grid")
})

output$vn_author <- renderVisNetwork({
    net <- .build_author_network()
    validate(need(nrow(net$nodes) > 0, "No Author network yet. Upload WoS and run routines first."))
    visNetwork(net$nodes, net$edges, height = "650px") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visInteraction(navigationButtons = TRUE)
  })

  output$tbl_author_nodes <- renderDT({
    net <- .build_author_network()
    tbl <- if (!is.null(net$modes) && is.data.frame(net$modes) && nrow(net$modes)) normalize_modes_for_table(net$modes) else net$nodes
    datatable(tbl, options = list(dom='t', pageLength=10, scrollX=TRUE))
  })

  output$tbl_author_edges <- renderDT({
    net <- .build_author_network()
    datatable(net$edges, options = list(dom='t', pageLength=10, scrollX=TRUE))
  })

  output$vn_journal_year <- renderVisNetwork({
    net <- .build_journal_year_network()
    validate(need(nrow(net$nodes) > 0, "No Journal/Year network yet. Upload WoS and run routines first."))
    visNetwork(net$nodes, net$edges, height = "520px") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visInteraction(navigationButtons = TRUE)
  })

  output$tbl_journal_year <- renderDT({
    net <- .build_journal_year_network()
    datatable(net$edges, options = list(dom='t', pageLength=10, scrollX=TRUE))
  })


# ---- Domain networks (WoS metatable -> FLCA/MS/SS), mirrored from PubMed visuals style ----
output$vn_year_articles <- renderVisNetwork({
  net <- .build_dom_network("year_articles")
  validate(need(nrow(net$nodes) > 0, "No Year/Articles network yet. Run routines first."))
  visNetwork(net$nodes, net$edges, height = "520px") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visInteraction(navigationButtons = TRUE)
})

output$vn_term_year <- renderVisNetwork({
  net <- .build_dom_network("term_year")
  validate(need(nrow(net$nodes) > 0, "No Term/Year network yet. Run routines first."))
  visNetwork(net$nodes, net$edges, height = "520px") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visInteraction(navigationButtons = TRUE)
})

output$vn_country <- renderVisNetwork({
  net <- .build_dom_network("country")
  validate(need(nrow(net$nodes) > 0, "No Country network yet. Run routines first."))
  visNetwork(net$nodes, net$edges, height = "520px") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visInteraction(navigationButtons = TRUE)
})

output$vn_institute <- renderVisNetwork({
  net <- .build_dom_network("institute")
  validate(need(nrow(net$nodes) > 0, "No Institute network yet. Run routines first."))
  visNetwork(net$nodes, net$edges, height = "520px") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visInteraction(navigationButtons = TRUE)
})

output$vn_department <- renderVisNetwork({
  net <- .build_dom_network("department")
  validate(need(nrow(net$nodes) > 0, "No Department network yet. Run routines first."))
  visNetwork(net$nodes, net$edges, height = "520px") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visInteraction(navigationButtons = TRUE)
})


output$vn_mesh <- renderVisNetwork({
  net <- .build_dom_network("mesh")
  validate(need(nrow(net$nodes) > 0, "No MeSH network yet. (Need PubMed terms loaded to map WoS keywords.)"))
  visNetwork(net$nodes, net$edges, height = "520px") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visInteraction(navigationButtons = TRUE)
})

output$vn_keyword <- renderVisNetwork({
  net <- .build_dom_network("keyword")
  validate(need(nrow(net$nodes) > 0, "No Keyword network yet. Run routines first."))
  visNetwork(net$nodes, net$edges, height = "520px") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visInteraction(navigationButtons = TRUE)
})



# --------------------------
# Kano (WoS): combo selection over metadomains (mirrors apppubmed.R behavior)
# --------------------------
.wos_domain_key <- function(dom){
  dom <- as.character(dom %||% "Author")
  if (dom == "Author")        return("dom_author")
  if (dom == "Journal/Year")  return("dom_journal_year")
  if (dom == "Year/Articles") return("dom_year_articles")
  if (dom == "Term/Year")     return("dom_term_year")
  if (dom == "Country")       return("dom_country")
  if (dom == "State/Province")return("dom_region")
  if (dom == "Institute")     return("dom_institute")
  if (dom == "Department")    return("dom_department")
  if (dom == "Keyword")       return("dom_keyword")
  if (dom == "MeSH")          return("dom_mesh")
  "dom_author"
}

.get_wos_domain <- function(dom){
  key <- .wos_domain_key(dom)
  res <- rv[[key]]
  if (is.list(res) && isTRUE(res$ok) && is.data.frame(res$modes)) {
    return(list(nodes = res$modes, edges = res$edges20))
  }
  # fallback for Author only (use FA/CA graph)
  if (identical(as.character(dom), "Author")) {
    net <- tryCatch(.build_author_network(), error = function(e) NULL)
    if (is.list(net) && is.data.frame(net$modes)) {
      return(list(nodes = net$modes, edges = net$edges_full))
    }
  }
  list(nodes = data.frame(), edges = data.frame())
}

.kano_scatter <- function(nd, xcol, ycol, title_txt, label_size = 4, xlab = NULL, ylab = NULL){
  if (!is.data.frame(nd) || nrow(nd) < 2) {
    plot.new(); text(0.5,0.5,"No data (run routines first)"); return(invisible())
  }
  # normalize name column
  if (!("name" %in% names(nd)) && ("id" %in% names(nd))) nd$name <- as.character(nd$id)
  if (!("name" %in% names(nd))) nd$name <- as.character(seq_len(nrow(nd)))

  # numeric conversions
  xx <- suppressWarnings(as.numeric(nd[[xcol]]))
  yy <- suppressWarnings(as.numeric(nd[[ycol]]))
  ok <- is.finite(xx) & is.finite(yy)
  if (sum(ok) < 2) { plot.new(); text(0.5,0.5,"Kano: no finite points"); return(invisible()) }

  # bubble size
  vv <- if ("value" %in% names(nd)) suppressWarnings(as.numeric(nd$value)) else rep(1, nrow(nd))
  vv[!is.finite(vv)] <- 1
  cex <- pmax(0.7, sqrt(pmax(1, vv))/4)

  plot(xx[ok], yy[ok],
       xlab = xlab %||% xcol, ylab = ylab %||% ycol,
       main = title_txt, pch = 19, cex = cex[ok])

  # quadrant lines (Kano-style)
  abline(v = median(xx[ok], na.rm = TRUE), h = median(yy[ok], na.rm = TRUE), lty = 2)

  # labels
  if (!is.null(label_size) && is.finite(label_size) && label_size > 0) {
    text(xx[ok], yy[ok], labels = nd$name[ok], pos = 4, cex = (label_size/10) + 0.3)
  }
  invisible()
}

output$kano_aac_header <- renderText({
  d <- .get_wos_domain(input$kano_domain)
  nd <- d$nodes
  if (!is.data.frame(nd) || !nrow(nd)) return("No data yet.")
  # show a lightweight header (avoid depending on AAC functions)
  dom <- input$kano_domain %||% "Author"
  sprintf("%s: n=%d (Top20).", dom, nrow(nd))
})

output$kano_vv2_plot <- renderPlot({
  d <- .get_wos_domain(input$kano_domain)
  nd <- d$nodes
  ed <- d$edges

  validate(need(is.data.frame(nd) && nrow(nd) > 1, "Run routines first (need Top20 nodes)."))

  # ensure value2 exists
  if (!("value2" %in% names(nd))) {
    if (is.data.frame(ed) && nrow(ed) && all(c("Leader","Follower","WCD") %in% names(ed))) {
      g0 <- tryCatch(igraph::graph_from_data_frame(ed[,c("Leader","Follower","WCD")], directed = FALSE), error=function(e) NULL)
      if (!is.null(g0)) {
        st <- tryCatch(igraph::strength(g0, weights = igraph::E(g0)$WCD), error=function(e) NULL)
        if (!is.null(st)) nd$value2 <- as.numeric(st[nd$name]) else nd$value2 <- nd$value
      } else nd$value2 <- nd$value
    } else {
      nd$value2 <- nd$value
    }
  }

  # prefer pretty kano.R functions if available
  if (exists("kano_plot", mode="function")) {
    # kano_plot expects 'carac' sometimes
    if (!("carac" %in% names(nd))) nd$carac <- 1L
    p <- tryCatch(kano_plot(nd, ed, label_size = input$kano_label_size, title_txt = sprintf("Kano (%s): value vs value2", input$kano_domain)), error=function(e) NULL)
    if (!is.null(p)) { try(print(p), silent=TRUE); return(invisible()) }
  }

  .kano_scatter(nd, xcol = "value2", ycol = "value",
               title_txt = sprintf("Kano (%s): value vs value2", input$kano_domain),
               label_size = input$kano_label_size, xlab="value2", ylab="value")
}, height = 650)

output$kano_ss_astar_plot <- renderPlot({
  d <- .get_wos_domain(input$kano_domain)
  nd <- d$nodes
  validate(need(is.data.frame(nd) && nrow(nd) > 1, "Run routines first (need Top20 nodes)."))

  # normalize SS and a*
  if (!("ssi" %in% names(nd)) && ("SSi" %in% names(nd))) nd$ssi <- nd$SSi
  if (!("a_star1" %in% names(nd)) && ("a_star" %in% names(nd))) nd$a_star1 <- nd$a_star

  # fallback if a_star1 missing
  if (!("a_star1" %in% names(nd))) {
    if ("a_i" %in% names(nd)) nd$a_star1 <- 1/(1+as.numeric(nd$a_i)) else nd$a_star1 <- NA_real_
  }

  if (exists("kano_plot_ss_astar", mode="function")) {
    if (!("carac" %in% names(nd))) nd$carac <- 1L
    p <- tryCatch(kano_plot_ss_astar(nd, d$edges, label_size = input$kano_label_size,
                                    title_txt = sprintf("Kano (%s): SS vs a*", input$kano_domain)),
                  error=function(e) NULL)
    if (!is.null(p)) { try(print(p), silent=TRUE); return(invisible()) }
  }

  .kano_scatter(nd, xcol = "a_star1", ycol = "ssi",
               title_txt = sprintf("Kano (%s): SS vs a*", input$kano_domain),
               label_size = input$kano_label_size, xlab="a*", ylab="SS")
}, height = 650)



# =========================
# Sankey / SSplot / SS Kano (combo selection)
# =========================

.as_sankey_edges <- function(ed){
  if (!is.data.frame(ed) || !nrow(ed)) return(NULL)
  # accept either Leader/Follower/WCD or from/to/weight
  if (all(c("Leader","Follower","WCD") %in% names(ed))) {
    data.frame(source=as.character(ed$Leader), target=as.character(ed$Follower),
               value=suppressWarnings(as.numeric(ed$WCD)), stringsAsFactors = FALSE)
  } else if (all(c("from","to","weight") %in% names(ed))) {
    data.frame(source=as.character(ed$from), target=as.character(ed$to),
               value=suppressWarnings(as.numeric(ed$weight)), stringsAsFactors = FALSE)
  } else if (all(c("source","target","value") %in% names(ed))) {
    data.frame(source=as.character(ed$source), target=as.character(ed$target),
               value=suppressWarnings(as.numeric(ed$value)), stringsAsFactors = FALSE)
  } else NULL
}

.sankeymatic_lines <- function(ed){
  sed <- .as_sankey_edges(ed)
  if (is.null(sed) || !nrow(sed)) return(character())
  sed$value[!is.finite(sed$value)] <- 0
  sed <- sed[sed$value > 0, , drop=FALSE]
  if (!nrow(sed)) return(character())
  sprintf("%s [%s] %s", sed$source, format(round(sed$value, 2), nsmall=2), sed$target)
}

output$sankeymatic_text <- renderText({
  dom <- input$sankey_domain %||% "Author"
  d <- .get_wos_domain(dom)
  ed <- d$edges
  lines <- .sankeymatic_lines(ed)
  if (!length(lines)) return(sprintf("No Sankey edges yet for '%s'. Run routines first.", dom))
  paste(lines, collapse = "
")
})

output$sankey_ui <- renderUI({
  dom <- input$sankey_domain %||% "Author"
  d <- .get_wos_domain(dom)
  ed <- d$edges
  sed <- .as_sankey_edges(ed)
  if (is.null(sed) || !nrow(sed)) {
    return(tags$div(class="small-note", sprintf("Sankey: no edges yet for domain '%s'. Run routines first.", dom)))
  }
  minw <- suppressWarnings(as.numeric(input$sankey_min_weight %||% 1))
  if (is.finite(minw)) sed <- sed[is.finite(sed$value) & sed$value >= minw, , drop=FALSE]
  if (!nrow(sed)) return(tags$div(class="small-note", "Sankey: no edges after filtering by Min weight."))

  if (!requireNamespace("networkD3", quietly = TRUE)) {
    return(tags$div(class="small-note", "networkD3 not installed; in-app Sankey disabled. Use SankeyMATIC text below."))
  }

  nodes <- data.frame(name = unique(c(sed$source, sed$target)), stringsAsFactors = FALSE)
  idx <- setNames(seq_len(nrow(nodes)) - 1L, nodes$name)
  links <- data.frame(source = unname(idx[sed$source]),
                      target = unname(idx[sed$target]),
                      value  = sed$value, stringsAsFactors = FALSE)

  
  # map node -> cluster (carac) to keep colors consistent with Network/Kano
  if (!is.null(d$nodes) && is.data.frame(d$nodes) && nrow(d$nodes)) {
    nd2 <- d$nodes
    if (!("name" %in% names(nd2))) {
      if ("label" %in% names(nd2)) nd2$name <- as.character(nd2$label)
      else if ("id" %in% names(nd2)) nd2$name <- as.character(nd2$id)
    }
    if (!("carac" %in% names(nd2)) && "cluster" %in% names(nd2)) nd2$carac <- nd2$cluster
    if (!("carac" %in% names(nd2))) nd2$carac <- "C1"
    nd2$name <- gsub("<br>.*$", "", as.character(nd2$name))
    nd2$carac <- as.character(nd2$carac)
    grp_map <- setNames(nd2$carac, nd2$name)
    nodes$group <- unname(grp_map[nodes$name])
    nodes$group[is.na(nodes$group) | !nzchar(nodes$group)] <- "C1"
  } else {
    nodes$group <- "C1"
  }
  pal <- .cluster_palette_map(nodes$group)
  colourScale <- networkD3::JS(sprintf(
    "d3.scaleOrdinal().domain([%s]).range([%s])",
    paste(sprintf("'%s'", names(pal)), collapse = ","),
    paste(sprintf("'%s'", unname(pal)), collapse = ",")
  ))
networkD3::sankeyNetwork(Links = links, Nodes = nodes,
                           Source = "source", Target = "target",
                           Value = "value", NodeID = "name", NodeGroup = "group",
                           colourScale = colourScale,
                           fontSize = 12, nodeWidth = 20)
})

output$ssplot_panel <- renderPlot({

  dom <- input$ss_domain %||% "Author"
  d <- .get_wos_domain(dom)
  nd <- d$nodes

  if (is.null(nd) || !is.data.frame(nd) || nrow(nd) == 0) {
    plot.new()
    text(0.5, 0.5, paste0("SSplot: domain '", dom, "' not ready (run routines first)."), cex = 0.95)
    return(invisible(NULL))
  }

  sil_df <- nd
  if (!("name" %in% names(nd))) nd$name <- if ("label" %in% names(nd)) as.character(nd$label) else as.character(nd[[1]])
  nd$name <- trimws(as.character(nd$name))
  nd$name <- sub("#[0-9]+$", "", nd$name)
  nd$name <- sub("#C[0-9]+$", "", nd$name)
  nd$label <- nd$name
  if ("value" %in% names(nd)) {
    nd$value <- suppressWarnings(as.numeric(nd$value))
    nd$value[!is.finite(nd$value) | nd$value < 1] <- 1
  }
  if ("value2" %in% names(nd)) {
    nd$value2 <- suppressWarnings(as.numeric(nd$value2))
    nd$value2[!is.finite(nd$value2)] <- 0
  }

  # Required columns for renderSSplot.R
  if (!("name" %in% names(sil_df))) {
    if ("label" %in% names(sil_df)) sil_df$name <- as.character(sil_df$label) else sil_df$name <- as.character(sil_df[[1]])
  }
  sil_df$name <- trimws(as.character(sil_df$name))
  sil_df$name <- sub("#[0-9]+$", "", sil_df$name)
  sil_df$name <- sub("#C[0-9]+$", "", sil_df$name)
  sil_df$label <- sil_df$name
  if (!("carac" %in% names(sil_df))) sil_df$carac <- "C1"

  # sil_width from ssi/SSi
  if (!("sil_width" %in% names(sil_df))) {
    if ("ssi" %in% names(sil_df)) sil_df$sil_width <- suppressWarnings(as.numeric(sil_df$ssi))
    else if ("SSi" %in% names(sil_df)) sil_df$sil_width <- suppressWarnings(as.numeric(sil_df$SSi))
    else sil_df$sil_width <- 0
  }
  sil_df$sil_width <- suppressWarnings(as.numeric(sil_df$sil_width))
  sil_df$sil_width[!is.finite(sil_df$sil_width)] <- 0

  # renderSSplot.R expects numeric cluster ids, not "C#" strings
  .normClusterNum <- function(x){
    suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(x))))
  }
  sil_df$carac <- .normClusterNum(sil_df$carac)
  sil_df$carac[is.na(sil_df$carac)] <- 1L

  # preserve nearest-cluster ids if they already exist on nodes; otherwise fall back to own cluster
  if (!("neighborC" %in% names(sil_df)) || all(is.na(.normClusterNum(sil_df$neighborC)))) {
    if ("neighborC" %in% names(nd)) sil_df$neighborC <- nd$neighborC else sil_df$neighborC <- sil_df$carac
  }
  sil_df$neighborC <- .normClusterNum(sil_df$neighborC)
  sil_df$neighborC[is.na(sil_df$neighborC)] <- sil_df$carac[is.na(sil_df$neighborC)]

  # Force SSplot inputs to metatable document counts so SSplot n matches Summary counts
  mt_counts <- tryCatch(.domain_doc_counts_from_metatable(rv$metatable, dom), error = function(e) setNames(numeric(), character()))
  if (length(mt_counts)) {
    nd_name_clean0 <- trimws(sub("#[0-9]+$", "", sub("#C[0-9]+$", "", as.character(nd$name))))
    hit0 <- unname(mt_counts[nd_name_clean0])
    hit0[!is.finite(hit0)] <- NA_real_
    if (length(hit0) == nrow(nd)) nd$value <- ifelse(is.finite(hit0) & hit0 > 0, hit0, nd$value)
  }

  # Force SSplot inputs to app-defined values: value=document count, value2=sum(edge)
  if (is.data.frame(nd) && nrow(nd)) {
    nd_name_clean <- trimws(sub("#[0-9]+$", "", sub("#C[0-9]+$", "", as.character(nd$name))))
    sil_name_clean <- trimws(sub("#[0-9]+$", "", sub("#C[0-9]+$", "", as.character(sil_df$name))))
    ixv <- match(sil_name_clean, nd_name_clean)
    if ("value" %in% names(nd)) sil_df$value <- suppressWarnings(as.numeric(nd$value[ixv]))
    if ("value2" %in% names(nd)) sil_df$value2 <- suppressWarnings(as.numeric(nd$value2[ixv]))
    if ("carac" %in% names(nd)) sil_df$carac <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(nd$carac[ixv]))))
  }
  sil_df$value[!is.finite(sil_df$value) | sil_df$value < 1] <- 1
  sil_df$value2[!is.finite(sil_df$value2)] <- 0
  if ("value" %in% names(nd)) {
    nd$value <- suppressWarnings(as.numeric(nd$value))
    nd$value[!is.finite(nd$value) | nd$value < 1] <- 1
  }
  if ("value2" %in% names(nd)) {
    nd$value2 <- suppressWarnings(as.numeric(nd$value2))
    nd$value2[!is.finite(nd$value2)] <- 0
  }

  # clusterwise summary for SSplot; also expose overall totals for renderers that use them
  clv <- unique(stats::na.omit(sil_df$carac))
  results <- do.call(rbind, lapply(clv, function(cc){
    sub <- sil_df[sil_df$carac == cc, , drop = FALSE]
    qwc <- suppressWarnings(as.numeric(sub$Qw_cluster %||% sub$Qw %||% NA_real_))
    quc <- suppressWarnings(as.numeric(sub$Qu_cluster %||% sub$Qu %||% NA_real_))
    data.frame(
      Cluster = paste0("C", cc),
      SS = mean(sub$sil_width, na.rm = TRUE),
      Qw = if (all(!is.finite(qwc))) 0 else mean(qwc[is.finite(qwc)], na.rm = TRUE),
      Qu = if (all(!is.finite(quc))) 0 else mean(quc[is.finite(quc)], na.rm = TRUE),
      n = sum(suppressWarnings(as.numeric(sub$value)), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  results$SS[!is.finite(results$SS)] <- 0
  results$Qw[!is.finite(results$Qw)] <- 0
  results$Qu[!is.finite(results$Qu)] <- 0
  results$SS_total <- mean(sil_df$sil_width, na.rm = TRUE)
  results$Qw_total <- sum(results$Qw, na.rm = TRUE)
  results$Qu_total <- sum(results$Qu, na.rm = TRUE)
  results <- rbind(
    results,
    data.frame(
      Cluster = "OVERALL",
      SS = results$SS_total[1],
      Qw = results$Qw_total[1],
      Qu = results$Qu_total[1],
      n = sum(suppressWarnings(as.numeric(sil_df$value)), na.rm = TRUE),
      SS_total = results$SS_total[1],
      Qw_total = results$Qw_total[1],
      Qu_total = results$Qu_total[1],
      stringsAsFactors = FALSE
    )
  )

  if (exists("render_panel", mode = "function")) {
    tryCatch(
      render_panel(
        sil_df = sil_df,
        nodes0 = nd,
        results = results,
        nodes = nd,
        font_scale = (input$ss_font_scale %||% 1.3)
      ),
      error = function(e){
        plot.new()
        text(0.5, 0.5, paste("SSplot render error:", e$message), cex = 1.0)
      }
    )
  } else if (exists(".render_panel_ss", mode = "function")) {
    tryCatch(
      .render_panel_ss(
        sil_df = sil_df,
        results = results,
        nodes0 = nd,
        nodes = nd,
        font_scale = (input$ss_font_scale %||% 1.3)
      ),
      error = function(e){
        plot.new()
        text(0.5, 0.5, paste("SSplot render error:", e$message), cex = 1.0)
      }
    )
  } else {
    plot.new()
    text(0.5, 0.5, "SSplot: renderSSplot.R not loaded (missing render_panel())", cex = 0.95)
  }

}, height = 820)


output$ss_kano_plot <- renderPlot({
  dom <- input$ss_kano_domain %||% "Author"
  d <- .get_wos_domain(dom)
  nd <- d$nodes
  validate(need(is.data.frame(nd) && nrow(nd) > 1, "Run routines first (need Top20 nodes)."))

  # normalize SS / a*
  if (!("ssi" %in% names(nd)) && ("SSi" %in% names(nd))) nd$ssi <- nd$SSi
  if (!("a_star1" %in% names(nd)) && ("a_star" %in% names(nd))) nd$a_star1 <- nd$a_star
  if (!("a_star1" %in% names(nd)) && ("a_i" %in% names(nd))) nd$a_star1 <- 1/(1+as.numeric(nd$a_i))

  if (exists("kano_plot_ss_astar", mode="function")) {
    if (!("carac" %in% names(nd))) nd$carac <- 1L
    p <- tryCatch(kano_plot_ss_astar(nd, d$edges, label_size = input$ss_kano_label_size,
                                    title_txt = sprintf("SS Kano (%s): SS vs a*", dom)),
                  error=function(e) NULL)
    if (!is.null(p)) { try(print(p), silent=TRUE); return(invisible()) }
  }

  .kano_scatter(nd, xcol = "a_star1", ycol = "ssi",
               title_txt = sprintf("SS Kano (%s): SS vs a*", dom),
               label_size = input$ss_kano_label_size, xlab="a*", ylab="SS")
}, height = 700)


output$taaa_table <- renderTable({
  validate(need(!is.null(rv$taaa_df) && is.data.frame(rv$taaa_df) && nrow(rv$taaa_df) > 0,
                "Run WoS (Upload) first to generate the TAAA table."))
  rv$taaa_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$taaa_theme_freq_table <- renderTable({
  validate(need(!is.null(rv$taaa_freq_df) && is.data.frame(rv$taaa_freq_df) && nrow(rv$taaa_freq_df) > 0,
                "Run WoS (Upload) first to generate the TAAA theme-frequency table."))
  rv$taaa_freq_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$taaa_profile_map_table <- renderTable({
  validate(need(!is.null(rv$taaa_profile_map_df) && is.data.frame(rv$taaa_profile_map_df) && nrow(rv$taaa_profile_map_df) > 0,
                "Shown when a Profile column exists and the number of themes equals max(Profile) after Hungarian mapping."))
  rv$taaa_profile_map_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$taaa_kappa_table <- renderTable({
  validate(need(!is.null(rv$taaa_kappa_df) && is.data.frame(rv$taaa_kappa_df) && nrow(rv$taaa_kappa_df) > 0,
                "Shown when a Profile column exists and the number of themes equals max(Profile) after Hungarian mapping."))
  rv$taaa_kappa_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$taaa_confusion_table <- renderTable({
  validate(need(!is.null(rv$taaa_conf_df) && is.data.frame(rv$taaa_conf_df) && nrow(rv$taaa_conf_df) > 0,
                "Shown when a Profile column exists and the number of themes equals max(Profile) after Hungarian mapping."))
  rv$taaa_conf_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$lotka_plot <- renderPlot({
  res <- lotka_res_reactive()
  validate(need(!is.null(res), "Run routines first to generate the Lotka tab."))
  .plot_lotka_result(res)
}, height = 460)

output$lotka_test_table <- renderTable({
  res <- lotka_res_reactive()
  validate(need(!is.null(res), "Run routines first to generate the Lotka tab."))
  .lotka_summary_table(res)
}, striped = TRUE, bordered = TRUE, spacing = "s", width = "100%")

output$lotka_observed_expected_table <- renderTable({
  res <- lotka_res_reactive()
  validate(need(!is.null(res), "Run routines first to generate the Lotka tab."))
  df <- as.data.frame(res$observed_expected_table, stringsAsFactors = FALSE)
  if ("expected_raw" %in% names(df)) df$expected_raw <- round(as.numeric(df$expected_raw), 4)
  if ("expected" %in% names(df)) df$expected <- round(as.numeric(df$expected), 4)
  if ("residual" %in% names(df)) df$residual <- round(as.numeric(df$residual), 4)
  df
}, striped = TRUE, bordered = TRUE, spacing = "s", width = "100%")

output$lotka_chisq_table <- renderTable({
  res <- lotka_res_reactive()
  validate(need(!is.null(res), "Run routines first to generate the Lotka tab."))
  df <- as.data.frame(res$test_table, stringsAsFactors = FALSE)
  if ("expected" %in% names(df)) df$expected <- round(as.numeric(df$expected), 4)
  df
}, striped = TRUE, bordered = TRUE, spacing = "s", width = "100%")

  # -----------------------------
  # Maps (World / China / USA) - WoS standalone
  # -----------------------------
  .mt <- function() {
    if (!is.null(rv$metatable) && is.data.frame(rv$metatable) && nrow(rv$metatable) > 0) return(rv$metatable)
    NULL
  }

  .collect_terms <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- x[nzchar(x)]
    x
  }

  .country_counts <- reactive({
    mt <- .mt()
    if (is.null(mt)) return(data.frame(country=character(), n=numeric(), stringsAsFactors = FALSE))

    cols <- intersect(c("FAcountry","CAcountry","Country","country","FA_Country","CA_Country"), names(mt))
    if (length(cols) == 0) return(data.frame(country=character(), n=numeric(), stringsAsFactors = FALSE))

    v <- unlist(lapply(cols, function(cc) .collect_terms(mt[[cc]])), use.names = FALSE)
    if (!length(v)) return(data.frame(country=character(), n=numeric(), stringsAsFactors = FALSE))

    tb <- sort(table(v), decreasing = TRUE)
    data.frame(country = names(tb), n = as.numeric(tb), stringsAsFactors = FALSE)
  })

  .split_terms_cell <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return(character(0))
    s <- trimws(as.character(x))
    if (!nzchar(s)) return(character(0))
    parts <- unlist(strsplit(s, "[;|\n]+", perl = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    unique(parts[nzchar(parts)])
  }

  .region_pairs <- reactive({
    mt <- .mt()
    if (is.null(mt) || !nrow(mt)) return(data.frame(country=character(), region=character(), stringsAsFactors = FALSE))

    pair_cols <- list(
      c("FAcountry","FAregion"),
      c("CAcountry","CAregion"),
      c("Country","Region"),
      c("country","region")
    )
    pair_cols <- pair_cols[vapply(pair_cols, function(z) all(z %in% names(mt)), logical(1))]
    if (!length(pair_cols)) return(data.frame(country=character(), region=character(), stringsAsFactors = FALSE))

    out <- vector("list", 0)
    for (pair in pair_cols) {
      ccol <- pair[1]; rcol <- pair[2]
      for (i in seq_len(nrow(mt))) {
        cs <- .split_terms_cell(mt[[ccol]][i])
        rs <- .split_terms_cell(mt[[rcol]][i])
        if (!length(cs) || !length(rs)) next
        if (length(cs) == length(rs)) {
          out[[length(out)+1]] <- data.frame(country = cs, region = rs, stringsAsFactors = FALSE)
        } else if (length(cs) == 1) {
          out[[length(out)+1]] <- data.frame(country = rep(cs, length(rs)), region = rs, stringsAsFactors = FALSE)
        } else if (length(rs) == 1) {
          out[[length(out)+1]] <- data.frame(country = cs, region = rep(rs, length(cs)), stringsAsFactors = FALSE)
        } else {
          n <- min(length(cs), length(rs))
          out[[length(out)+1]] <- data.frame(country = cs[seq_len(n)], region = rs[seq_len(n)], stringsAsFactors = FALSE)
        }
      }
    }

    if (!length(out)) return(data.frame(country=character(), region=character(), stringsAsFactors = FALSE))
    df <- do.call(rbind, out)
    df$country <- .std_country(trimws(as.character(df$country)))
    df$region  <- trimws(as.character(df$region))
    df <- df[nzchar(df$country) & nzchar(df$region), , drop=FALSE]
    if (!nrow(df)) return(data.frame(country=character(), region=character(), stringsAsFactors = FALSE))
    df
  })

  # ---- World map ----
  output$country_map <- shiny::renderPlot({
    counts <- .country_counts()
    if (!is.data.frame(counts) || nrow(counts) == 0) {
      plot.new(); text(0.5, 0.5, "World map: no country data yet.\nRun WoS routines to build metatable."); return()
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      plot.new(); text(0.5, 0.5, "World map requires ggplot2."); return()
    }
    world <- ggplot2::map_data("world")
    counts$country2 <- tolower(trimws(counts$country))
    world$region2   <- tolower(trimws(world$region))

    agg <- dplyr::as_tibble(counts) %>%
      dplyr::group_by(country2) %>% dplyr::summarise(n = sum(n, na.rm = TRUE), .groups="drop")

    world2 <- dplyr::left_join(world, agg, by = c("region2" = "country2"))
    world2$n0 <- world2$n
    world2$n0[is.na(world2$n0)] <- 0
    world2$fillv <- ifelse(world2$n0 <= 0, NA, log10(world2$n0 + 1))

    ggplot2::ggplot(world2, ggplot2::aes(long, lat, group = group)) +
      ggplot2::geom_polygon(ggplot2::aes(fill = fillv), color = "white", linewidth = 0.1) +
      ggplot2::coord_fixed(1.3) +
      ggplot2::theme_void() +
      ggplot2::labs(fill = "Count (log10)") +
      ggplot2::scale_fill_gradient(na.value = "grey95")
  })

  # ---- USA map ----
  usa_state_counts <- reactive({
    rp <- .region_pairs()
    if (!nrow(rp)) return(data.frame(state=character(), code=character(), count=numeric(), stringsAsFactors = FALSE))

    usa_keys <- c("usa","us","u.s.","u.s.a.","united states","united states of america")
    rp$country2 <- tolower(trimws(rp$country))
    sp <- rp[rp$country2 %in% usa_keys, , drop=FALSE]
    if (!nrow(sp)) return(data.frame(state=character(), code=character(), count=numeric(), stringsAsFactors = FALSE))

    st_raw <- trimws(as.character(sp$region))
    st_raw <- st_raw[nzchar(st_raw)]
    if (!length(st_raw)) return(data.frame(state=character(), code=character(), count=numeric(), stringsAsFactors = FALSE))

    name_to_code <- setNames(state.abb, tolower(state.name))
    code <- toupper(st_raw)
    code <- ifelse(nchar(code) == 2, code, name_to_code[tolower(st_raw)])
    code[is.na(code)] <- ""
    code <- code[nchar(code) == 2]
    if (!length(code)) return(data.frame(state=character(), code=character(), count=numeric(), stringsAsFactors = FALSE))

    tb <- sort(table(code), decreasing = TRUE)
    data.frame(
      state = state.name[match(names(tb), state.abb)],
      code  = names(tb),
      count = as.numeric(tb),
      stringsAsFactors = FALSE
    )
  })

  output$usa_country_map <- renderPlot({
    df <- usa_state_counts()
    if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("maps", quietly = TRUE)) {
      plot.new(); text(0.5, 0.5, "USA map requires packages: ggplot2 + maps", cex=1.1); return(invisible())
    }
    usa_poly <- ggplot2::map_data("state")
    if (!is.data.frame(df) || !nrow(df)) {
      gg <- ggplot2::ggplot(usa_poly, ggplot2::aes(long, lat, group = group)) +
        ggplot2::geom_polygon(fill = "grey90", color = "white", linewidth = 0.2) +
        ggplot2::coord_fixed(1.3) + ggplot2::theme_void() + ggplot2::ggtitle("USA (States) - no state data yet")
      print(gg); return(invisible())
    }
    df$region <- tolower(df$state)
    usa2 <- dplyr::left_join(usa_poly, df, by = "region")
    usa2$fillv <- ifelse(is.na(usa2$count) | usa2$count <= 0, NA_real_, log10(usa2$count + 1))
    gg <- ggplot2::ggplot(usa2, ggplot2::aes(long, lat, group = group)) +
      ggplot2::geom_polygon(ggplot2::aes(fill = fillv), color = "white", linewidth = 0.2) +
      ggplot2::coord_fixed(1.3) + ggplot2::theme_void() +
      ggplot2::scale_fill_gradient(na.value = "grey95") +
      ggplot2::labs(fill = "Count (log10)") +
      ggplot2::ggtitle(sprintf("USA (States) - pre-FLCA | Total=%s", sum(df$count, na.rm=TRUE)))
    print(gg)
  })

  output$tbl_usa_counts <- DT::renderDT({
    df <- usa_state_counts()
    if (!is.data.frame(df) || nrow(df) == 0) {
      return(DT::datatable(data.frame(Message="No USA state terms found in FAregion/CAregion.", stringsAsFactors=FALSE),
                           options=list(dom='t'), rownames=FALSE))
    }
    DT::datatable(df[, c("state","code","count")], options=list(pageLength=20), rownames=FALSE)
  })

  # ---- China map ----
  # Normalize region keys to hchinamap's expected Chinese province/SAR names.
  .to_hchinamap_name <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x)] <- ""
    if (!length(x)) return(x)

    norm <- tolower(x)
    norm <- gsub("[[:space:][:punct:]]+", " ", norm, perl = TRUE)
    norm <- trimws(norm)

    # direct alias map -> hchinamap join keys
    alias <- c(
      "beijing"="北京", "peking"="北京", "北京"="北京", "北京市"="北京",
      "tianjin"="天津", "天津"="天津", "天津市"="天津",
      "hebei"="河北", "河北"="河北", "河北省"="河北",
      "shanxi"="山西", "山西"="山西", "山西省"="山西",
      "inner mongolia"="内蒙古", "nei mongol"="内蒙古", "inner mongol"="内蒙古", "内蒙古"="内蒙古", "內蒙古"="内蒙古", "内蒙古自治区"="内蒙古", "內蒙古自治區"="内蒙古",
      "liaoning"="辽宁", "辽宁"="辽宁", "遼寧"="辽宁", "辽宁省"="辽宁", "遼寧省"="辽宁",
      "jilin"="吉林", "吉林"="吉林", "吉林省"="吉林",
      "heilongjiang"="黑龙江", "黑龙江"="黑龙江", "黑龍江"="黑龙江", "黑龙江省"="黑龙江", "黑龍江省"="黑龙江",
      "shanghai"="上海", "上海"="上海", "上海市"="上海",
      "jiangsu"="江苏", "江苏"="江苏", "江蘇"="江苏", "江苏省"="江苏", "江蘇省"="江苏",
      "zhejiang"="浙江", "浙江"="浙江", "浙江省"="浙江",
      "anhui"="安徽", "安徽"="安徽", "安徽省"="安徽",
      "fujian"="福建", "福建"="福建", "福建省"="福建",
      "jiangxi"="江西", "江西"="江西", "江西省"="江西",
      "shandong"="山东", "山东"="山东", "山東"="山东", "山东省"="山东", "山東省"="山东",
      "henan"="河南", "河南"="河南", "河南省"="河南",
      "hubei"="湖北", "湖北"="湖北", "湖北省"="湖北",
      "hunan"="湖南", "湖南"="湖南", "湖南省"="湖南",
      "guangdong"="广东", "canton"="广东", "广东"="广东", "廣東"="广东", "广东省"="广东", "廣東省"="广东",
      "guangxi"="广西", "广西"="广西", "廣西"="广西", "guangxi zhuang autonomous region"="广西", "广西壮族自治区"="广西", "廣西壯族自治區"="广西",
      "hainan"="海南", "海南"="海南", "海南省"="海南",
      "chongqing"="重庆", "重庆"="重庆", "重慶"="重庆", "重庆市"="重庆", "重慶市"="重庆",
      "sichuan"="四川", "四川"="四川", "四川省"="四川",
      "guizhou"="贵州", "贵州"="贵州", "貴州"="贵州", "贵州省"="贵州", "貴州省"="贵州",
      "yunnan"="云南", "云南"="云南", "雲南"="云南", "云南省"="云南", "雲南省"="云南",
      "tibet"="西藏", "xizang"="西藏", "西藏"="西藏", "西藏自治区"="西藏", "西藏自治區"="西藏",
      "shaanxi"="陕西", "陕西"="陕西", "陝西"="陕西", "陕西省"="陕西", "陝西省"="陕西",
      "gansu"="甘肃", "甘肃"="甘肃", "甘肅"="甘肃", "甘肃省"="甘肃", "甘肅省"="甘肃",
      "qinghai"="青海", "青海"="青海", "青海省"="青海",
      "ningxia"="宁夏", "宁夏"="宁夏", "寧夏"="宁夏", "ningxia hui autonomous region"="宁夏", "宁夏回族自治区"="宁夏", "寧夏回族自治區"="宁夏",
      "xinjiang"="新疆", "新疆"="新疆", "xinjiang uygur autonomous region"="新疆", "xinjiang uygur autonomous region"="新疆", "新疆维吾尔自治区"="新疆", "新疆維吾爾自治區"="新疆",
      "taiwan"="台湾", "taiwan province"="台湾", "台灣"="台湾", "臺灣"="台湾", "台湾"="台湾", "臺北"="台湾", "台北"="台湾", "taipei"="台湾", "new taipei"="台湾", "taichung"="台湾", "tainan"="台湾", "kaohsiung"="台湾",
      "hong kong"="香港", "hong kong sar"="香港", "hongkong"="香港", "香港"="香港", "香港特别行政区"="香港", "香港特別行政區"="香港",
      "macau"="澳门", "macao"="澳门", "macau sar"="澳门", "澳門"="澳门", "澳门"="澳门", "澳门特别行政区"="澳门", "澳門特別行政區"="澳门"
    )

    out <- unname(alias[norm])

    # fallback: strip common Chinese suffixes and try again
    need <- is.na(out) | !nzchar(out)
    if (any(need)) {
      raw <- x[need]
      raw2 <- gsub("(省|市|特别行政区|特別行政區|壮族自治区|壯族自治區|回族自治区|回族自治區|维吾尔自治区|維吾爾自治區|自治区|自治區)$", "", raw, perl = TRUE)
      raw2 <- trimws(raw2)
      norm2 <- tolower(gsub("[[:space:][:punct:]]+", " ", raw2, perl = TRUE))
      out2 <- unname(alias[norm2])
      out2[is.na(out2) | !nzchar(out2)] <- raw2[is.na(out2) | !nzchar(out2)]
      out[need] <- out2
    }

    # final cleanup to hchinamap province names
    out <- gsub("市$", "", out, perl = TRUE)
    out <- gsub("省$", "", out, perl = TRUE)
    out <- gsub("特别行政区$|特別行政區$", "", out, perl = TRUE)
    out <- gsub("壮族自治区$|壯族自治區$|回族自治区$|回族自治區$|维吾尔自治区$|維吾爾自治區$|自治区$|自治區$", "", out, perl = TRUE)
    out <- trimws(out)
    out
  }

  china_counts <- reactive({
    rp <- .region_pairs()
    if (!nrow(rp)) return(data.frame(name_cn=character(), value=numeric(), stringsAsFactors = FALSE))
    rp$country2 <- tolower(trimws(rp$country))
    china_keys <- c("china","peoples r china","people's r china","people s r china","people's republic of china","people s republic of china","pr china","p.r. china","prc","taiwan","hong kong","hongkong","macau","macao")
    df <- rp[rp$country2 %in% china_keys, , drop=FALSE]
    if (!nrow(df)) return(data.frame(name_cn=character(), value=numeric(), stringsAsFactors = FALSE))

    reg <- trimws(as.character(df$region))
    reg <- reg[nzchar(reg)]
    if (!length(reg)) return(data.frame(name_cn=character(), value=numeric(), stringsAsFactors = FALSE))

    mapped <- .to_hchinamap_name(reg)
    mapped <- mapped[nzchar(mapped)]
    if (!length(mapped)) return(data.frame(name_cn=character(), value=numeric(), stringsAsFactors = FALSE))

    valid_names <- c("北京","天津","河北","山西","内蒙古","辽宁","吉林","黑龙江","上海","江苏","浙江","安徽","福建","江西","山东","河南","湖北","湖南","广东","广西","海南","重庆","四川","贵州","云南","西藏","陕西","甘肃","青海","宁夏","新疆","台湾","香港","澳门")
    mapped <- mapped[mapped %in% valid_names]
    if (!length(mapped)) return(data.frame(name_cn=character(), value=numeric(), stringsAsFactors = FALSE))

    tb <- sort(table(mapped), decreasing = TRUE)
    data.frame(name_cn = names(tb), value = as.numeric(tb), stringsAsFactors = FALSE)
  })

  output$china_map <- shiny::renderUI({
    df <- china_counts()

    if (!requireNamespace("hchinamap", quietly = TRUE)) {
      return(shiny::tags$div(class = "alert alert-warning",
                             "China map requires package: hchinamap. Please install.packages('hchinamap')."))
    }
    if (!is.data.frame(df) || nrow(df) == 0) {
      return(shiny::tags$div(class = "alert alert-info",
                             "China map: no China province/SAR terms found after hchinamap name normalization."))
    }

    hchinamap::hchinamap(
      name   = as.character(df$name_cn),
      value  = as.numeric(df$value),
      region = "China",
      width  = "100%",
      height = "650px",
      title  = "Map of China"
    )
  })

  output$tbl_china_counts <- DT::renderDT({
    df <- china_counts()
    if (!is.data.frame(df) || nrow(df)==0) df <- data.frame(region=character(), count=numeric())
    else df <- df %>% dplyr::transmute(region = name_cn, count = value) %>% dplyr::arrange(dplyr::desc(count))
    DT::datatable(df, options=list(pageLength=10), rownames=FALSE)
  })


}

shinyApp(ui = ui, server = server)
