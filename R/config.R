#' Resolve a path from a config file
#'
#' Relative paths are taken relative to the project root; absolute paths
#' (e.g., a scratch directory set in a machine-specific config section) are
#' returned as is.
#'
#' @param path character path
#'
#' @returns absolute path

resolve_path <- function(path) {
  if (grepl("^(/|~|[A-Za-z]:[/\\\\])", path)) path.expand(path) else here::here(path)
}

#' Path to a canonical data source
#'
#' Looks up `name` in `config/paths.yml` with `config::get()`, so the active
#' section (environment variable `R_CONFIG_ACTIVE`, defaulting to `default`)
#' is respected, then fills any `{placeholder}` with the matching argument.
#'
#' @param name name of the data source, e.g., "stimuli"
#' @param ... values for placeholders, e.g., `region = "V1"`
#' @param file path to the paths file, defaults to `config/paths.yml`
#' @param config name of the config section, defaults to the active one
#'
#' @returns absolute path to the data source
#'
#' @examples
#' data_path("mua", region = "IT")
#' data_path("mua", region = "IT", config = "server")

data_path <- function(
  name,
  ...,
  file = here::here("config", "paths.yml"),
  config = Sys.getenv("R_CONFIG_ACTIVE", "default")
) {
  paths <- config::get(file = file, config = config)
  if (!name %in% names(paths))
    stop("Unknown data source '", name, "'. See ", file, ".")

  path <- paths[[name]]
  fills <- list(...)
  for (key in names(fills))
    path <- gsub(paste0("{", key, "}"), fills[[key]], path, fixed = T)
  if (grepl("\\{[^}]+\\}", path))
    stop("Unfilled placeholder in '", path, "'.")

  resolve_path(path)
}

#' Write a provenance record next to an output
#'
#' Saves `<output>.json` with the git commit, whether the working tree had
#' uncommitted changes, the script, its options, the input paths, the active
#' config section, and the R version, so that an output can be traced back to the exact code and
#' settings that produced it.
#'
#' @param output path of the output file being described
#' @param script path of the script that produced it
#' @param opts list of options the script ran with
#' @param inputs named list or vector of input paths
#'
#' @returns path of the provenance file, invisibly

write_provenance <- function(output, script, opts, inputs) {
  git <- function(...) {
    res <- tryCatch(
      system2("git", c("-C", here::here(), ...), stdout = T, stderr = F),
      error = function(e) NA_character_,
      warning = function(w) NA_character_)
    if (length(res) == 0) "" else res
  }
  rel <- function(p) sub(paste0("^", here::here(), "/?"), "", p)

  record <- list(
    output = rel(output),
    script = rel(script),
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    git_commit = git("rev-parse", "HEAD")[1],
    git_dirty = any(nzchar(git("status", "--porcelain", "--", "R", "scripts", "config"))),
    host = Sys.info()[["nodename"]],
    config_active = Sys.getenv("R_CONFIG_ACTIVE", "default"),
    r_version = R.version.string,
    options = opts,
    inputs = lapply(as.list(inputs), rel))

  f <- paste0(output, ".json")
  jsonlite::write_json(record, f, auto_unbox = T, pretty = T)
  invisible(f)
}
