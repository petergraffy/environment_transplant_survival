resolve_saf_dir <- function(required = TRUE, release = Sys.getenv("SAF_RELEASE", "q2_2025")) {
  release <- tolower(release)

  if (release %in% c("q1_2026", "2026q1", "2603", "most_recent")) {
    env_candidates <- c(
      Sys.getenv("SAF_Q1_2026_DIR", unset = NA_character_),
      Sys.getenv("SAF_2603_DIR", unset = NA_character_),
      Sys.getenv("SAF_DIR", unset = NA_character_)
    )

    userprofile <- Sys.getenv("USERPROFILE", unset = NA_character_)
    home <- Sys.getenv("HOME", unset = NA_character_)

    path_candidates <- c(
      env_candidates,
      if (!is.na(userprofile)) file.path(userprofile, "Box", "Q1 2026 SAF"),
      if (!is.na(home)) file.path(home, "Box", "Q1 2026 SAF"),
      if (!is.na(home)) file.path(home, "Library", "CloudStorage", "Box-Box", "Q1 2026 SAF")
    )

    path_candidates <- unique(path_candidates[!is.na(path_candidates) & nzchar(path_candidates)])
    found <- path_candidates[dir.exists(path_candidates)]

    if (length(found) > 0L) {
      return(normalizePath(found[[1]], winslash = "/", mustWork = TRUE))
    }

    if (isTRUE(required)) {
      stop(
        "Could not find SAF Q1 2026 data directory. Set SAF_Q1_2026_DIR to the folder ",
        "containing pubsaf2603 and Supplement2603. Checked: ",
        paste(path_candidates, collapse = "; "),
        call. = FALSE
      )
    }

    return(NA_character_)
  }

  env_candidates <- c(
    Sys.getenv("SAF_Q2_2025_DIR", unset = NA_character_),
    Sys.getenv("SAF_DIR", unset = NA_character_)
  )

  userprofile <- Sys.getenv("USERPROFILE", unset = NA_character_)
  home <- Sys.getenv("HOME", unset = NA_character_)

  path_candidates <- c(
    env_candidates,
    if (!is.na(userprofile)) file.path(userprofile, "Box", "SAF Q2 2025"),
    if (!is.na(home)) file.path(home, "Box", "SAF Q2 2025"),
    if (!is.na(home)) file.path(home, "Library", "CloudStorage", "Box-Box", "SAF Q2 2025")
  )

  path_candidates <- unique(path_candidates[!is.na(path_candidates) & nzchar(path_candidates)])
  found <- path_candidates[dir.exists(path_candidates)]

  if (length(found) > 0L) {
    return(normalizePath(found[[1]], winslash = "/", mustWork = TRUE))
  }

  if (isTRUE(required)) {
    stop(
      "Could not find SAF Q2 2025 data directory. Set SAF_Q2_2025_DIR to the folder ",
      "containing pubsaf2506 and SupplementalData2506. Checked: ",
      paste(path_candidates, collapse = "; "),
      call. = FALSE
    )
  }

  NA_character_
}

get_saf_paths <- function(required = TRUE, release = Sys.getenv("SAF_RELEASE", "q2_2025")) {
  release <- tolower(release)
  saf_dir <- resolve_saf_dir(required = required, release = release)

  if (release %in% c("q1_2026", "2026q1", "2603", "most_recent")) {
    paths <- list(
      release = "q1_2026",
      saf_dir = saf_dir,
      pubsaf_dir = file.path(saf_dir, "pubsaf2603"),
      supp_dir = file.path(saf_dir, "Supplement2603"),
      canzip_file = file.path(saf_dir, "Supplement2603", "canzip2603.sas7bdat"),
      data_dictionary = file.path(saf_dir, "dataDictionary.html"),
      linking_diagram = file.path(saf_dir, "SAFsLinkingDiagram.pdf")
    )

    if (isTRUE(required)) {
      required_dirs <- c(paths$pubsaf_dir, paths$supp_dir)
      missing_dirs <- required_dirs[!dir.exists(required_dirs)]
      if (length(missing_dirs) > 0L) {
        stop("SAF directory is missing required subdirectories: ", paste(missing_dirs, collapse = "; "), call. = FALSE)
      }
    }

    return(paths)
  }

  paths <- list(
    release = "q2_2025",
    saf_dir = saf_dir,
    pubsaf_dir = file.path(saf_dir, "pubsaf2506"),
    supp_dir = file.path(saf_dir, "SupplementalData2506"),
    canzip_file = file.path(saf_dir, "SupplementalData2506", "canzip2506.sas7bdat"),
    data_dictionary = file.path(saf_dir, "dataDictionary.html"),
    linking_diagram = file.path(saf_dir, "SAFsLinkingDiagram.pdf")
  )

  if (isTRUE(required)) {
    required_dirs <- c(paths$pubsaf_dir, paths$supp_dir)
    missing_dirs <- required_dirs[!dir.exists(required_dirs)]
    if (length(missing_dirs) > 0L) {
      stop("SAF directory is missing required subdirectories: ", paste(missing_dirs, collapse = "; "), call. = FALSE)
    }
  }

  paths
}

assert_saf_files <- function(paths, include_stathist = FALSE) {
  required_files <- c(
    file.path(paths$pubsaf_dir, c("cand_kipa.sas7bdat", "cand_liin.sas7bdat", "cand_thor.sas7bdat")),
    paths$canzip_file
  )

  if (isTRUE(include_stathist)) {
    required_files <- c(
      required_files,
      file.path(paths$pubsaf_dir, c("stathist_kipa.sas7bdat", "stathist_liin.sas7bdat", "stathist_thor.sas7bdat"))
    )
  }

  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("SAF directory is missing required files: ", paste(missing_files, collapse = "; "), call. = FALSE)
  }

  invisible(required_files)
}
