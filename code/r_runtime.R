ensure_user_library <- function() {
  r_minor <- strsplit(R.version$minor, "[.]", fixed = FALSE)[[1]][[1]]
  r_libs_user <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  candidate_libs <- character()
  if (!is.na(r_libs_user) && nzchar(r_libs_user)) {
    r_libs_user <- gsub("\\\\", "/", r_libs_user)
    r_version <- paste0(R.version$major, ".", r_minor)
    r_libs_user <- gsub("%v", r_version, r_libs_user, fixed = TRUE)
    r_libs_user <- gsub("%V", paste0(R.version$major, ".", R.version$minor), r_libs_user, fixed = TRUE)
    r_libs_user <- gsub("%p", R.version$platform, r_libs_user, fixed = TRUE)
    candidate_libs <- c(candidate_libs, r_libs_user)
  }

  local_appdata <- Sys.getenv("LOCALAPPDATA", unset = NA_character_)
  if (is.na(local_appdata) || !nzchar(local_appdata)) {
    userprofile <- Sys.getenv("USERPROFILE", unset = NA_character_)
    if (is.na(userprofile) || !nzchar(userprofile)) userprofile <- path.expand("~")
    local_appdata <- file.path(userprofile, "AppData", "Local")
  }

  user_lib <- paste0(gsub("\\\\", "/", local_appdata), "/R/win-library/", R.version$major, ".", r_minor)
  userprofile <- Sys.getenv("USERPROFILE", unset = NA_character_)
  if (!is.na(userprofile) && nzchar(userprofile)) {
    candidate_libs <- c(candidate_libs, paste0(gsub("\\\\", "/", userprofile), "/AppData/Local/R/win-library/", R.version$major, ".", r_minor))
  }
  candidate_libs <- c(candidate_libs, user_lib)
  candidate_libs <- unique(candidate_libs[dir.exists(candidate_libs)])
  .libPaths(c(candidate_libs, .libPaths()))

  invisible(.libPaths())
}
