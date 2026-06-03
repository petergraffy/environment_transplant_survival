ensure_user_library <- function() {
  r_libs_user <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  if (!is.na(r_libs_user) && nzchar(r_libs_user)) {
    r_libs_user <- gsub("\\\\", "/", r_libs_user)
    .libPaths(c(r_libs_user, .libPaths()))
    return(invisible(.libPaths()))
  }

  local_appdata <- Sys.getenv("LOCALAPPDATA", unset = NA_character_)
  if (is.na(local_appdata) || !nzchar(local_appdata)) {
    userprofile <- Sys.getenv("USERPROFILE", unset = NA_character_)
    if (is.na(userprofile) || !nzchar(userprofile)) userprofile <- path.expand("~")
    local_appdata <- file.path(userprofile, "AppData", "Local")
  }

  r_minor <- strsplit(R.version$minor, "[.]", fixed = FALSE)[[1]][[1]]
  user_lib <- paste0(gsub("\\\\", "/", local_appdata), "/R/win-library/", R.version$major, ".", r_minor)
  .libPaths(c(user_lib, .libPaths()))

  invisible(.libPaths())
}
