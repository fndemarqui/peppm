




#usethis::use_travis()
devtools::install(quick = TRUE)
devtools::document()

usethis::use_build_ignore("inst/script_development.R")
usethis::use_build_ignore(".travis.yml")
# usethis::use_cran_comments()

devtools::check()
devtools::test()


# devtools::check_win_devel()
# rhub::check_for_cran()


devtools::build()
devtools::build_manual()
devtools::release_checks()
#devtools::spell_check()


devtools::submit_cran()

