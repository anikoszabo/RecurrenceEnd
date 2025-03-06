library(devtools)

use_package("nspmix")
use_package("reda")
use_package("survival")
use_build_ignore("RE_develop.R")
use_gpl3_license()
use_github()

re <- as.package(".")
load_all(re)
document(re)

check(re)
