#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

braak <- suppressWarnings(as.numeric(args[1]))
cerad <- suppressWarnings(as.numeric(args[2]))
cdr <- suppressWarnings(as.numeric(args[3]))

status <- "undetermined"

if (cdr %in% c(0, 0.5)){
    if (cerad %in% c(1, 4)){
      if (braak %in% c(0, 1, 2)){
        status <- "Not"
      }else if(braak %in% c(3, 4)){
        status <- "Not"
      }else if(braak %in% c(5, 6)){
        status <- "Not"
      }
    } 
  } else if (cdr %in% c(1, 2)){
    if (cerad %in% c(1, 4)){
      if (braak %in% c(0, 1, 2)){
        status <- "Low"
      }else if(braak %in% c(3, 4)){
        status <- "Low"
      }else if(braak %in% c(5, 6)){
        status <- "Low"
      }
    }else if (cerad %in% c(2, 3)){
      if (braak %in% c(0, 1, 2)){
        status <- "Low"
      }else if(braak %in% c(3, 4)){
        status <- "Intermediate"
      }else if(braak %in% c(5, 6)){
        status <- "Intermediate"
      }
    }
  }else if (cdr %in% c(3)){
    if (cerad %in% c(1, 2, 3, 4)){
      if (braak %in% c(0, 1, 2)){
        status <- "Low"
      }else if(braak %in% c(3, 4)){
        status <- "Intermediate"
      }else if(braak %in% c(5, 6)){
        status <- "Intermediate"
      }
    }
  }else if (cdr %in% c(4, 5)){
    if (cerad %in% c(1, 4)){
      if (braak %in% c(0, 1, 2)){
        status <- "Low"
      }else if(braak %in% c(3, 4)){
        status <- "Intermediate"
      }else if(braak %in% c(5, 6)){
        status <- "Intermediate"
      }
    }else if (cerad %in% c(2, 3)){
      if (braak %in% c(0, 1, 2)){
        status <- "Low"
      }else if(braak %in% c(3, 4)){
        status <- "Intermediate"
      }else if(braak %in% c(5, 6)){
        status <- "High"
      }
    }
  }
}

cat(status, "\n")
