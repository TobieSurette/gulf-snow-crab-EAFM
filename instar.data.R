library(gulf.data)
library(gulf.stats)

# Define data:
s <- read.scsset(years, survey = "regular", valid = 1)
b <- read.scsbio(years, survey = "regular", sex = sex)
b$tow.id <- tow.id(b)
if (sex == 1) b$maturity <- morphometric.maturity(b) else b$maturity <- is.mature(b)
b <- b[which(b$carapace.width >= 2), ]
import(s, fill = 0) <- freq(b[which(!b$maturity), ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
fi <- apply(s[fvars], 2, mean)
fi[setdiff(as.character(seq(5, xlim[2]-20, by = step)), names(fi))] <- 0
fi <- fi[order(as.numeric(names(fi)))]

s <- read.scsset(years+1, survey = "regular", valid = 1)
b <- read.scsbio(years+1, survey = "regular", sex = sex)
b$tow.id <- tow.id(b)
if (sex == 1) b$maturity <- morphometric.maturity(b) else b$maturity <- is.mature(b)
b <- b[which(b$carapace.width >= 35), ]
ix <- which(b$maturity & is.new.shell(b))
import(s, fill = 0) <- freq(b[ix, ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
fm <- apply(s[fvars], 2, mean)
#fm <- length(ix) * (fm / sum(fm))
fm[setdiff(as.character(seq(10, xlim[2], by = step)), names(fm))] <- 0
fm <- fm[order(as.numeric(names(fm)))]

# Define data:
data <- list(x_imm = as.numeric(names(fi)),
             f_imm = as.numeric(fi),
             x_mat = as.numeric(names(fm)),
             f_mat = as.numeric(fm),
             dx = step)
