# Define immature data:
s <- read.scsset(years, survey = "regular", valid = 1)
s$grid <- deg2grid(lon(s), lat(s))
b <- read.scsbio(years, survey = "regular", sex = sex)
b$tow.id <- tow.id(b)
if (sex == 1) b$maturity <- morphometric.maturity(b) else b$maturity <- is.mature(b)
b <- b[which(b$carapace.width >= 2), ]
import(s, fill = 0) <- freq(b[which(!b$maturity), ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
f_imm <- NULL
for (i in 1:length(years)){
  ss <- s[year(s) == years[i], c("grid", fvars)]
  if (years[i] <= 2010) ss <- aggregate(ss[fvars], by = ss["grid"], mean)
  tmp <- apply(ss[fvars], 2, mean)
  tmp[setdiff(as.character(seq(5, xlim[2]-20, by = step)), names(tmp))] <- 0
  tmp <- tmp[order(as.numeric(names(tmp)))]
  tmp <- t(tmp)
  rownames(tmp) <- years[i]
  f_imm <- rbind(f_imm, tmp)
}

# Define skip data:
s <- read.scsset(years, survey = "regular", valid = 1)
s$grid <- deg2grid(lon(s), lat(s))
b <- read.scsbio(years, survey = "regular", sex = sex)
b$tow.id <- tow.id(b)
if (sex == 1) b$maturity <- morphometric.maturity(b) else b$maturity <- is.mature(b)
b <- b[which(b$carapace.width >= 2), ]
import(s, fill = 0) <- freq(b[which(!b$maturity & !is.new.shell(b)), ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
f_skp <- NULL
for (i in 1:length(years)){
  ss <- s[year(s) == years[i], c("grid", fvars)]
  if (years[i] <= 2010) ss <- aggregate(ss[fvars], by = ss["grid"], mean)
  tmp <- apply(ss[fvars], 2, mean)
  tmp[setdiff(as.character(seq(5, xlim[2]-20, by = step)), names(tmp))] <- 0
  tmp <- tmp[order(as.numeric(names(tmp)))]
  tmp <- t(tmp)
  rownames(tmp) <- years[i]
  f_skp <- rbind(f_skp, tmp)
}

# Define mature recruitment data:
s <- read.scsset(years, survey = "regular", valid = 1)
s$grid <- deg2grid(lon(s), lat(s))
b <- read.scsbio(years, survey = "regular", sex = sex)
b$tow.id <- tow.id(b)
if (sex == 1) b$maturity <- morphometric.maturity(b) else b$maturity <- is.mature(b)
b <- b[which(b$carapace.width >= 35), ]
ix <- which(b$maturity & is.new.shell(b))
import(s, fill = 0) <- freq(b[ix, ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
f_rec <- NULL
for (i in 1:length(years)){
  ss <- s[year(s) == years[i], c("grid", fvars)]
  if (years[i] <= 2010) ss <- aggregate(ss[fvars], by = ss["grid"], mean)
  tmp <- apply(ss[fvars], 2, mean)
  tmp[setdiff(as.character(seq(3, xlim[2]-20, by = step)), names(tmp))] <- 0
  tmp <- tmp[order(as.numeric(names(tmp)))]
  tmp <- t(tmp)
  rownames(tmp) <- years[i]
  f_rec <- rbind(f_rec, tmp)
}

# Define mature residual data:
s <- read.scsset(years, survey = "regular", valid = 1)
s$grid <- deg2grid(lon(s), lat(s))
b <- read.scsbio(years, survey = "regular", sex = sex)
b$tow.id <- tow.id(b)
if (sex == 1) b$maturity <- morphometric.maturity(b) else b$maturity <- is.mature(b)
b <- b[which(b$carapace.width >= 35), ]
ix <- which(b$maturity & !is.new.shell(b))
import(s, fill = 0) <- freq(b[ix, ], by = key(s), step = step)
fvars <- names(s)[gsub("[0-9.]", "", names(s)) == ""]
s[fvars] <- 1000000 * s[fvars] / repvec(s$swept.area, ncol = length(fvars))
f_res <- NULL
for (i in 1:length(years)){
  ss <- s[year(s) == years[i], c("grid", fvars)]
  if (years[i] <= 2010) ss <- aggregate(ss[fvars], by = ss["grid"], mean) 
  tmp <- apply(ss[fvars], 2, mean)
  tmp[setdiff(as.character(seq(3, xlim[2]-20, by = step)), names(tmp))] <- 0
  tmp <- tmp[order(as.numeric(names(tmp)))]
  tmp <- t(tmp)
  rownames(tmp) <- years[i]
  f_res <- rbind(f_res, tmp)
}

# Total matures:
f_mat <- f_rec + f_res

# Define data:
data <- list(x_imm    = as.numeric(repvec(as.numeric(colnames(f_imm)), nrow = nrow(f_imm))),
             f_imm    = as.numeric(f_imm),
             year_imm = as.numeric(repvec(as.numeric(rownames(f_imm)), ncol = ncol(f_imm))) - min(years),
            
              x_skp    = as.numeric(repvec(as.numeric(colnames(f_skp)), nrow = nrow(f_skp))),
             f_skp    = as.numeric(f_rec),
             year_skp = as.numeric(repvec(as.numeric(rownames(f_skp)), ncol = ncol(f_skp))) - min(years),
             
             x_rec    = as.numeric(repvec(as.numeric(colnames(f_rec)), nrow = nrow(f_rec))),
             f_rec    = as.numeric(f_rec),
             year_rec = as.numeric(repvec(as.numeric(rownames(f_rec)), ncol = ncol(f_rec))) - min(years),
             
             x_res    = as.numeric(repvec(as.numeric(colnames(f_res)), nrow = nrow(f_res))),
             f_res    = as.numeric(f_res),
             year_res = as.numeric(repvec(as.numeric(rownames(f_res)), ncol = ncol(f_res))) - min(years),
             
             x_mat    = as.numeric(repvec(as.numeric(colnames(f_mat)), nrow = nrow(f_mat))),
             f_mat    = as.numeric(f_mat),
             year_mat = as.numeric(repvec(as.numeric(rownames(f_mat)), ncol = ncol(f_mat))) - min(years),
             
             delta_x  = step)

save(data, file = "females.2006-2020.rdata")

