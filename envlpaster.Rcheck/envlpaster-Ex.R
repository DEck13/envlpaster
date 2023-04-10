pkgname <- "envlpaster"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('envlpaster')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("fitness_boot")
### * fitness_boot

flush(stderr()); flush(stdout())

### Name: fitness_boot
### Title: fitness_boot
### Aliases: fitness_boot

### ** Examples

## Not run: # see Web-based Supplementary Materials for ``Enveloping the aster model.''



cleanEx()
nameEx("manifold1Dplus")
### * manifold1Dplus

flush(stderr()); flush(stdout())

### Name: manifold1Dplus
### Title: manifold1Dplus
### Aliases: manifold1Dplus

### ** Examples

## Not run: 
##D library(envlpaster)
##D  data(simdata30nodes)
##D  data <- simdata30nodes.asterdata
##D  nnode <- length(vars)
##D  xnew <- as.matrix(simdata30nodes[,c(1:nnode)])
##D  m1 <- aster(xnew, root, pred, fam, modmat)
##D  avar <- m1$fisher
##D  beta <- m1$coef
##D  U <- beta \##D 
##D  manifold1Dplus(M = avar, U = U, u = 1)
## End(Not run)



cleanEx()
nameEx("scanner")
### * scanner

flush(stderr()); flush(stdout())

### Name: scanner
### Title: scanner
### Aliases: scanner

### ** Examples

## Not run: 
##D library(envlpaster)
##D data(simdata30nodes)
##D data <- simdata30nodes.asterdata
##D nnode <- length(vars)
##D xnew <- as.matrix(simdata30nodes[,c(1:nnode)])
##D m1 <- aster(xnew, root, pred, fam, modmat)
##D avar <- m1$fisher
##D beta <- m1$coef
##D scanner(M = avar, coef = beta, u = 1)
## End(Not run)



cleanEx()
nameEx("secondboot")
### * secondboot

flush(stderr()); flush(stdout())

### Name: secondboot
### Title: secondboot
### Aliases: secondboot

### ** Examples

### Web-based Supplementary Materials for ``Enveloping the aster model.'' ###



cleanEx()
nameEx("selection")
### * selection

flush(stderr()); flush(stdout())

### Name: selection
### Title: selection
### Aliases: selection

### ** Examples

## Not run: 
##D set.seed(13)
##D library(envlpaster)
##D library(aster2)
##D data(generateddata)
##D m.null <- aster(resp ~ 0 + varb, fam = fam, pred = pred,
##D                 varvar = varb, idvar = id, root = root, data = redata)
##D m1 <- aster(resp ~ 0 + varb + mass + timing,
##D             fam = fam, pred = pred, varvar = varb, idvar = id,
##D             root = root, data = redata)
##D m2 <- aster(resp ~ 0 + varb + mass + timing +
##D               I(mass^2) + I(timing^2) + I(mass*timing),
##D             fam = fam, pred = pred, varvar = varb, idvar = id,
##D             root = root, data = redata)
##D anova.table <- anova(m.null,m1,m2); anova.table
##D beta <- m1$coef
##D a <- grepl( "offsp", names(beta))
##D a <- a + grepl( "surviv", names(beta))
##D b <- which(a == 1)
##D target <- c(1:length(beta))[-b]
##D nnode <- ncol(m1$x)
##D data.aster <- asterdata(data, vars, pred, rep(0,nnode),
##D                         fam, families = list("bernoulli", "poisson",
##D                                              fam.zero.truncated.poisson()))
##D selection(parm  = beta, index = target, model = m1,
##D           data = data.aster, alpha = 0.05, type = "canonical",
##D           method = "eigen")
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
