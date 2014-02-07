# A sloppy script to compare a few ways to model the histograms.
# Not necessary for running the pipeline, but good to have around in case
# I want to reproduce the comparison or add new models to compare.

# Given a list of frequencies, allele lengths (NOT normalized to the major
# length and the length of the major allele, plot the raw data as a histogram
# and several model fits for the data as colored lines.
mu <- function(f, a) {
    sum(f*a) / sum(f)
}

v <- function(f, a) {
    sum(f*(a - mu(f,a))**2) / (sum(f) - 1)
}

f <- function(f, a, major) {
    require('plotrix')
    xs <- 0:max(a)
    barp(f/sum(f), x=a, xlim=c(min(xs), max(xs)), ylog=TRUE, main=sprintf("Model fits for major allele len=%d", major), xlab="Observed allele lengths (units)", ylab="Sex chromosome data frequencies and model probabilities")
    m <- mu(f, a)
    v <- v(f, a)
    print(sprintf("data mean=%f", m))
    print(sprintf("data variance=%f", v))
    lines(dpois(xs, lambda=m), x=xs, col="green", type="b")
    lines(dnorm(xs, mean=m, sd=sqrt(v)), x=xs, col="red", type="b")
    
    # build the double negbin distn
    #   in gain, a=0 means a +1 gain, a=1 means +2 gain, etc.
    gain <- data.frame(f=f[a>major], a=a[a>major] - major - 1)
    #   in loss, a=0 means a -1 loss, a=2 means -2 loss, etc.
    loss <- data.frame(f=f[a<major], a=major - a[a<major] - 1)
    id <- data.frame(f=f[a==major], a=a[a==major] - major)
    mgain <- mu(gain$f, gain$a)
    vgain <- v(gain$f, gain$a)
    disp.gain <- mgain^2 / (vgain - mgain)
    mloss <- mu(loss$f, loss$a)
    vloss <- v(loss$f, loss$a)
    disp.loss <- mloss^2 / (vloss - mloss)
    p.id <- sum(id$f) / (sum(gain$f) + sum(loss$f) + sum(id$f))
    p.gain <- sum(gain$f) / (sum(gain$f) + sum(loss$f) + sum(id$f))
    p.loss <- sum(loss$f) / (sum(gain$f) + sum(loss$f) + sum(id$f))
    print(sprintf("double negbin: p.id=%f, p.gain=%f, p.loss=%f, mgain=%f, vgain=%f, disp.gain=%f, mloss=%f, vloss=%f, disp.loss=%f", p.id, p.gain, p.loss, mgain, vgain, disp.gain, mloss, vloss, disp.loss))
    double.negbin.pmf <- function(x) {
        if(x > major) {
            p.gain*dnbinom(x - major - 1, mu=mgain, size=disp.gain)
        } else if(x < major) {
            p.loss*dnbinom(major - x - 1, mu=mloss, size=disp.loss)
        } else {
            p.id
        }
    }
    lines(sapply(xs, double.negbin.pmf), x=xs, col="blue", type="b")
    
    # bin+negbin model is a binomial loss model and the same negbin gain model as double negbin
    bin.p.loss <- sum(loss$f*loss$a) / (sum(loss$f)*major)
    print(sprintf("bin.nebgin: bin.p.loss=%f", bin.p.loss))
    bin.negbin.pmf <- function(x) {
        if(x > major) {
            p.gain*dnbinom(x - major - 1, mu=mgain, size=disp.gain)
        } else if (x < major) {
            p.loss*dbinom(major - x - 1, size=major, prob=bin.p.loss)
        } else {
            p.id
        }
    }
    lines(sapply(xs, bin.negbin.pmf), x=xs, col="orange", type="b")

    legend("topright", legend=c("poisson", "normal", "double negbin", "bin.negbin"), fill=c("green", "red", "blue", "orange"), title="model")
}

sqr <- function(x) x^2
cub <- function(x) x^3
quar <- function(x) x^4
quin <- function(x) x^5

model.mean <- function(est, data, title='') {
    df <- est
    df$mean <- abs(df$mean) - 1
    N <- sapply(data, sum)
    plot(x=df$major, y=df$mean, main=title)
    segments(df$major, df$mean+sqrt(df$var/N), df$major, df$mean-sqrt(df$var/N))
    errbar.width = 0.1
    segments(df$major-errbar.width, df$mean+sqrt(df$var/N), df$major+errbar.width, df$mean+sqrt(df$var/N))
    segments(df$major-errbar.width, df$mean-sqrt(df$var/N), df$major+errbar.width, df$mean-sqrt(df$var/N))

    plot.model <- function(g, col) {
        m = lm(mean ~ g(major), data=df, weights=sapply(data, sum))
        lines(x=df$major, y=sapply(df$major, function(x) sum(coef(m)*c(1,g(x)))), col=col)
    }

    plot.model(identity, "red")
    plot.model(sqrt, "green")
    plot.model(log, "blue")
    plot.model(sqr, "orange")
    #plot.model(cub, "purple")
    #plot.model(quar, "grey")
    #plot.model(exp, "cyan")
}
    
model.var <- function(est, data, title) {
    df <- est
    plot(x=df$major, y=df$var, main=title)

    plot.model <- function(g, col) {
        m = lm(var ~ g(major), data=df, weights=sapply(data, sum))
        lines(x=df$major, y=sapply(df$major, function(x) sum(coef(m)*c(1,g(x)))), col=col)
    }

    plot.model(identity, "red")
    plot.model(sqrt, "green")
    plot.model(log, "blue")
    plot.model(sqr, "orange")
    #plot.model(cub, "purple")
    #plot.model(quar, "grey")
    #plot.model(exp, "cyan")
}
    
# Optional, generate the plots.  the above function is generally useful
# without this data
#load("hists3.RData")
load("hists_cigar_bin0.RData")
#f(mono.15.freq, mono.15.allele, 15)
layout(matrix(1:4, ncol=2))
binned.data <- list(mono.binned.data, di.binned.data, tri.binned.data, tetra.binned.data)
gains <- list(mono.gain.ests, di.gain.ests, tri.gain.ests, tetra.gain.ests)
losses <- list(mono.loss.ests, di.loss.ests, tri.loss.ests, tetra.loss.ests)
titles <- c("mono", "di", "tri", "tetra")
for (i in 1:4) {
    par(mar=c(3,4,2,1))
    model.mean(gains[[i]], binned.data[[i]], titles[i])
}
#model.mean(di.gain.ests, di.binned.data)
#model.var(di.loss.ests, di.binned.data)
