# Compte-rendu
# Module ADC : Classification supervisée
# Binome : Simon Besson-Girard & Nicolas Wattiez
# Semestre d'automne 2014 - Master ecosciences-microbiologie (MIV)

#install.packages("kernlab")
#install.packages("ROCR")
#install.packages("plotrix")
library(MASS)
library(kernlab)
library(ROCR)
library(plotrix)

########################################################################
#	SCRIPT PREMIERE PARTIE : CORRECTION DE L'ENSEIGNANT
########################################################################

### Q1.6
generateMultivariateGaussian <- function(n) {
    A <- matrix(c(2,1,1,0),nrow=2,ncol=2)
    Z <- matrix(c(rnorm(2 * n)),ncol=2)
    mu <- c(1,1)
    Z %*% A + mu
}

#generateMultivariateGaussian(10)#

#plot(generateMultivariateGaussian(1000))#

### Q 2.1
generateSimpleDataset <- function(n) {
    np <- ceiling(n/2)
    nn <- floor(n/2)
    pos <- mvrnorm(np,mu=c(0,0),Sigma=matrix(c(1,0.5,0.5,1),2,2))
    neg <- mvrnorm(nn,mu=c(6,2),Sigma=matrix(c(0.2,0.4,0.4,4),2,2))
    x <- as.data.frame(rbind(pos,neg))
    cbind(x, Y = c(rep(1,np),rep(-1,nn)))
}

#generateSimpleDataset(11)#
    
### Q 2.2
plot.dataset <- function(d) {
    lim <- c(min(d$V1,d$V2), max(d$V1,d$V2))
    plot(c(), c(),
         xlim=lim, # même échelle en abscisses et ordonnées
         ylim=lim,
         xlab=expression(x[1]),
         ylab=expression(x[2]))
    
    points(d[d$Y==1,c("V1","V2")],col='red',lw=2,pch="*")
    points(d[d$Y==-1,c("V1","V2")],col='blue',lw=2,pch="o")
}

#plot.dataset(generateSimpleDataset(1000))#

### Q 2.3

### Une alternative utile en pratique, non demandée dans le devoir
#generateDifficultDatasetAlt <- function(n) {
#    w <- 2 * runif(2) - 1
#    points <- 4 * cbind(runif(n),runif(n)) - 2
#    p <- 1 / (1 + exp(3 * (points %*% w)))
#    Y <- Vectorize(function(p) rbinom(1,1,p))(p)
#    cbind(as.data.frame(points), Y = 2*Y - 1)
#}
#plot.dataset(generateDifficultDatasetAlt(200))#

generateDifficultDatasetAlt <- function(d,n) {
    set.seed(42)
    w <- 2 * runif(d) - 1
    w <- w / sqrt(sum(w^2))
    points <- sapply(1:n,function(i) 4 * runif(d) - 2)
    p <- 1 / (1 + exp(0.5 * (t(points) %*% w)))
    Y <- Vectorize(function(p) rbinom(1,1,p))(p)
    cbind(as.data.frame(t(points)), data.frame(Y = 2*Y - 1))
}

generateDifficultDataset <- function(n) {
    np <- ceiling(n/2)
    nn <- floor(n/2)
    pos <- mvrnorm(np,mu=c(0,0),Sigma=matrix(c(4,2,2,4),2,2))
    neg <- mvrnorm(nn,mu=c(4,2),Sigma=matrix(c(0.2,-0.4,-0.4,4),2,2))
    x <- as.data.frame(rbind(pos,neg))
    cbind(x, Y = c(rep(1,np),rep(-1,nn)))
}

#plot.dataset(generateDifficultDataset(200))#

### Q 3.4
mediator.demo <- function(d) {
    plot.dataset(d)

    c.pos <- t(colMeans(d[d$Y==1,c("V1","V2")]))
    c.neg <- t(colMeans(d[d$Y==-1,c("V1","V2")]))
    
    points(c.pos, pch=7, lwd=2, col="red")
    points(c.neg, pch=7,lwd=2,col="blue")

    w <- c.pos - c.neg
    b <- (sum(c.neg ^ 2) - sum(c.pos ^ 2)) / 2
    abline(b = - w[1] / w[2], a = - b / w[2])
}

#mediator.demo(generateDifficultDataset(500))#

### Q 3.5
mediatorHyperplane <- function(d) {
    c.pos <- t(colMeans(d[d$Y==1, colnames(d) != "Y"]))
    c.neg <- t(colMeans(d[d$Y==-1,colnames(d) != "Y"]))

    w <- c.pos - c.neg
    b <- (sum(c.neg ^ 2) - sum(c.pos ^ 2)) / 2
    f <- function(x) {
        as.matrix(x[,colnames(d) != "Y"]) %*% t(w) + b # alternativement : utiliser apply
        # On prévient la possibilité que x contienne
        # une colonne Y, ce qui est le cas dans les tests que nous
        # faisons plus bas
    }
    pred <- function(x) {
        y <- f(x) >= 0
        y - !y
    }
    list(w = w, b = b, f = f, pred = pred)
}

#train_set <- generateDifficultDataset(20)#
#test_set <- generateDifficultDataset(100)#
#preds <- mediatorHyperplane(train_set)$pred(test_set)#

### Q 4.1
#confusion_matrix <- table(preds, test_set$Y)

### Q 4.2
# prend une matrice de confusion en argument
accuracy <- function(M) {
    (M[2,2] + M[1,1]) / sum(M)
}

sensitivity <- function(M) {
    M[2,2] / sum(M[,2])
}

specificity <- function(M) {
    M[1,1] / sum(M[,1])
}

precision <- function(M) {
    M[2,2] / sum(M[2,])
}

erreur <- function(M) {
	(M[1,2] + M[2,1]) / sum(rowSums(M))
}

#accuracy(confusion_matrix)#
#sensitivity(confusion_matrix)#
#specificity(confusion_matrix)#
#precision(confusion_matrix)#

evaluation <- function(preds,labels) {
    confusion_matrix <- table(preds, labels)
    list(
        accuracy = accuracy(confusion_matrix),
        sensitivity = sensitivity(confusion_matrix),
        specificity = specificity(confusion_matrix),
        precision = precision(confusion_matrix),
        error = erreur(confusion_matrix)
    )
}

### Q 4.3

# l'option pred permet de préciser si l'on souhaite des prédictions ou
# juste des scores
holdOut <- function(f, d_train, d_test, pred = TRUE) {
    if(pred) f(d_train)$pred(d_test)
    else f(d_train)$f(d_test)
}

#d_train <- generateDifficultDataset(200)#
#d_test <- generateDifficultDataset(100)#
#preds <- holdOut(mediatorHyperplane,d_train,d_test)#
#evaluation(preds,d_test$Y)#

### Q4.4
#rocr_pred <- prediction(preds,d_test$Y)#
# la fonction performance calcule des courbes de performance, de type
# ROC. Nous lui donnons ici des prédictions au lieu de scores, les
# courbes ont donc trois points, obtenus pour les seuils -1, 0 et
# 1. C'est à chaque fois le deuxième point qui nous intéresse.
#performance(rocr_pred,"acc")@y.values[[1]][2]#
#performance(rocr_pred,"sens")@y.values[[1]][2]#
#performance(rocr_pred,"spec")@y.values[[1]][2]#
#performance(rocr_pred,"prec")@y.values[[1]][2]#

### Q 4.5

# ATTENTION: pour la courbe ROC, il faut utiliser les scores de
# prédiction !
#pred_scores <- holdOut(mediatorHyperplane,d_train,d_test,pred=F)#
#plot(performance(prediction(pred_scores,d_test$Y),"sens","spec"))#


CrossValidation <- function(f,d,N,pred = TRUE) {
    n <- dim(d)[1]
    permutation <- sample(1:n)
    d <- d[permutation,]

    fold <- function(i) {
        a <- round(n * (i - 1) / N + 1)
        b <- round(n * i / N)
        test.idx <- a:b
        train_set <- d[- test.idx,]
        test_set <- d[test.idx,]
        classifier <- f(train_set)
        if(pred) classifier$pred(test_set)
        else classifier$f(test_set)
    }

    preds <- unlist(lapply(1:N,fold))
    
    # pour finir, on remet les prédictions dans le bon ordre
    preds[order(permutation)]
}

### Q 4.7
#dataset <- generateDifficultDataset(100)#
#preds <- CrossValidation(mediatorHyperplane,dataset,7)#
#evaluation(preds,dataset$Y)#

error.wrt.n <- function(n,gen) {
    f <- function() {
        d <- gen(n)
        preds <- CrossValidation(mediatorHyperplane,d,5)
        evaluation(preds,d$Y)$error
    }
    k <- 30
    error <- replicate(k,f())
    # on renvoie la moyenne sur 30 essais et l'IC à 95%
    error_bar <- mean(error)
    eps <- qnorm(0.975) * sd(error) / sqrt(k)
    c(error_bar, error_bar - eps, error_bar + eps)
}

#error.wrt.n(20,generateDifficultDataset)#

plot.error.wrt.n <- function(gen) {
    points <- seq(from=5, to=50, by=5)
    data <- sapply(points, function(n) error.wrt.n(n,gen))
    plotCI(points,data[1,],li=data[2,],ui=data[3,])
}

#plot.error.wrt.n(generateDifficultDataset)#

# La même chose mais avec la méthode du hold-out
error.wrt.n <- function(n,gen) {
    f <- function() {
        d_test <- gen(1000)
        preds <- holdOut(mediatorHyperplane,gen(n),d_test)
        evaluation(preds,d_test$Y)$accuracy
    }
    k <- 100
    error <- replicate(k,f())
    # on renvoie la moyenne sur 30 essais et l'IC à 95%
    error_bar <- mean(error)
    eps <- qnorm(0.975) * sd(error) / sqrt(k)
    c(error_bar, error_bar - eps, error_bar + eps)
}

#error.wrt.n(20,generateDifficultDataset)#

plot.error.wrt.n <- function(gen) {
    points <- seq(from=5, to=50, by=5)
    data <- sapply(points, function(n) error.wrt.n(n,gen))
    plotCI(points,data[1,],li=data[2,],ui=data[3,])
}

#plot.error.wrt.n(generateDifficultDataset)#


########################################################################
#	5
########################################################################

###
SVM.accuracy.wrt.C <- function(d) {
	C.values <- sapply(c(seq(-10,10,le=500)), function (x) 2^x)
	svm.cross <- sapply(C.values,function(x) cross(ksvm(Y~.,data=d,type='C-svc',kernel='vanilladot',C=x,cross=20)) )
	svm.error <- sapply(C.values,function(x) error(ksvm(Y~.,data=d,type='C-svc',kernel='vanilladot',C=x,cross=20)) ) 
	plot(C.values,svm.cross,type='o',xlab="Valeurs de C",ylab="Erreur",cex=.5,ylim=c(min(c(svm.error,svm.cross)),max(c(svm.cross,svm.error))))# ,log="x")
	points(C.values,svm.error,type='o',col='blue',cex=.5)
	legend("topright",c("Erreur de validation croisée","Erreur à l'apprentissage"),fill=c("black","blue"))
}

#SVM.accuracy.wrt.C(generateDifficultDatasetAlt(100,30))#

###
selectC <- function(d) {
	C.values <- sapply(c(seq(-10,10,by=1)), function (x) 2^x)
	svm.cross <- sapply(C.values,function(x) cross(ksvm(Y~.,data=d,type='C-svc',kernel='vanilladot',C=x,cross=20)) )
	tab <- as.data.frame(cbind(C.values,svm.cross))
	return(tab$C.values[tab$svm.cross==min(svm.cross)])
}

#(resultat.selectC <- selectC(generateDifficultDatasetAlt(100,20)))#

###
compare.SVM.mH <- function(nbjeux,fonctionquigenere,taillejeu) {
	f <- function() {
		d <- fonctionquigenere(taillejeu)
		C.value <- selectC(d)
		error.SVM <- mean(sapply(C.values,function(x) cross(ksvm(Y~.,data=d,type='C-svc',kernel='vanilladot',C=C.value,cross=5))) )
		error.mediatorHyperplane <- error.wrt.n(taillejeu,fonctionquigenere)[1]
		return(c(error.SVM,error.mediatorHyperplane))
	}
	moyennes <- replicate(nbjeux,f())
	return(c(mean(moyennes[1,]),(mean(moyennes[2,]))))
}

###
plot.ROC.SVM <- function() {
	predictor <- prediction(f.Difficult(dataset.generateDifficultDataset[,1:2]),generateDifficultDataset(100)$y)
	perfector <- performance(predictor, measure = "sens", x.measure = "spec") 
	plot(perfector)
}

###
read.prostate.dataset <- function(nomfichier) {
	prostate <- read.table(nomfichier)
	colnames(prostate)[1] <- "Y"
	for (i in 1:nrow(prostate)) {
		if (prostate[i,1]==0) {prostate[i,1] <- -1}
	}
	return(prostate)
}

###
read.prostate.dataset <- function(nomfichier) {
	prostate <- read.table(nomfichier)
	colnames(prostate)[1] <- "Y"
	for (i in 1:nrow(prostate)) {if (prostate[i,1]==0) {prostate[i,1] <- -1} }
	prostate.moyenne <- sapply((1:ncol(prostate)),function(x) mean(prostate[,x]) )
	prostate.ecart.type <- sapply((1:ncol(prostate)),function(x) sd(prostate[,x]) )
	for (c in 2:ncol(prostate)) {for (l in 1:nrow(prostate)) {prostate[l,c] <- (prostate[l,c]-prostate.moyenne[c])/prostate.ecart.type[c]} }
	return(prostate)
}

###
illustration <- function() {
as.vector(t(pred.prostate.mediatorHyperplane <- mediatorHyperplane(prostate)$pred(prostate)))
 #selectC vaut 2
 
plot(0,0,tck=1,las=1,xlim=c(0,6),ylim=c(0,6),col.axis="white",col.lab="white",pch=46,bty="n",bty="o",fg="grey")
arrows(0,0,2,0,col="red")
arrows(0,0,0,2,col="red")
arrows(3,0,3,1,col="blue")
arrows(3,0,4,0,col="blue")
text(0,2.2,"V1")
text(2.2,0,"V2")
text(3,1.2,"V3")
text(4.2,0,"V4")
arrows(0,0,1,2,col="black")
arrows(3,0,4,2,col="black")
text(1.1,2.1,"u")
text(4.1,2.1,"v")
}
