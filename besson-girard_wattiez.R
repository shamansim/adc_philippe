# Compte-rendu
# Module ADC : Classification supervisée
# Binome : Simon Besson-Girard & Nicolas Wattiez
# Semestre d'automne 2014 - Master ecosciences-microbiologie (MIV)

library(MASS)
install.packages("kernlab")
library(kernlab)
install.packages("ROCR")
library(ROCR)

########################################################################
#	1
########################################################################

generation_multivariee<-function(n,mu,sigma){
	#on cree deux vecteurs qui suivent chacun une loi normale centree reduite
	z1<-rnorm(n,0,1)
	z2<-rnorm(n,0,1)
	#on rassemble les deux variables pour creer un vecteur aleatoire  
	Z<-rbind(z1,z2)
	#on applique la factorisation de Cholesky
	A<-chol(sigma)
	#on applique la formule du 1.3
	d<-t(A)%*%Z+mu
	d<-(d)
	return(d)
}

########################################################################
#	2
########################################################################

#génere un jeu de données
generateSimpleDataset<-function(n) {
	#on definit deux fonctions qui fabriquent un echantillon suivant deux loi gaussiennes lineairement separables
	positifs.echantillon<-function(n) mvrnorm(n,mu=c(0,0),Sigma=matrix(c(1,.5,.5,1),2,2))
	negatifs.echantillon<-function(n) mvrnorm(n,mu=c(3,2),Sigma=matrix(c(.2,.4,.4,4),2,2))
	#on definit deux echantillons ayant chacun la moitie de l'effectif total desire
	positifs<-positifs.echantillon(n/2)
	negatifs<-negatifs.echantillon(n/2)
	colnames(positifs)<-c("x1","x2")
	#soit x le tableau de donnees qui pour chaque ligne indique la provenance de l'echantillon
	x<-as.data.frame(rbind(positifs,negatifs))
	cbind(x,y=c(rep(1,n/2),rep(-1,n/2)))
}

#on génère un jeu de données qui ne variera pas
dataset.generateSimpleDataset=generateSimpleDataset(100)#

#affiche le graphique coloré par groupe d'un jeu de données
plot.dataset<-function(d) {
	#parametres du plot ajustes aux donnees
	lim<- c(min(d$x1,d$x2),max(d$x1,d$x2))
	plot(c(),c(),
		xlim=lim,
		ylim=lim,
		xlab=expression(x[1]),
		ylab=expression(x[2]))
	#ajout des points avec un motif et une couleur selon leur echantillon de provenance
	points(d[d$y==1,c("x1","x2")],col="red",lw=2,pch="+")
	points(d[d$y==-1,c("x1","x2")],col="blue",lw=2,pch="o")
}

#variante de generateSimpleDataset
generateDifficultDataset<-function(n) {
	#on definit deux fonctions qui fabriquent un echantillon suivant deux loi gaussiennes specifiees dans l'enonce
	pos.sample<-function(n) mvrnorm(n,mu=c(0,0),Sigma=matrix(c(4,2,2,4),2,2))
	neg.sample<-function(n) mvrnorm(n,mu=c(4,2),Sigma=matrix(c(.2,-.4,-.4,4),2,2))
	#on definit deux echantillons ayant chacun la moitie de l'effectif total desire
	pos<-pos.sample(n/2)
	neg<-neg.sample(n/2)
	colnames(pos)<-c("x1","x2")
	#soit x le tableau de donnees qui pour chaque ligne indique la provenance de l'echantillon
	x<-as.data.frame(rbind(pos,neg))
	cbind(x,y=c(rep(1,n/2),rep(-1,n/2)))
}

#on génère un jeu de données qui ne variera pas
dataset.generateDifficultDataset=generateDifficultDataset(100)#

########################################################################
#	3
########################################################################
	
#ajoute les centres d'inertie
plot.add.inertie<-function(d) {
	#pour le groupe des positifs, calcul de la moyenne des abscisses puis des ordonnees
	centre.pos.abs=mean(d$x1[d$y==1])
	centre.pos.ord=mean(d$x2[d$y==1])
	#pour le groupe des negatifs, calcul de la moyenne des abscisses puis des ordonnees
	centre.neg.abs=mean(d$x1[d$y==-1])
	centre.neg.ord=mean(d$x2[d$y==-1])
	#dessine les centres
	points(centre.pos.abs,centre.pos.ord,col="green",pch=16)
	points(centre.neg.abs,centre.neg.ord,col="green",pch=16)
}

#calcule le w et le b et trace la mediatrice
plot.add.mediatrice<-function(d) {
	#
	c.pos<-t(colMeans(d[d$y==1,c("x1","x2")]))
	c.neg<-t(colMeans(d[d$y==-1,c("x1","x2")]))
	#calcul de w et de b
	w<-c.pos-c.neg
	b<-(sum(c.neg^2) - sum(c.pos^2))/2
	#on trace la mediatrice
	abline(b=-w[1]/w[2],a=-b/w[2])
}

#fonction mediator.demo(d)
mediator.demo<-function(d){
	plot.dataset(d)
	plot.add.inertie(d)
	plot.add.mediatrice(d)
} 


mediatorHyperplane<-function(d){
			#pour s'adapter aux dimensions on recupere le nombre de dimensions+1
			l=length(colnames(d))
			#on suppose que l'etiquette s'appelle forcement y et que c'est la derniere colonne
			c.pos<-t(colMeans(d[d$y==1,colnames(d[-l])]))
			c.neg<-t(colMeans(d[d$y==-1,colnames(d[-l])]))
			w<-c.pos-c.neg
			b<-(sum(c.neg^2) - sum(c.pos^2))/2
			#f calcule <wx>+b pour chaque valeur du vecteur
			f<-function(v){
				return(rowSums(as.matrix(v)%*%t(as.matrix(w)))+b)
			}
			#pour retourner la fonction pred dans le return, on cree le vecteur reponse
			h=numeric() 
			#le return consiste a renvoyer un vecteur contenant lui-meme les definitions de fonctions que l'on souhaite retourner mais qui sont parametrees par les resultats internes de mediatorHyperplane()
			return(
				c(
					w,
					b,
					#fonction f
					function(v){return(rowSums(as.matrix(v)%*%t(as.matrix(w)))+b)},
					#fonction pred
					function(v){
						for (i in f(v)) {
							if (i>0) {
								h=append(h,1,length(h))
							} 
							else {
								h=append(h,-1,length(h))
							}
						}
						return(h)
					}
				)
			)
		}

resultat.mediatorHyperplane.Simple<-mediatorHyperplane(dataset.generateSimpleDataset)#
resultat.mediatorHyperplane.Difficult<-mediatorHyperplane(dataset.generateDifficultDataset)#
f.Simple<-resultat.mediatorHyperplane.Simple[4][[1]]#
pred.Simple<-resultat.mediatorHyperplane.Simple[5][[1]]#
f.Difficult<-resultat.mediatorHyperplane.Difficult[4][[1]]#
pred.Difficult<-resultat.mediatorHyperplane.Difficult[5][[1]]#

########################################################################
#	4
########################################################################

#affiche la table de confusion
#les lignes correspondent à la realité et les colonnes aux predictions
matriceconfusion<-table(dataset.generateDifficultDataset$y,pred.Difficult(dataset.generateDifficultDataset[,1:2]))

#taux d'erreur
calcul.tauxerreur<-function(matricedeconfusion){
	return((matricedeconfusion[1,2]+matricedeconfusion[2,1])/(sum(rowSums(matricedeconfusion))))
}

#sensibilite
calcul.sensibilite<-function(matricedeconfusion){
	return(matricedeconfusion[1,1]/(matricedeconfusion[1,1]+matricedeconfusion[2,1]))
}

#specificite
calcul.specifite<-function(matricedeconfusion){
	return(matricedeconfusion[2,2]/(matricedeconfusion[2,2]+matricedeconfusion[1,2]))
}

#précision
calcul.precision<-function(matricedeconfusion){
	return(matricedeconfusion[1,1]/(matricedeconfusion[1,1]+matricedeconfusion[1,2]))
}

calcul.tauxerreur(matriceconfusion)#
calcul.sensibilite(matriceconfusion)#
calcul.specificite(matriceconfusion)#
calcul.precision(matriceconfusion)#

library(ROCR)

dataset.generateDifficultDataset<-generateDifficultDataset(100)
predictor <- prediction(pred.Difficult(dataset.generateDifficultDataset[,1:2]),dataset.generateDifficultDataset$y)
matriceconfusion<-table(dataset.generateDifficultDataset$y,pred.Difficult(dataset.generateDifficultDataset[,1:2]))
#taux d'erreur selon ROC
(perfector <- performance(predictor, "err")@y.values[[1]][2])
#taux d'erreur selon nous
calcul.tauxerreur(matriceconfusion)
#sensibilite selon ROC
(perfector <- performance(predictor, "sens")@y.values[[1]][2])
#sensibilite selon nous
calcul.sensibilite(matriceconfusion)
#specificite selon ROC
(perfector <- performance(predictor, "spec")@y.values[[1]][2])
#specificite selon nous
calcul.specifite(matriceconfusion)
#precision selon ROC
(perfector <- performance(predictor, "prec")@y.values[[1]][2])
#precision selon nous
calcul.precision(matriceconfusion)

#on trace la courbe ROC
dataset.generateDifficultDataset<-generateDifficultDataset(100)#
predictor <- prediction(f.Difficult(dataset.generateDifficultDataset[,1:2]),generateDifficultDataset(100)$y)
perfector <- performance(predictor, measure = "sens", x.measure = "spec") 
plot(perfector)

CrossValidation <- function(mediator,d,N){
	tauxderreur<-c()
	d$nom<-rownames(d)
	#on melange le tableau
	e <- d[sample(1:nrow(d),nrow(d)),]
	#si congrue a 0 modulo N
	if (nrow(e)%%N==0){
		for (i in 1:N){
			test<-e[((((i-1)/N)*nrow(e))+1):((i/N)*nrow(e)),]
			#train est le complementaire de test
			train<-e[e$nom!=test$nom,]
			#on parametre le f avec le train
			pred.train<-mediator(train[,1:3])[5][[1]]
			#on cree une matrice de confusion
			matriceconfusion<-table(test$y,pred.train(test[,1:2]))
			print(matriceconfusion)#
			tauxderreur<-c(tauxderreur,calcul.tauxerreur(matriceconfusion))
			#print(i)#
			#print(tauxderreur)#
		}
	}
	#si congrue a autre chose que 0 modulo N
	else{
		for (i in 1:N){
			nblignesaretirer<-nrow(e)%%N
			#on divise test d'une autre maniere
			if (i==N){
				test<-e[((((i-1)/N)*(nrow(e)-nblignesaretirer))+1):((i/N)*(nrow(e))),]
				print(test)#
			}
			else{
				test<-e[((((i-1)/N)*(nrow(e)-nblignesaretirer))+1):((i/N)*(nrow(e)-nblignesaretirer)),]
				print(test)#
			}
			#train est le complementaire de test
			train<-e[e$nom!=test$nom,]
			#on parametre le f avec le train
			pred.train<-mediator(train[,1:3])[5][[1]]
			print(pred.train(test[,1:2]))#
			#on cree une matrice de confusion
			matriceconfusion<-table(test$y,pred.train(test[,1:2]))
			print(matriceconfusion)#
			tauxderreur<-c(tauxderreur,calcul.tauxerreur(matriceconfusion))
		}
	}
	return(mean(tauxderreur))
}

(cv.congru0=CrossValidation(mediatorHyperplane,dataset.generateDifficultDataset,5))
(cv.pascongrue0=CrossValidation(mediatorHyperplane,dataset.generateDifficultDataset,3))

CrossValidation(mediatorHyperplane,generateDifficultDataset(40),5)
CrossValidation(mediatorHyperplane,generateDifficultDataset(50),5)
CrossValidation(mediatorHyperplane,generateDifficultDataset(60),5)
CrossValidation(mediatorHyperplane,generateDifficultDataset(100),5)
CrossValidation(mediatorHyperplane,generateDifficultDataset(1000),5)
CrossValidation(mediatorHyperplane,generateDifficultDataset(10000),5)
CrossValidation(mediatorHyperplane,generateDifficultDataset(100000),5)
#1000000 prend du temps...
#CrossValidation(mediatorHyperplane,generateDifficultDataset(1000000),10)#
