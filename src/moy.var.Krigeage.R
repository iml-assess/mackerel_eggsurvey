moy.var.Krigeage<-function(xCoords, xValues, kCoords, vario){
	#--------------------------------------------------------------------------#
	# Merci à Denis Marcotte pour sont aide dans la définition des formules:
	#	La variance d'estimation globale est: C(v,v)-lambda_v'*C(x,v)
	#		C(v,v) 		est la moyenne de toutes les covariances entre les 
	#					p(x0,y0) points krigés;
	#		lambda_v 	est le vecteur de poids du krigeage global obtenu de:
	#					C(x,x) * lambda_v = C(x,v)
	#		C(x,x) 		est la matrice contenant les covariances entre les 
	#					points de données brutes puis compléter par la rangée 
	#					et la colonne de 1 pour la contrainte de la somme des 
	#					poids égale à 1)
	#		C(x,v) 		est le vecteur contenant la moyenne des covariances 
	#					entre chacun des points de données brutes et les points 
	#					krigés. Ce vecteur est complété par 1 pour la contrainte.
	#		Z			Moyenne de krigeage
	#		VarKrigeage	Variance globale de krigeage
	#--------------------------------------------------------------------------#
	require(geoR)
	
	# Vérification
	if(!is.matrix(xCoords) | ncol(xCoords)!=2) stop("xCoords doit etre une matrice de deux colonnes!\n")
	if(!is.matrix(kCoords) | ncol(kCoords)!=2) stop("kCoords doit etre une matrice de deux colonnes!\n")
	if(class(vario)[1]!="variogramModel") stop("vario doit etre de la classe variogramModel!\n")
	if(!is.vector(xValues) | length(xValues)!=nrow(xCoords)) stop("xValues is not a vector or is not the same length as xCoords!\n")

	nb.x=nrow(xCoords)
	nb.k=nrow(kCoords)
	
	# Extraction des données du variogramme
	if(nrow(vario)==1){
		if(vario[1,1]=="Nug") model="pure.nugget"			
		else if(vario[1,1]=="Exp") model="exponential"			
		else if(vario[1,1]=="Sph") model="spherical"			
		else if(vario[1,1]=="Gau") model="gaussian"			
		else if(vario[1,1]=="Mat") model="matern"			
		else if(vario[1,1]=="Cir") model="circular"			
		else if(vario[1,1]=="Pow") model="power"
		else(stop("Le modele fournit n'est pas reconnu!\n"))
		pars=c(vario[1,2],vario[1,3])
		nugget=0
	}else if(nrow(vario)==2){
		tmp=which(vario[,1]=="Nug")
		if(length(tmp)==0) 	stop("Le modele fournit n'est pas reconnu!\n")
		nugget=vario[tmp,2]
		vario=vario[-tmp,]
		if(vario[1,1]=="Nug") model="pure.nugget"			
		else if(vario[1,1]=="Exp") model="exponential"			
		else if(vario[1,1]=="Sph") model="spherical"			
		else if(vario[1,1]=="Gau") model="gaussian"			
		else if(vario[1,1]=="Mat") model="matern"			
		else if(vario[1,1]=="Cir") model="circular"			
		else if(vario[1,1]=="Pow") model="power"
		else( stop("Le modele fournit n'est pas reconnu!\n"))
		pars=c(vario[1,2],vario[1,3])
	}else{ 
		stop("vario model can have only two rows!\n")
	}
	
	# Matrice de covariance entre les points de données brutes
	Cxx=varcov.spatial(xCoords, cov.model=model, cov.pars=pars, nugget=nugget)
	Cxx=cbind(Cxx$varcov, rep(1, nb.x))
	Cxx=rbind(Cxx,c(rep(1,nb.x),0))
	
	# Matrice de covariance entre les points krigés
	Cx0x0=varcov.spatial(kCoords, cov.model=model, cov.pars=pars, nugget=nugget)
	
	# Matrice de covariance entre les données et le points krigés
	Cxx0=varcov.spatial(cbind(c(xCoords[,1],kCoords[,1]),c(xCoords[,2],kCoords[,2])), 
	cov.model=model, cov.pars=pars, nugget=nugget)
	Cxx0=t(Cxx0$varcov[c(1:nb.k)+nb.x, c(1:nb.x)])

	# Système de krigeage global
	# C(v,v)
	Cvv=mean(Cx0x0$varcov)
	# C(x,v)
	Cxv=apply(Cxx0,1,mean)
	Cxv=matrix(c(Cxv,1),ncol=1)
	
	# Résolution du système de krigeage
	# Poids
	lambda_v=solve(Cxx,Cxv, tol=0.000000000000000000000000000000000000000001)
	
	# Estimation de la variance de krigeage globale
	VarKrigeage= Cvv-t(lambda_v)%*%Cxv

	# Estimation de la moyenne de krigeage globale
	Z=sum(xValues*lambda_v[1:(length(lambda_v)-1)])
	
	return(list(Z=Z,VarKrigeage=VarKrigeage))
}
	