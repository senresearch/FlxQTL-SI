using Statistics
@everywhere using FlxQTL



####### generating Sigma 
# a=0.4; b=0.07  (default)
# A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);B=b*ones(3,3);
# Σtrue=[ A B;B A] 

######### generate extreme cases of Kc
#generate off Kc for model misspecification (r<0.2 for pdf)
function H(r,m)
    (1+r)*Matrix(1.0I,m,m)-r*ones(m,m)
end


#generate AR Kc
function Kc_ar(ρ,m)
   K=zeros(m,m)
    for j=1:m
        for i=1:m
            K[j,i]=ρ^abs(i-j)
        end
    end
    
    return K
end
################

##check the variance Ysim0
function checkVar(n,m,nchr,itr,Ysim0,K,Kc,τ2true,Σtrue)
     y=zeros(n*m,itr);
     for j=1:itr
       y[:,j]=vec(Ysim0[:,n*(j-1)+1:n*(j-1)+n])
     end
# relative error
    for k=1:nchr
     println(norm(cov(y')-(kron(τ2true*Kc,K[:,:,k])+kron(Σtrue,Matrix(1.0I,n,n))))/
       norm(kron(τ2true*Kc,K[:,:,k])+kron(Σtrue,Matrix(1.0I,n,n))))
    end
end

#compute type I, II errors
function tErrors(itr,LOD0,LOD1,idxTrue;erate=[95,99])  
    maxlod0=zeros(itr);idxlrg=copy(maxlod0)
      for i=1:itr
     maxlod0[i]=maximum(LOD0[:,i])
      end
     α=percentile(maxlod0,erate)
     power=zeros(length(idxTrue))
     for j=1:length(idxTrue)
     power[j] = mean((LOD1[idxTrue[j],:].>α[1]))
      end
   return α,power
end


#find a median of each parameter estimation
function Parameters(nchr::Int64,itr,par_est)
        ##get τ2
Τ0=zeros(nchr,itr); tau1=[]
for j=1:nchr*itr
tau=par_est[j].τ2
tau1=[tau1;tau]
end
for i=1:itr
   Τ0[:,i]= tau1[nchr*(i-1)+1:nchr*(i-1)+nchr]
end

#get  Σ   
m=size(par_est[1].Σ,1)
Σ0=zeros(m,m,nchr,itr);
    Sig=par_est[1].Σ
for i=2:nchr*itr
  Sig=[Sig;par_est[i].Σ]
end
      
for j=1:itr
    for i=1:nchr
    sig1=Sig[(nchr*m)*(j-1)+1:(nchr*m)*(j-1)+nchr*m,:]
     Σ0[:,:,i,j]=sig1[m*(i-1)+1:m*(i-1)+m,:]  
    end
end
    return mean(Τ0,dims=2),median(Τ0,dims=2), mean(Σ0,dims=4),median(Σ0,dims=4)
end


# cr : Int64, cross (1 for genotypes, >1 type of crosses)
# X_sml : # of small effect qtl
# X_lrg : # of large effct qtl
#q= size(Z,2): mxq
# bs : a vector for sampling B
function Bgenerator(cr::Int64,X_lrg::Int64,q::Int64,bs::Array{Float64,1})
    
   if(cr!=1)
         p=cr
    else
         p=cr+1
    end
        
     Bfix=zeros(q,p,X_lrg)
     B=reshape(sample(bs,(p-1)*q),q,:)
    display(B) # to check the Bfix
     for i=1:X_lrg
#             Bfix[:,:,i]= hcat(zeros(q),B.^(X_lrg-i))
             Bfix[:,:,i]= hcat(zeros(q),B.^i)
    end
        return Bfix
    
end


# a random generator  of qtl large effect indices
function Lqtl_index(X_lrg::Int64,p::Int64)

    Lqidx=sample(1:p,X_lrg;replace=false);sort!(Lqidx)   
        return Lqidx
end


# generative Y

# w : # of weather covariates
# markers :standardized for genotypes only / Z: normalized
# Climate: a climate matrix (or data)
function Ygenerator(itr,τ2true,Σtrue,Bfix,cr,Z,Lqidx::Array{Int64,1},X_sml::Int64,w::Int64,Climate,markers)

m=size(Z,1);q=size(Z,2);n=size(markers.X,2);
p=Int(size(markers.X,1)/cr); q1=size(Climate,1)  #C : a q1 x m climate data
X_lrg=length(Lqidx)
Z=Z./norm(Z)
# generating fixed effects (large effects)
#Lqidx=sample(1:p,X_lrg;replace=false);sort!(Lqidx)  
Sqidx=sample(setdiff(1:p,Lqidx),X_sml;replace=false);sort!(Sqidx)
   
        Yfix=zeros(m,n,X_lrg)
    if(cr!=1)
        X=FlxQTL.mat2array(markers.X,cr)
            for j=1:X_lrg
                Yfix[:,:,j]=(Z*Bfix[:,:,j])*X[Lqidx[j],:,:]
            end
    else
        XL=markers.X[Lqidx,:]
        if(length(Lqidx)==1)
            Xstd=(XL.-mean(XL))./std(XL)
        else
            Xstd=(XL.-mean(XL,dims=2))./std(XL,dims=2)
        end
            for j=1:X_lrg
              Yfix[:,:,j]=(Z*Bfix[:,:,j])*vcat(ones(1,n),Xstd[[j],:])
            end
    end
        Y=sum(Yfix,dims=3)

#generating small effects and weather indices
    if (Climate!=Matrix(1.0I,m,m))
    widx=sample(1:q1,w;replace=false);sort!(widx)
    C=Climate[widx,:]
    else
    C=Climate
    end 

b=MvNormal(τ2true*Matrix(1.0I,w*X_sml*cr,w*X_sml*cr));       
e=MvNormal(kron(Matrix(1.0I,n,n),Σtrue));
#####
    if (cr!=1)
      Xs=FlxQTL.array2mat(cr,X[Sqidx,:,:])                                            
      Y_sim= @distributed (hcat) for i=1:itr
         Brand=reshape(rand(b),w,:)/sqrt(w*X_sml*cr)
         Yrand=(C'*Brand)*Xs
         E=reshape(rand(e),m,:)/sqrt(m*n)
    ##generate a phenotype matrix
         Y0=(Yrand.-mean(Yrand,dims=2))./std(Yrand,dims=2)+(E.-mean(E,dims=2))./std(E,dims=2)
         Y1=Y+Y0
     
         [Y0;Y1]
                          end
    else #cr=1
      Xs=markers.X[Sqidx,:]
      Y_sim= @distributed (hcat) for i=1:itr
         Brand=reshape(rand(b),w,:)/sqrt(w*X_sml*cr)
         Yrand=(C'*Brand)*Xs
         E=reshape(rand(e),m,:)/sqrt(m*n)
       ##generate a phenotype matrix
         Y0=(Yrand.-mean(Yrand,dims=2))./std(Yrand,dims=2)+(E.-mean(E,dims=2))./std(E,dims=2)
         Y1=Y+Y0
                        
          [Y0;Y1]
                          end                        
    end                          

     Ysim0=Y_sim[1:m,:]  # H0
    Ysim1=Y_sim[m+1:end,:] # H1
   
 return  Ysim0, Ysim1
   
end
    


#direct sampling Y w/o Y
    
function Ygen(itr,τ2true::Float64,Σtrue::Array{Float64,2},Bfix::Array{Float64,3},cr,X_lrg::Int64,Kg,Kc,XX::FlxQTL.Markers)
        
n=size(XX.X,2);m=size(Kc,1);
p=Int(size(XX.X,1)/cr);
 Chr0=unique(XX.chr);

# generating fixed terms
     #get large qtl indices   
        Chr=sample(Chr0,X_lrg;replace=false);sort!(Chr)
        println("large qtl are assinged to $(Chr).")
        Qidx=[]
        for j=1:X_lrg
            maridx=findall(XX.chr.==Chr[j])
            qidx=sample(maridx,1)
            Qidx=[Qidx;qidx]
        end
         println("Large effect QTL indices are $(Qidx).")   
        
  Yfix=zeros(m,n,X_lrg)
    if(cr!=1) # genotype probability
        X=FlxQTL.mat2array(XX.X,cr)
            for j=1:X_lrg
                Yfix[:,:,j]=(Bfix[:,:,j])*X[Qidx[j],:,:]
            end
        else # genotypes
        XL=XX.X[Qidx,:]
                
#         if(X_lrg==1)
#             Xstd=(XL.-mean(XL))./std(XL)
#         else
#             Xstd=(XL.-mean(XL,dims=2))./std(XL,dims=2)
#         end
            for j=1:X_lrg
              
              Yfix[:,:,j]=(Bfix[:,:,j])*vcat(ones(1,n),XL[[j],:])
            end
    end
        Y=sum(Yfix,dims=3)[:,:,1];  #Y=(Y.-mean(Y,dims=2))./std(Y,dims=2)

#direct random terms sampling
b=MvNormal(kron(Kg,τ2true*Kc));       
e=MvNormal(kron(Matrix(1.0I,n,n),Σtrue));
#####
    
      Y_sim= @distributed (hcat) for i=1:itr
         Yrand=reshape(rand(b),m,:)/sqrt(m*n)
         E=reshape(rand(e),m,:)/sqrt(m*n)
       ##generate a phenotype matrix
#          Y0=(Yrand.-mean(Yrand,dims=2))./std(Yrand,dims=2)+(E.-mean(E,dims=2))./std(E,dims=2)
          Y0= Yrand+E
          Y1=Y+Y0
              
          [Y0;Y1]
                          end                                                 

     Ysim0=Y_sim[1:m,:]  # H0
    Ysim1=Y_sim[m+1:end,:] # H1
   
 return Qidx, Ysim0, Ysim1
   
        
end
    
    
    

function simulation(itr,cr,Lqidx,Y_sim0,Y_sim1,Tg,Tc,λg,λc,markers,Z;ρ=0.001,itol=1e-4,tol=1e-4)
    est=[];q=size(Z,2); p,n=size(markers.X); X_lrg=length(Lqidx)
LOD0=zeros(p,itr); LOD1=copy(LOD0);
if (cr==1)
    p1=cr+1
else
    p1=cr   
end
    Bs=zeros(q,p1,X_lrg,itr); T=zeros(itr,2);
    
    for t=1:itr
    println("Start $(t)th simulation.")
    println("genome scan with small effects(H0).")
t0=@elapsed lod0,B0,est0=FlxQTL.geneScan(cr,Tg,Tc,λg,λc,Y_sim0[:,n*(t-1)+1:n*(t-1)+n] ,markers,Z,true;ρ=ρ,itol=itol,tol0=tol,tol=tol)
    println("genome scan with large effects(H1).")                
t1=@elapsed lod1,B1,est01=FlxQTL.geneScan(cr,Tg,Tc,λg,λc,Y_sim1[:,n*(t-1)+1:n*(t-1)+n],markers,Z,true;ρ=ρ,itol=itol,tol0=tol,tol=tol)
    LOD0[:,t]=lod0;LOD1[:,t]=lod1; est=[est;est0];T[t,:]=[t0 t1]
         for l=1:X_lrg                                           
                Bs[:,:,l,t]=B1[:,:,Lqidx[l]]
         end
    end

     return LOD0, LOD1, Bs, est, median(T,dims=1)   
end

#MVLMM
function simulation(itr,cr,Lqidx,Y_sim0,Y_sim1,Tg,λg,markers,Z;ρ=0.001,itol=1e-4,tol=1e-4)
    est=[];q=size(Z,2); p,n=size(markers.X); X_lrg=length(Lqidx)
LOD0=zeros(p,itr); LOD1=copy(LOD0);
if (cr==1)
    p1=cr+1
else
    p1=cr   
end
    Bs=zeros(q,p1,X_lrg,itr); T=zeros(itr,2);
    
    for t=1:itr
    println("Start $(t)th simulation.")
    println("genome scan with small effects(H0).")
t0=@elapsed lod0,B0,est0=FlxQTL.geneScan(cr,Tg,λg,Y_sim0[:,n*(t-1)+1:n*(t-1)+n],markers,Z,true;itol=itol,tol0=tol,tol=tol,ρ=ρ)
    println("genome scan with large effects(H1).")                
t1=@elapsed lod1,B1,est01=FlxQTL.geneScan(cr,Tg,λg,Y_sim1[:,n*(t-1)+1:n*(t-1)+n],markers,Z,true;itol=itol,tol0=tol,tol=tol,ρ=ρ)
    LOD0[:,t]=lod0;LOD1[:,t]=lod1; est=[est;est0];T[t,:]=[t0 t1]
         for l=1:X_lrg                                           
                Bs[:,:,l,t]=B1[:,:,Lqidx[l]]
         end
    end

     return LOD0, LOD1, Bs, est, median(T,dims=1)   
end

