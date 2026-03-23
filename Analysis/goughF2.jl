using RCall, StatsBase, Statistics, DelimitedFiles, Distributed,LinearAlgebra
#generating workers
addprocs(length(Sys.cpu_info())-3)

@everywhere using FlxQTL
setSeed(123)
##readding data
phen=readdlm("../processedData/goughF2_sex_imp_16weight.csv",',';skipstart=1)
## get marker info excluding X
minfo=readdlm("../processedData/goughF2_markers.csv",',';skipstart=1);
midx=findall(minfo[:,2].!="X");
mar=minfo[midx,:]

#traits
wght=readdlm("../processedData/goughF2_16wk_weight.txt",'\t')

#genotype probability data excluding X chromosome (1212 x 12777 markers (38331 prob))
genpr=readdlm("../processedData/goughF2_imp_genoprob.csv",',';skipstart=1);

###generate 4 df of cubic spline bases
R"library(splines)"
R"Ns4<-bs(c(1:16),df=4,intercept=T)"
@rget Ns4

#log transformation for traits
logY= log.(wght')
m,n=size(logY)

XX=Markers(mar[:,1],mar[:,2],mar[:,3],genpr') # for flxMLMM scan

#kinship
@time K=kinshipLoco(kinshipLin,XX,3);
@time T, Λ=K2eig(K,true);

#using bs= 4 df 
# multivariate linear (mixed) model:  flxMLMM vs flxMLM
#flxMLMM
@time Lod4,B4,est0= geneScan(3,T,Λ,logY,XX,Ns4,true;df_prior=Int64(ceil(1.9m)), Prior=3cov(logY,dims=2));  #H0: loglik estimated by MLMM
# flxMLM
X1= Markers(mar[:,1],mar[:,2],mar[:,3],genpr);  lnY=log.(wght);
@time mlmlod, mlmb,mlmns04=mlm1Scan(3,lnY,X1,Ns4);

open("../Result/goughF2_ns4_mlmm_mlm.txt","w")do io
    writedlm(io,[Lod4 mlmlod])
end

#### Permutation Tests
#permutation test for flxMLMM (LOCO supported) to get cutoff values
@time maxLODs,H1par_perm,cutoff=permutationTest(500,3,K,logY,XX;Z=Ns4,df_prior=Int64(ceil(1.9m)),Prior=3cov(logY,dims=2),LOCO_all=true);

writedlm("../Result/goughF2_maxlods_mlmm.txt",maxLODs)

f=open("../Result/goughF2_maxlods_mlmm.txt","a") #another 500 maxlods
writedlm(f,maxLODs)
close(f)

# permutation test for flxMLM
@time mmlod,h1permlm,mlm_cut=mlmTest(1000,3,lnY,X1,Ns4);
writedlm("../Result/goughF2_maxlods_mlm_$(mlm_cut).txt",mmlod)


############
## 1-d plot : Fig8.eps
# Lod4=readdlm("../Result/goughf2_lod_z(4).txt")
p=layers(mar[:,2],mar[:,3],Lod4)
plot1d(p;title="goughF2 with α=0.1,0.05",yint=[6.3869,7.7035])


##2d-scan

geno2d=readdlm("../processedData/goughF2_genoprob_2d_step0_auto.csv",',';skipstart=1);
marker2d=readdlm("../../fda/processedData/goughF2_markerinfo_2d_step0.csv",',';skipstart=1);
idx_2d=findall(marker2d[:,2].!="X");
mar2d=marker2d[idx_2d,:]

## dropping markers whose distance is less than 0.25 (filtering)
incl_idx=findall(diff(mar2d[:,3]).>=0.25)

#selecting markers
X0=mat2array(3,geno2d')
X0=X0[:,:,incl_idx]
X1=array2mat(3,X0)

#filtered marker info
mar2=mar2d[incl_idx,:]    
    
X2d=Markers(mar2[:,1],mar2[:,2],mar2[:,3],X1)  # for flxMLMM
X2D=Markers(mar2[:,1],mar2[:,2],mar2[:,3],X1') #  for flxMLM
# @time l2d,est2d=FlxQTL.gene2Scan(3,T,Matrix(1.0I,m,m),Λ,ones(m),Ystd,X2d,Z_b4,true);    

## 2-d scan for three chromosomes
# chr7=findall(X2d.chr.==7);X7=FlxQTL.Markers(mar2[chr7,1],mar2[chr7,2],mar2[chr7,3],X1[3*chr7[1]-2:chr7[end]*3,:])
# chr8=findall(X2d.chr.==8);X8=FlxQTL.Markers(mar2[chr8,1],mar2[chr8,2],mar2[chr8,3],X1[3*chr8[1]-2:chr8[end]*3,:])
# chr10=findall(X2d.chr.==10);X10=FlxQTL.Markers(mar2[chr10,1],mar2[chr10,2],mar2[chr10,3],X1[3*chr10[1]-2:chr10[end]*3,:])

@time mlmlod2,est02=mlm2Scan(3,lnY,X2D,Ns4); # for flxMLM
@time lod2d,es02d=gene2Scan(3,T,Λ,logY,X2d,true;Z=Ns4,df_prior=Int64(ceil(1.9m)),Prior=3cov(logY,dims=2)); # for flxMLMM
# saved in Sen's computer but also copied to /fda/analysis/result/


open("../Result/goughF2_2D_logY_ns4_mlm.txt","w")do io
    writedlm(io,mlmlod2)
end

open("../Result/goughF2_2D_logY_ns4_mlmm.txt","w")do io
    writedlm(io,lod2d)
end

###########
## 2d-plots for the three Chromosomes: Fig9 (three plots)
# ch7=readdlm("../Result/goughf2_2dlod_chr7.txt");
# ch8=readdlm("../Result/goughf2_2dlod_chr8.txt");
# ch10=readdlm("../Result/goughf2_2dlod_chr10.txt");

p7=layers(X7.chr,X7.pos,ch7)
p8=layers(X8.chr,X8.pos,ch8)
p10=layers(X10.chr,X10.pos,ch10)

plot2d(p7)
plot2d(p8)
plot2d(p10)

# ## effect plots for particular markers in three chromosomes
# using Gadfly

# # find markers that need to be inspected based on the paper
# #findall(c7.==maximum(c7[7:20,55:80])) the next bigest lod (11,55)
# a1=findall(c7.==maximum(c7)) #(59,57)
# a2=findall(c8.==maximum(c8)) #(40,14)
# a3=findall(c10.==maximum(c10)) #(84,22)
# ######
# ch71=findall(mar[:,3].==X7.pos[11]);ch72=findall(mar[:,3].==X7.pos[59])
# ch81=findall(mar[:,3].==X8.pos[14]);ch82=findall(mar[:,3].==X8.pos[40])
# ch101=findall(mar[:,3].==X10.pos[22]);ch102=findall(mar[:,3].==X10.pos[84])
# #####

#################
##### For plots not to run all code
# chr7=findall(mar2[:,2].==7);
# chr8=findall(mar2[:,2].==8)
# chr10=findall(mar2[:,2].==10)
# ch71=findall(mar[:,3].==mar2[chr7,end][11]);ch72=findall(mar[:,3].==mar2[chr7,end][59])
# ch81=findall(mar[:,3].==mar2[chr8,end][14]);ch82=findall(mar[:,3].==mar2[chr8,end][40])
# ch101=findall(mar[:,3].==mar2[chr10,end][22]);ch102=findall(mar[:,3].==mar2[chr10,end][84])


# # genotype data for effect plots
# gene=readdlm("../processedData/goughF2_imp_genotypes_auto.csv",',';skipstart=1); #genotypes = c("GG","GW","WW") =(1,2,3)
# maxidx=[ch71;ch72;ch81;ch82;ch101;ch102]

# gg=zeros(16,length(maxidx));gw=copy(gg);ww=copy(gg);
# for k=1:length(maxidx)
#  gg[:,k]=  mean(wght[gene[:,maxidx[k]].==1.0,:],dims=1)' 
#  gw[:,k]=  mean(wght[gene[:,maxidx[k]].==2.0,:],dims=1)'  
#  ww[:,k]=  mean(wght[gene[:,maxidx[k]].==3.0,:],dims=1)' 
# end

# ## Fig_S11 ~ Fig_S13
# #generate effect plots from raw data
# Plots=Array{Plot,1}(undef,length(maxidx));ticks=[Float64(k) for k=1:16]
# for j=1:length(maxidx)
#     Plots[j] =Gadfly.plot(layer(y=gg[:,j],Geom.point,Geom.line,
#                 Theme(default_color=colorant"green")),
#             layer(y=gw[:,j],Theme(default_color=colorant"purple"),
#            Geom.point,Geom.line),layer(y=ww[:,j],Geom.point,Geom.line),
#         Coord.cartesian(ymin=floor(minimum([gg gw ww])), ymax=ceil(maximum([gg gw ww]))),Guide.XLabel("t (week)"),Guide.xticks(ticks=ticks)
#         ,Guide.title("Chr$(mar[maxidx[j],2])@$(mar[maxidx[j],3])")
#         ,Guide.YLabel("Body weight",orientation=:vertical)
#         ,Guide.manual_color_key("Genotypes",["gg","gw","ww" ], ["green","purple","skyblue" ]) 
#         );
    
# end

# #subtract overall mean function (in a sense of fda)

# overall=mean(wght,dims=1)'
# Plot1=Array{Plot,1}(undef,length(maxidx));
# for j=1:length(maxidx)
#     Plot1[j] =Gadfly.plot(layer(y=gg[:,j].-overall,Geom.point,Geom.line,
#                 Theme(default_color=colorant"green")),
#             layer(y=gw[:,j].-overall,Theme(default_color=colorant"purple"),
#            Geom.point,Geom.line),layer(y=ww[:,j].-overall,Geom.point,Geom.line),
#         Coord.cartesian(ymin=floor(minimum([gg gw ww].-overall)), ymax=ceil(maximum([gg gw ww].-overall)))
#         ,Guide.title("Body weight-E(Body weight(t))")
#         ,Guide.YLabel(" ",orientation=:vertical)
#         ,Guide.XLabel("t (week)"),Guide.xticks(ticks=ticks)
#         ,Guide.manual_color_key("Genotypes",["gg","gw","ww" ], ["green","purple","skyblue" ]));
    
# end

# set_default_plot_size(25cm, 20cm)
# display(gridstack([Plots[1] Plots[2];Plot1[1] Plot1[2]]))
# display(gridstack([Plots[3] Plots[4];Plot1[3] Plot1[4]]))
# display(gridstack([Plots[5] Plots[6];Plot1[5] Plot1[6]]))



