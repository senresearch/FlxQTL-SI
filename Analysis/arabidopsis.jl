using DelimitedFiles, StatsBase, Statistics, LinearAlgebra, Distributed
#generating workers
addprocs(length(Sys.cpu_info())-3)

@everywhere using FlxQTL
setSeed(123)

## read imputed genotype data: a(=1):italian parent, b(=2)swedish parent
impgen = readdlm("../processedData/fullrank_imput.csv",',';skipstart=1);
impgen[impgen.==1.0].=0.0;impgen[impgen.==2.0].=1.0;
impgen1=copy(impgen);
impgen1[impgen1.==0.0].=-1.0;

#genotype data for 1D scan
labels=readdlm(string(@__DIR__,"/../processedData/marlabels_agren.csv"),',';skipstart=1);
XX=Markers(labels[:,1],labels[:,2],labels[:,3],impgen1');

## read phenotype data
pheno =readdlm(string(@__DIR__,"/../processedData/pheno2013_imp.csv"),',';header=true);
pheno=pheno[1][:,2:end-1];
Y=convert(Array{Float64,2},pheno');
Ystd=(Y.-mean(Y,dims=2))./std(Y,dims=2);

#site contrasts: italy (-1), sweden (1)
Z=hcat(ones(6),vcat(-ones(3),ones(3)));
m,q=size(Z)

#Kinship
Xkin=Markers(labels[:,1],labels[:,2],labels[:,3],impgen');
@time K=shrinkgLoco(kinshipMan,30,Xkin);
@time T, Λ=K2eig(K,true);


#multivariate genome scans
@time LODs,Bs,est0=geneScan(1,T,Λ,Ystd,XX,Z,true); 

open("../Result/agren_1dlod.txt","w")do io
    writedlm(io,LODs)
end

### rearrange effects data 
main=zeros(699,1); inter=copy(main); 
for l=1:699
    main[l]=Bs[1,2,l] 
    inter[l]=Bs[2,2,l]
end
#store effects in a file
open("../Result/agren_effects_main_inter_comparison.txt", "w") do io
           writedlm(io, [main inter])
       end


#italy
@time LODi,bi,est0i=geneScan(1,T,Λ,Ystd[1:3,:],XX,true);

#Sweden
@time LODw,bw,est0w=geneScan(1,T,Λ,Ystd[4:end,:],XX,true);

open("../Result/agren_sitewise_1dlod_it_sw_Z=I.txt","w") do io
    writedlm(io,[LODi LODw])
end

#run univariate genomescan
uni_lods=zeros(699,m);
for l=1:m
  @time uni_lods[:,l],_,_=geneScan(1,T,Λ,Ystd[[l],:],XX,true)
end

open("../Result/agren_univariate_1dlod_it_sw.txt","w")do io
    writedlm(io,uni_lods)
end

##2d-scan
labels1=readdlm("../processedData/agren_newmarkers_multipleqtl.csv",',';skipstart=1);
G1=readdlm("../processedData/agren_gene_imp_revision.csv",',';skipstart=1);
Xp=FlxQTL.Markers(labels1[:,1],labels1[:,2],labels1[:,3],G1');
#@time K2=shrinkgLoco(kinshipMan,100,Xp); ## sec
#T, Λ=EcmNestrvQTL.K2eig(K2,true);

@time LOD2,est02=gene2Scan(1,T,Λ,Ystd,Xp,true;Z=Z);
open("../Result/agren_2dlod.txt","w")do io
    writedlm(io,LOD2)
end


##permutation

@time MaxLod,H1par,cutoff=permutationTest(1000,1,K,Ystd,XX;Z=Z);
writedlm("../Result/agren-locoperm-maxlods.txt",MaxLod)

# ## Rearrange effects from permutations to get 95% band 
# bs=@distributed (vcat) for l=1:1000*699
#         [H1par[l].B[1,2] H1par[l].B[2,2]]
#  end    

# Bmain_per=zeros(699,1000); Binter_per=zeros(699,1000);
# for k=1:1000
#      Bmain_per[:,k]=bs[699*(k-1)+1:699*(k-1)+699,1]
#      Binter_per[:,k]=bs[699*(k-1)+1:699*(k-1)+699,2]
# end

# # two-sided bands for main and interaction effects, respectively
# perm_eff=Array{Float64,2}(undef,699,4);
# for i=1:699
#     perm_eff[i,1]=percentile(Bmain_per[i,:],2.5)
#     perm_eff[i,2]=percentile(Bmain_per[i,:],97.5)
#     perm_eff[i,3]=percentile(Binter_per[i,:],2.5)
#     perm_eff[i,4]=percentile(Binter_per[i,:],97.5)
# end

####### Effects plots
#using Gadfly
include("effectplot.jl")

#Fig_S9.png (main), Fig_S10.png (interaction)
pmain_perm=Mlayers(XX.chr,XX.pos,hcat(main[:,1],perm_eff[:,1:2]))
pint_perm=Mlayers(XX.chr,XX.pos,hcat(inter[:,1],perm_eff[:,3:4]))

effectplot(pmain_perm,true;title="Main effect with 95% band", ylabel="effect size") 
effectplot(pint_perm,true;title="Interaction effect with 95% band", ylabel="effect size") 

# peff=Mlayers(XX.chr,XX.pos,hcat(main[:,1],inter[:,1]))
# effect2plot(peff;title="Effect Plots",ylabel="Effect size",legend=["Main","Interaction"])

### Fig5.eps
effects=readdlm("../Result/agren_effects_main_inter_comparison.txt")
pe=layers(labels[:,2],labels[:,end],effects[:,[1,4]])
plot1d(pe;title="Effect plots (Kc=I)",ylabel="Effects size" ,Legend=["Main", "Interaction"],loc="upper left")

#############
## 1d plot
lod1=readdlm("../Result/agren_1dlod.txt")
# # Fig_S4.eps
# p0=layers(labels[:,2],labels[:,end],lod1)
# plot1d(p0;title= "LOD scores by different Kc's ",Legend=["Kc=soil daily range temperature","Kc=air daily range temperature","Kc=weekly drought index", "Kc=Avg(soil,air,drought)","Kc=I"])

# p1=layers(labels[:,2],labels[:,end],lod1[:,[1,end]])
# plot1d(p1;title= "LOD scores with 95% cutoffs",yint=[4.05876 4.17416],yint_color=["red","blue"],Legend=["Kc=soil daily range temperature with blue cuttoff ","Kc=I with red cutoff"])

# Fig_S6.eps
lod11=readdlm("../Result/agren_sitewise_multi_1dlod_it_soil&i_sw.txt")
p11=layers(labels[:,2],labels[:,end],lod11[:,[2,4]])
plot1d(p11;title="LODs by sitewise multivariate genome scan (Kc=I)",Legend=["Italy for 3-year fitness","Sweden for 3-year fitness"])

# Fig_S7.eps (Italy), Fig_S8.eps (Sweden)
lod=readdlm("../Result/agren_univariate_1dlod_it_sw.txt")
pt=layers(labels[:,2],labels[:,end],lod[:,1:3]);pw=layers(labels[:,2],labels[:,end],lod[:,4:end])
plot1d(pt;title="LODs by sitewise univariate genome scan (Italy)",Legend=["Italy 2009","Italy 2010","Italy 2011"])
plot1d(pw;title="LODs by sitewise univariate genome scan (Sweden)",Legend=["Sweden 2009","Sweden 2010","Sweden 2011"],loc="upper left")


# pi=layers(labels[:,2],labels[:,end],hcat(lod[:,1:3],lod11[:,2]));ps=layers(labels[:,2],labels[:,end],hcat(lod[:,4:end],lod11[:,4]))
# plot1d(pi;title="Italy : LODs by sitewise-multivariate & univariate genome scan",Legend=["Italy 2009","Italy 2010","Italy 2011","Italy (3 years) by Kc=I"])
# plot1d(ps;title="Sweden : LODs by sitewise-multivariate & univariate genome scan",Legend=["Sweden 2009","Sweden 2010","Sweden 2011","Sweden (3 years) by Kc=I"],loc="upper left")


# ## 2d plot : Fig_S5
# lod2=readdlm("../Result/agren_2dlod_soil.txt")
# p2=layers(labels1[:,2],labels1[:,end],lod2)
# plot2d(p2)

