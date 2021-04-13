using DelimitedFiles, StatsBase, Statistics, LinearAlgebra, Distributed
using RCall
using Gadfly

@everywhere using FlxQTL
FlxQTL.setSeed(10,200)

## read imputed genotype data: a(=1):italian parent, b(=2)swedish parent
impgen = readdlm("processedData/fullrank_imput.csv",',';skipstart=1);
impgen[impgen.==1.0].=0.0;impgen[impgen.==2.0].=1.0;
impgen1=copy(impgen);
impgen1[impgen1.==0.0].=-1.0;

#genotype data for 1D scan
labels=readdlm(string(@__DIR__,"/../processedData/marlabels_agren.csv"),',';skipstart=1);
XX=FlxQTL.Markers(labels[:,1],labels[:,2],labels[:,3],impgen1');

## read phenotype data
pheno =readdlm(string(@__DIR__,"/../processedData/pheno2013_imp.csv"),',';header=true);
pheno=pheno[1][:,2:end-1];
Y=convert(Array{Float64,2},pheno');
Ystd=(Y.-mean(Y,dims=2))./std(Y,dims=2);

#site contrasts: italy (-1), sweden (1)
Z=hcat(ones(6),vcat(-ones(3),ones(3)));
m,q=size(Z)


#climate covariates: daily soil/air temperature ranges: it 7/2009-6/2012, sw 7/2009-6/2012
#soil daily range temperature
# 365 days by 6 envir
soilDailyRng=readdlm(string(@__DIR__,"/../processedData/agren_soilDailyRng_it_sw_09-12.txt"),'\t');
sRngStd=(soilDailyRng.-mean(soilDailyRng))./std(soilDailyRng);

#air daily range temperature
airDailyRng=readdlm("../processedData/agren_airDailyRng_it_sw_09-12.txt",'\t');
aRngStd=(airDailyRng.-mean(airDailyRng))./std(airDailyRng);

#satellite derived weekly measure of draught :52 wks by 6 envir
wDrg=readdlm("../processedData/agren_wkly_drought_it_sw_09-12.txt",'\t' )
DrgStd=(wDrg.-mean(wDrg))./std(wDrg);


#Kinship
Xkin=FlxQTL.Markers(labels[:,1],labels[:,2],labels[:,3],impgen');
@time K=FlxQTL.shrinkgLoco(FlxQTL.kinshipMan,100,Xkin);
@time T, Λ=FlxQTL.K2eig(K,true);

@time Ks=FlxQTL.kinshipMan(sRngStd);
@time Ka=FlxQTL.kinshipMan(aRngStd);
@time Kd=FlxQTL.kinshipMan(DrgStd);
 Kc= @. 1/3*(Ks+Ka+Kd)

@time Ts,λs=FlxQTL.K2eig(Ks)
@time Ta,λa=FlxQTL.K2eig(Ka)
@time Td,λd=FlxQTL.K2eig(Kd)
@time Tc,λc=FlxQTL.K2eig(Kc)

# ## OR using one function
# @time T,Λ,Ts,λs = FlxQTL.K2Eig(K,Ks,true)

#multivariate genome scans
@time LODs,Bs,est0s=FlxQTL.geneScan(1,T,Ts,Λ,λs,Ystd,XX,Z,true);
@time LODa,Ba,est0a=FlxQTL.geneScan(1,T,Ta,Λ,λa,Ystd,XX,Z,true);
@time LODd,Bd,est0d=FlxQTL.geneScan(1,T,Td,Λ,λd,Ystd,XX,Z,true);
@time LODc,Bc,est0c=FlxQTL.geneScan(1,T,Tc,Λ,λc,Ystd,XX,Z,true);
@time LODi,Bi,est0i=FlxQTL.geneScan(1,T,Matrix(1.0I,m,m),Λ,ones(m),Ystd,XX,Z,true);
@time Lodg, Bg,est0g=FlxQTL.geneScan(1,T,Λ,Ystd,XX,Z,true); #MLMM

open("../result/agren_1dlod_soil_air_drgt_avg_i.txt","w")do io
    writedlm(io,[LODs LODa LODd LODc LODi])
end

### rearrange effects data 
main=zeros(699,3); inter=copy(main); 
for l=1:699
    main[l,:]=[Bi[1,2,l] Bc[1,2,l] Bg[1,2,l]]
    inter[l,:]=[Bi[2,2,l] Bc[2,2,l] Bg[2,2,l]]
end
#store effects in a file
open("../result/agren_effects_main_inter_comparison.txt", "w") do io
           writedlm(io, [main inter])
       end


#italy
@time Kt=FlxQTL.kinshipMan(sRngStd[:,1:3]);
@time Tt,λt=FlxQTL.K2eig(Kt)
@time LODt,bt,est0t=FlxQTL.geneScan(1,T,Tt,Λ,λt,Ystd[1:3,:],XX,Matrix(1.0I,3,3),true);
@time LODti,bti,est0it=FlxQTL.geneScan(1,T,Matrix(1.0I,3,3),Λ,ones(3),Ystd[1:3,:],XX,Matrix(1.0I,3,3),true);

#Sweden
@time Kw=FlxQTL.kinshipMan(sRngStd[:,4:end]);
@time Tw,λw=FlxQTL.K2eig(Kw)
@time LODw,bw,est0w=FlxQTL.geneScan(1,T,Tw,Λ,λw,Ystd[4:end,:],XX,Matrix(1.0I,3,3),true);
@time LODwi,bwi,est0wi=FlxQTL.geneScan(1,T,Matrix(1.0I,3,3),Λ,ones(3),Ystd[4:end,:],XX,Matrix(1.0I,3,3),true);

open("../result/agren_sitewise_multi_1dlod_it_soil&i_sw.txt","w")do io
    writedlm(io,[LODt LODti LODw LODwi])
end

#run univariate genomescan
 uni_h0=Array{Any,2}(undef,5,m);uni_lods=zeros(699,m);uni_Bs=zeros(1,2,699,6);
for l=1:m
  @time uni_lods[:,l],uni_Bs[:,:,:,l],uni_h0[:,l]=FlxQTL.geneScan(1,T,ones(1,1),Λ,ones(1),Ystd[[l],:],XX,ones(1,1),true)
end

open("../result/agren_univariate_1dlod_it_sw.txt","w")do io
    writedlm(io,uni_lods)
end

##2d-scan
labels1=readdlm("../processedData/agren_newmarkers_multipleqtl.csv",',';skipstart=1);
G1=readdlm("../processedData/agren_gene_imp_revision.csv",',';skipstart=1);
Xp=FlxQTL.Markers(labels1[:,1],labels1[:,2],labels1[:,3],G1');
#@time K2=shrinkgLoco(kinshipMan,100,Xp); ## sec
#T, Λ=EcmNestrvQTL.K2eig(K2,true);

@time LOD2,est02=FlxQTL.gene2Scan(1,T,Ts,Λ, λs,Ystd,Xp,Z,true);
open("../result/agren_2dlod_soil.txt","w")do io
    writedlm(io,LOD2)
end


##permutation
@time K0=FlxQTL.shrinkg(FlxQTL.kinshipMan,100,Xkin.X)
# @time MaxLod,H1par,cutoff=FlxQTL.permTest(1000,1,K0,Ks,Ystd,XX,Z);
# open("../result/agren_soil_maxlod_cutoff=$(cutoff).txt","w")do io
#     writedlm(io,MaxLod)
# end

@time MaxLod,H1par,cutoff=FlxQTL.permTest(1000,1,K0,Matrix(1.0I,m,m),Ystd,XX,Z);
# open("../result/agren_i_maxlod_cutoff=$(cutoff).txt","w")do io
#     writedlm(io,MaxLod)
# end

## Rearrange effects from permutations to get 95% band 
bs=@distributed (vcat) for l=1:1000*699
        [H1par[l].B[1,2] H1par[l].B[2,2]]
 end    

Bmain_per=zeros(699,1000); Binter_per=zeros(699,1000);
for k=1:1000
     Bmain_per[:,k]=bs[699*(k-1)+1:699*(k-1)+699,1]
     Binter_per[:,k]=bs[699*(k-1)+1:699*(k-1)+699,2]
end

# two-sided bands for main and interaction effects, respectively
perm_eff=Array{Float64,2}(undef,699,4);
for i=1:699
    perm_eff[i,1]=percentile(Bmain_per[i,:],2.5)
    perm_eff[i,2]=percentile(Bmain_per[i,:],97.5)
    perm_eff[i,3]=percentile(Binter_per[i,:],2.5)
    perm_eff[i,4]=percentile(Binter_per[i,:],97.5)
end

####### Effects plots
include("effectplot.jl")

pmain_perm=Mlayers(XX.chr,XX.pos,hcat(main[:,1],perm_eff[:,1:2]))
pint_perm=Mlayers(XX.chr,XX.pos,hcat(inter[:,1],perm_eff[:,3:4]))

effectplot(pmain_perm,true;title="Main effect with 95% band", ylabel="effect size")
effectplot(pint_perm,true;title="Interaction effect with 95% band", ylabel="effect size")

peff=Mlayers(XX.chr,XX.pos,hcat(main[:,1],inter[:,1]))
effect2plot(peff;title="Effect Plots",ylabel="Effect size",legend=["Main","Interaction"])

#############
## 1d plot
lod1=readdlm("../result/agren_1dlod_soil_air_drgt_avg_i.txt")

p0=layers(labels[:,2],labels[:,end],lod1)
plot1d(p0;title= "LOD scores by different Kc's ",Legend=["Kc=soil daily range temperature","Kc=air daily range temperature","Kc=weekly drought index", "Kc=Avg(soil,air,drought)","Kc=I"])

# p1=layers(labels[:,2],labels[:,end],lod1[:,[1,end]])
# plot1d(p1;title= "LOD scores with 95% cutoffs",yint=[4.05876 4.17416],yint_color=["red","blue"],Legend=["Kc=soil daily range temperature with blue cuttoff ","Kc=I with red cutoff"])

lod11=readdlm("../result/agren_sitewise_multi_1dlod_it_soil&i_sw.txt")
p11=layers(labels[:,2],labels[:,end],lod11[:,[2,4]])
plot1d(p11;title="LODs by sitewise multivariate genome scan (Kc=I)",Legend=["Italy for 3-year fitness","Sweden for 3-year fitness"])

lod=readdlm("../result/agren_univariate_1dlod_it_sw.txt")
pt=layers(labels[:,2],labels[:,end],lod[:,1:3]);pw=layers(labels[:,2],labels[:,end],lod[:,4:end])
plot1d(pt;title="LODs by sitewise univariate genome scan (Italy)",Legend=["Italy 2009","Italy 2010","Italy 2011"])
plot1d(pw;title="LODs by sitewise univariate genome scan (Sweden)",Legend=["Sweden 2009","Sweden 2010","Sweden 2011"],loc="upper left")

pi=layers(labels[:,2],labels[:,end],hcat(lod[:,1:3],lod11[:,2]));ps=layers(labels[:,2],labels[:,end],hcat(lod[:,4:end],lod11[:,4]))
plot1d(pi;title="Italy : LODs by sitewise-multivariate & univariate genome scan",Legend=["Italy 2009","Italy 2010","Italy 2011","Italy (3 years) by Kc=I"])
plot1d(ps;title="Sweden : LODs by sitewise-multivariate & univariate genome scan",Legend=["Sweden 2009","Sweden 2010","Sweden 2011","Sweden (3 years) by Kc=I"],loc="upper left")


## 2d plot
lod2=readdlm("../result/agren_2dlod_soil.txt")
p2=layers(labels1[:,2],labels1[:,end],lod2)
plot2d(p2)

