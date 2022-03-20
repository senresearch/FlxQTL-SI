# julia --machine-file=$PBS_NODEFILE
# addprocs to different nodes : addprocs([("node1",100)])
# to see if julia communicates with different nodes
#for i in workers()
#    host, pid = fetch(@spawnat i (gethostname(), getpid()))
#    println("Hello from process $(pid) on host $(host)!")
#end

using Distributed, StatsBase, Statistics, Distributions,DelimitedFiles, LinearAlgebra
@everywhere using FlxQTL
@everywhere include("simulation.jl")

@time FlxQTL.setSeed(123);
## read imputed genotype data: a(=1):italian parent, b(=2)swedish parent
impgen = readdlm("../processedData/fullrank_imput.csv",',';skipstart=1);
impgen[impgen.==1.0].=0.0;impgen[impgen.==2.0].=1.0;
gen1=copy(impgen);
gen1[gen1.==0.0].=-1.0;

#genotype data for 1D scan
labels=readdlm("../processedData/marlabels_agren.csv",',';skipstart=1);
chrname=convert(Array{String,1},labels[:,1])
chr=convert(Array{Any,1},labels[:,2]);
pos=convert(Array{Float64,1},labels[:,3]);
XX=EcmNestrvQTL.Markers(chrname,chr,pos,gen1');

Z=hcat(ones(6),vcat(-ones(3),ones(3)))
m,q=size(Z)


#climate covariates: daily soil/air temperature ranges: it 7/2009-6/2012, sw 7/2009-6/2012
# 365 days by 6 envir
soilDailyRng=readdlm("../processedData/agren_soilDailyRng_it_sw_09-12.txt",'\t');
#soilDailyAvgs=readdlm("../processedData/agren_soilDailyAvgs_it_sw_09-12.txt",'\t');

#standardization
sRngStd=(soilDailyRng.-mean(soilDailyRng))./std(soilDailyRng);
#sAvgStd=(soilDailyAvgs.-mean(soilDailyAvgs))./std(soilDailyAvgs);

#satellite derived weekly measure of draught :52 wks by 6 envir
#wDrg=readdlm("../processedData/agren_wkly_drought_it_sw_09-12.txt",'\t');
#DrgStd=(wDrg.-mean(wDrg))./std(wDrg);

# compute K
kinX=FlxQTL.Markers(chrname,chr,pos,impgen');
@time K=FlxQTL.shrinkgLoco(FlxQTL.kinshipMan,100,kinX);
@time T, Λ=FlxQTL.K2eig(K,true);

#sampling weather indices
#q1=365 ; q1=52
w=30
p,n=size(XX.X)
@time Ks=FlxQTL.kinshipMan(sRngStd);
@time Ts,λs=FlxQTL.K2eig(Ks);

itr=1000; cr=1;X_sml=25 ;X_lrg=1;
τ2true=2.0.^-[10:-1:0.1;]

Bfix=zeros(2,2,4);
Bfix[:,:,1]=[-1.0 -1.0;1.0 1.0]
for j=2:4
    Bfix[:,:,j]=Bfix[:,:,1].*sqrt(2)^(j-1)
end

# Σ:  (a,b)=(0.4,0.07),(0.3,0.07)
a=0.4; b=0.07
A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);B=b*ones(3,3);
Σtrue=[ A B;B A]



idxL=[336]

#generating data sets
for j=1:4
 for i=1:length(τ2true)
    @time Ys0, Ys1 = Ygenerator(itr,τ2true[i],Σtrue,Bfix[:,:,j],cr,Z,idxL,X_sml,w,sRngStd,XX)
    open("/lustre/haven/user/hkim89/dataSets/agren_H0_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt","w")do io
        writedlm(io,Ys0)
    end
    open("/lustre/haven/user/hkim89/dataSets/agren_H1_B$(Bfix[:,:,j])_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt","w")do io
        writedlm(io,Ys1)
    end
 end
end





############checking Z's
## Z : site contrasts
for j=1:4
f=open("/lustre/haven/user/hkim89/results/agrenZ_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt","w")

  for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H0_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H1_B$(Bfix[:,:,j])_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Ts,Λ,λs,XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])
        f=open("/lustre/haven/user/hkim89/results/agrenZ_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)

#     Pwer[i,:]=[τ2true[i] pwr]; Mtime[i,:]=medtime
    end
end


# Z=I
for j=1:4
f=open("/lustre/haven/user/hkim89/results/agrenZi_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt","w")

  for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H0_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H1_B$(Bfix[:,:,j])_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
#simulation soil daily range 
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Ts,Λ,λs,XX,Matrix(1.0I,m,m));
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])
        f=open("/lustre/haven/user/hkim89/results/agrenZi_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)

#     Pwer[i,:]=[τ2true[i] pwr]; Mtime[i,:]=medtime
    end
end



## MVLMM case
for j=1:4
f=open(realpath(string(@__DIR__,"/../results"))*"/agren_mvlmm_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt","w")

for i=1:length(τ2true)
Ys0=readdlm(realpath(string(@__DIR__,"/../dataSets"))*"/agren_H0_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
Ys1=readdlm(realpath(string(@__DIR__,"/../dataSets"))*"/agren_H1_B$(Bfix[:,:,j])_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
#simulation 
@time Lod0,Lod1, B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Λ,XX,Matrix(1.0I,m,m))
     α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])
        f=open(realpath(string(@__DIR__,"/../results"))*"/agren_mvlmm_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)

#     Pwer[i,:]=[τ2true[i] pwr]; Mtime[i,:]=medtime
    end
end


######################




##checking Kc=I
f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bfix[:,:,1])_Σa$(a)_b$(b).txt","w")
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H0_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H1_B$(Bfix[:,:,1])_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Matrix(1.0I,m,m),Λ,ones(m),XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])
        f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bfix[:,:,1])_Σa$(a)_b$(b).txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)

#     Pwer[i,:]=[τ2true[i] pwr]; Mtime[i,:]=medtime
    end

#Kc !=I
f=open("/lustre/haven/user/hkim89/results/agren_pwr_B$(Bfix[:,:,1])_Σa$(a)_b$(b).txt","w")
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H0_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H1_B$(Bfix[:,:,1])_τ2true$(τ2true[i])_Σa$(a)_b$(b).txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Ts,Λ,λs,XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])
        f=open("/lustre/haven/user/hkim89/results/agren_pwr_B$(Bfix[:,:,1])_Σa$(a)_b$(b).txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)

#     Pwer[i,:]=[τ2true[i] pwr]; Mtime[i,:]=medtime
    end




