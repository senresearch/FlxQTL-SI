# julia --machine-file=$PBS_NODEFILE
# addprocs to different nodes : addprocs([("node1",100)])
# to see if julia communicates with different nodes
#for i in workers()
#    host, pid = fetch(@spawnat i (gethostname(), getpid()))
#    println("Hello from process $(pid) on host $(host)!")
#end

using Distributed, StatsBase, Statistics, Distributions, DelimitedFiles, LinearAlgebra

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

#K
kinX=FlxQTL.Markers(chrname,chr,pos,impgen');
@time K=FlxQTL.shrinkgLoco(FlxQTL.kinshipMan,100,kinX);
@time T, Λ=FlxQTL.K2eig(K,true);

#climate covariates: daily soil/air temperature ranges: it 7/2009-6/2012, sw 7/2009-6/2012
# 365 days by 6 envir
soilDailyRng=readdlm("../processedData/agren_soilDailyRng_it_sw_09-12.txt",'\t');

#standardization
sRngStd=(soilDailyRng.-mean(soilDailyRng))./std(soilDailyRng);
@time Ks=FlxQTL.kinshipMan(sRngStd);
@time Ts,λs=FlxQTL.K2eig(Ks);


#sampling weather indices
#q1=365 ; q1=52
w=30
p,n=size(XX.X)


itr=1000; cr=1;X_sml=25 ;X_lrg=1;
τ2true=2.0.^-[10:-1:0.1;]

idxL=[336]

# Σ:  (a,b)=(0.4,0.07): default
a=0.4; b=0.07
A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);B=b*ones(3,3);
Σtrue=[ A B;B A]



### power varying with B for each τ2 
Btrue= [0.25  -0.25 ;0.25   0.25]
for l=1:10
B1=Btrue*sqrt(2)^(l-1) ; 
    f= open("../results/agren_pwr_$(l)th_B$(B1).txt","w")  
    for i=1:length(τ2true)
         # generate datasets (H0, H1)
    @time Ys0, Ys1 = Ygenerator(itr,τ2true[i],Σtrue,B1,cr,Z,idxL,X_sml,w,sRngStd,XX)
     #simulation soil daily range
    @time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Ts,Λ,λs,XX,Z);
       α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])
          f=open("../results/agren_pwr_$(l)th_B$(B1).txt","a")
             writedlm(f, [τ2true[i] pwr])
           close(f)

    end
end


## power varying with B that narrows down
δ=sqrt(2).^[3.:0.5:5.;]
B1=zeros(2,2,length(δ));
for j=1:length(δ)
    B1[:,:,j]=round.(Btrue*δ[j],digits=5)
end
    
#generating data sets
for j=1:length(δ)
 for i=1:length(τ2true)
    @time Ys0, Ys1 = Ygenerator(itr,τ2true[i],Σtrue,B1[:,:,j],cr,Z,idxL,X_sml,w,sRngStd,XX)
    open("/lustre/haven/user/hkim89/dataSets/agren_H0_τ2true$(τ2true[i]).txt","w")do io
        writedlm(io,Ys0)
    end
    open("/lustre/haven/user/hkim89/dataSets/agren_H1_τ2true$(τ2true[i])_B$(B1[:,:,j]).txt","w")do io
        writedlm(io,Ys1)
    end
 end
end


## Kc
for j=1:length(δ)
f=open("/lustre/haven/user/hkim89/results/agren1_pwr_$(j)th_B$(B1[:,:,j]).txt","w")
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H0_τ2true$(τ2true[i]).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H1_τ2true$(τ2true[i])_B$(B1[:,:,j]).txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Ts,Λ,λs,XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])
        f=open("/lustre/haven/user/hkim89/results/agren1_pwr_$(j)th_B$(B1[:,:,j].txt","a")
        writedlm(f,[τ2true[i] pwr ])
       close(f)
    end
end

##checking Kc=I
for j=1:length(δ)
f=open("/lustre/haven/user/hkim89/results/agren2_I_pwr_$(j)th_B$(B1[:,:,j]).txt","w")
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H0_τ2true$(τ2true[i]).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_H1_τ2true$(τ2true[i])_B$(B1[:,:,j]).txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Matrix(1.0I,m,m),Λ,ones(m),XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])
        f=open("/lustre/haven/user/hkim89/results/agren2_I_pwr_$(j)th_B$(B1[:,:,j]).txt","a")
        writedlm(f,[τ2true[i] pwr ])
       close(f)

    end
end




#####










