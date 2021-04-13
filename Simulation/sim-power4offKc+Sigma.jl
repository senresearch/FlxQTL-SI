using Distributed, StatsBase, Statistics, Distributions,DelimitedFiles, LinearAlgebra
@everywhere using FlxQTL
@everywhere include(simulation.jl)
FlxQTL.setSeed(10,150)

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

#K
kinX=FlxQTL.Markers(chrname,chr,pos,impgen');
@time K=FlxQTL.shrinkgLoco(FlxQTL.kinshipMan,100,kinX);
@time T, Λ=FlxQTL.K2eig(K,true);

#soil daily range temperature (Kc)
soilDailyRng=readdlm("../processedData/agren_soilDailyRng_it_sw_09-12.txt",'\t');
#standardization
sRngStd=(soilDailyRng.-mean(soilDailyRng))./std(soilDailyRng);
@time Ks=FlxQTL.kinshipMan(sRngStd);
@time Ts,λs=FlxQTL.K2eig(Ks);

Z=hcat(ones(6),vcat(-ones(3),ones(3)))
m,q=size(Z)

itr=1000; cr=1;X_sml=25 ;X_lrg=1;w=30
τ2true=2.0.^-[10:-1:0.1;]

idxL=[336]
Bt=[sqrt(2) -sqrt(2);sqrt(2) sqrt(2)]


#### Σ 
# Σ:  default case
a=0.4; b=0.07
A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);B=b*ones(3,3);
Σtrue=[ A B;B A]

#generating data sets
 for i=1:length(τ2true)
    @time Ys0, Ys1 = Ygenerator(itr,τ2true[i],Σtrue,Bt,cr,Z,idxL,X_sml,w,sRngStd,XX)
    open("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_null_traits_Σa$(a).txt","w")do io
        writedlm(io,Ys0)
    end
    open("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_traits_Σa$(a).txt","w")do io
        writedlm(io,Ys1)
    end
end

#genomscan with Kc=I
f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bt).txt","w")
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_null_traits_Σa$(a).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_traits_Σa$(a).txt")
#simulation 
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Matrix(1.0I,m,m),Λ,ones(m),XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
        f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bt).txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)
    end
 

#genomescan with Kc
f=open("/lustre/haven/user/hkim89/results/agren_pwr_B$(Bt).txt","w")
    for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_null_traits_Σa$(a).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_traits_Σa$(a).txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Ts,Λ,λs,XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
       f= open("/lustre/haven/user/hkim89/results/agren_pwr_B$(Bt).txt","a")
        writedlm(f, [τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)
    end


# Σ:  
a=0.3; b=0.07
A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);B=b*ones(3,3);
Σtrue=[ A B;B A]

#generating data sets
 for i=1:length(τ2true)
    @time Ys0, Ys1 = Ygenerator(itr,τ2true[i],Σtrue,Bt,cr,Z,idxL,X_sml,w,sRngStd,XX)
    open("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_null_traits_Σa$(a).txt","w")do io
        writedlm(io,Ys0)
    end
    open("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_traits_Σa$(a).txt","w")do io
        writedlm(io,Ys1)
    end
end

#genomscan with Kc=I
#Pwer=Array{Float64,2}(undef,length(τ2true),2); Mtime=copy(Pwer)
f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bt)_Σa$(a).txt","w")
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_null_traits_Σa$(a).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_traits_Σa$(a).txt")
#simulation 
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Matrix(1.0I,m,m),Λ,ones(m),XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
        f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bt)_Σa$(a).txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)
      
#     Pwer[i,:]=[τ2true[i] pwr]; Mtime[i,:]=medtime
    end
 

#genomescan with Kc
f=open("/lustre/haven/user/hkim89/results/agren_pwr_B$(Bt)_Σa$(a).txt","w")
    for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_null_traits_Σa$(a).txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_traits_Σa$(a).txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Ts,Λ,λs,XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
       f= open("/lustre/haven/user/hkim89/results/agren_pwr_B$(Bt)_Σa$(a).txt","a")
        writedlm(f, [τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)
    end



#### another case
a=0.3; b=0.0; c=0.1
A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);C=c*ones(3,3)+(1-c)*Matrix(1.0I,3,3);B=zeros(3,3)
Σtrue=[A B;B C]

#generating data sets
 for i=1:length(τ2true)
    @time Ys0, Ys1 = Ygenerator(itr,τ2true[i],Σtrue,Bt,cr,Z,idxL,X_sml,w,sRngStd,XX)
    open("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_Σ_a$(a)bc$(c)_null.txt","w")do io
        writedlm(io,Ys0)
    end
    open("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_Σ_a$(a)bc$(c)_traits.txt","w")do io
        writedlm(io,Ys1)
    end
end

#genomscan with Kc=I
f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bt)_Σ_a$(a)bc$(c).txt","w")
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_Σ_a$(a)bc$(c)_null.txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_Σ_a$(a)bc$(c)_traits.txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Matrix(1.0I,m,m),Λ,ones(m),XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
        f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bt)_Σ_a$(a)bc$(c).txt","a")
        writedlm(f,[τ2true[i] pwr Lod1[idxL]])
       close(f)
      
    end
 


#genomescan with Kc
f=open("/lustre/haven/user/hkim89/results/agren_pwr_B$(Bt)_Σ_a$(a)bc$(c).txt","w")
    for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_Σ_a$(a)bc$(c)_null.txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_Σ_a$(a)bc$(c)_traits.txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Ts,Λ,λs,XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
       f= open("/lustre/haven/user/hkim89/results/agren_pwr_B$(Bt))_Σ_a$(a)bc$(c).txt","a")
        writedlm(f, [τ2true[i] pwr Lod1[idxL]])
       close(f)

    end



#generating data sets
 for i=1:length(τ2true)
    @time Ys0, Ys1 = Ygenerator(itr,τ2true[i],Σtrue,Bt,cr,Z,idxL,X_sml,w,sRngStd,XX)
    open("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_Σ=I_null_traits.txt","w")do io
        writedlm(io,Ys0)
    end
    open("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_Σ=I_traits.txt","w")do io
        writedlm(io,Ys1)
    end
end

#sk017
#genomscan with Kc=I ,Σ=I
f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bt)_Σ=I.txt","w") 
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_Σ=I_null_traits.txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_Σ=I_traits.txt")
#simulation 
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Matrix(1.0I,m,m),Λ,ones(m),XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
        f=open("/lustre/haven/user/hkim89/results/agrenI_pwr_B$(Bt)_Σ=I.txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)
      
    end


##### extreme cases for Kc (analyses using the same data just above ): 
#results in ACF were deleted. could not generate result plots.

Kh=H(0.19,m)
@time Th,λh=EcmNestrvQTL.K2eig(Kh);

Kc=Kc_ar(0.9,m)
@time Tc,λc=EcmNestrvQTL.K2eig(Kc);
#shuffle
c=[2;3;6;5;1;4;]
Kp=Kc[c,c]
@time Tp,λp=EcmNestrvQTL.K2eig(Kp);
Σtrue=Matrix(1.0I,m,m)

#sk020
#genomscan with Kc=H(0.19) ,Σ=I

f=open("/lustre/haven/user/hkim89/results/agren_H(0.19)_pwr_B$(Bt)_Σ=I.txt","w") 
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_Σ=I_null_traits.txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_Σ=I_traits.txt")
#simulation 
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Th,Λ,λh,XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
        f=open("/lustre/haven/user/hkim89/results/agren_H(0.19)_pwr_B$(Bt)_Σ=I.txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)

    end

#sk021:Kc=ρ, Σ=I
#genomscan with Kc=ρ ,Σ=I
f=open("/lustre/haven/user/hkim89/results/agren_ρ$(ρ)_pwr_B$(Bt)_Σ=I.txt","w") #ρ=0.9 
for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_Σ=I_null_traits.txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_Σ=I_traits.txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Tc,Λ,λc,XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
        f=open("/lustre/haven/user/hkim89/results/agren_ρ$(ρ)_pwr_B$(Bt)_Σ=I.txt","a")
        writedlm(f,[τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)
      
    end
 

#sk023:Kc=permuted ρ, Σ=I
#genomescan with Kc=ρ permuted ,Σ=I

f=open("/lustre/haven/user/hkim89/results/agren_ρ$(ρ)permute_pwr_B$(Bt)_Σ=I.txt","w")
    for i=1:length(τ2true)
Ys0=readdlm("/lustre/haven/user/hkim89/dataSets/agren_τ2true$(τ2true[i])_Σ=I_null_traits.txt")
Ys1=readdlm("/lustre/haven/user/hkim89/dataSets/agren_B$(Bt)_τ2true$(τ2true[i])_Σ=I_traits.txt")
#simulation soil daily range
@time Lod0,Lod1,B,sdr0,medtime=simulation(itr,cr,idxL,Ys0,Ys1,T,Tp,Λ,λp,XX,Z);
    α,pwr=tErrors(itr,Lod0,Lod1,idxL;erate=[95])  
       f= open("/lustre/haven/user/hkim89/results/agren_ρ$(ρ)permute_pwr_B$(Bt)_Σ=I.txt","a")
        writedlm(f, [τ2true[i] pwr α median(Lod1[idxL,:])])
       close(f)
    end



