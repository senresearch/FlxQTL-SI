#  implemented in linux machine & saved in linux & imac
using Distributed, Statistics, StatsBase, DelimitedFiles, LinearAlgebra
#gererating workers to match the size of threads used by GEMMA
addprocs(32)
@everywhere using FlxQTL
@time setSeed(123)

# data:  unzipped and saved to a csv format
@time anno=readdlm("../testdata/mouse_hs1940_anno.csv",',');
@time gen=readdlm("../testdata/mouse_hs1940_geno.csv",',');
##using imputed phenotypes (1940x3)
@time phen=readdlm("../testdata/mouse_hs1940_imp_pheno.txt");

#FlxQTL
@time rmid=getGenoidx(gen[:,4:end] ,0.01); #10783
X19=Markers(anno[rmid,1],anno[rmid,3],anno[rmid,4],gen[rmid,4:end])

Y=convert(Array{Float64,2},phen');
# m,n=size(Y)

# computing kinship (no loco) and save it to a file (using actual data from GEMMA)
time gemma -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940.pheno.txt -a ../testdata/mouse_hs1940.anno.txt -gk -o kinship -no-check

K=readdlm("../runtest/output/kinship.cXX.txt");
@time T,λ =K2eig(K+0.00001I); #2.15s

#####whole genescan for comparison w/o loco
@time lnp,B,est0= gene1Scan(1,T,λ,Y,X19;LogP=true);
open("../runtest/HS1940mouse3traits-allscan_logP.txt","w")do io
    writedlm(io, lnp)
end

time gemma -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_imp_pheno.txt -n 1 2 3 -k ./output/kinship.cXX.txt -lmm 2 -o lrt3_ch19-v3 -no-check

######  computation time comparison by increasing trait data size


#generating simulated phenotype data

#generating a set of 3-30 traits increased by 3
function sim_Data(rawData, itr)
    n,m=size(rawData);y=[]
    for i=1:itr #how many datasets?
     y=rawData[shuffle(1:n),:]
     
      for j=2:10
         indx=shuffle(1:n)
         y=[y rawData[indx,:]]
          if (rank(y)!= (m+3*(j-1)))
             y[:,(end-2):end]= rawData[shuffle(1:n),:]
             print(rank(y)==size(y,2))
          end
        end
      writedlm("../testdata/ms_HS1940_simulated$(i).txt",y)
    end
end

###### computation times
# measure runtime per dataset.
function runTime(Y,nset,itr) 

     rt=[];mtime=zeros(nset)
    for j=1:nset

        for l=1:itr+1
         rt0 = @elapsed lod,B,est0 = gene1Scan(1,T,λ,Y[1:3j,:],X19);
         rt=[rt;rt0]
        end
        mtime[j]=mean(rt[2:end]) #discard the first
    end
    return mtime
end

#generating 5 simulated datsets per trait set
sim_Data(phen,5)
#check if all data are correctly generated:
for j=1:5
    y=readdlm("../testdata/ms_HS1940_simulated$(j).txt")
    print([size(y) rank(y)])
end

    Y0=readdlm("../testdata/ms_HS1940_simulated1.txt");
    Y0=convert(Array{Float64,2},Y0');
    meantimes = runTime(Y0,10,5)
    writedlm("../Result/time_comparision_flxqtl.txt",[meantimes'])


j=2:5
    Y0=readdlm("../testdata/ms_HS1940_simulated$(j).txt")
    Y0=convert(Array{Float64,2},Y0');
    meantimes = runTime(Y0,10,5)
    f=open("../Result/time_comparision_flxqtl.txt","a")
     writedlm(f,[meantimes'])
     close(f)
end
#  mtim0= runTime(Y);
#  median(mtim0,dims=2)
#  open("mouse_hs1940_runtime4chr_3traits.txt","w")do io
#         writedlm(io,mtim0)
#         end


#     med3tr= runTime(Y)
#     writedlm("../Result/mouse_hs1940_flxqtl_64runtimes_$(m)traits_chr1-3-7-11-15-19.txt",med3tr)
#     display([m mean(med3tr,dims=2)'])
#     writedlm("../Result/mouse_hs1940_flxqtl-mean_time_3-12traitsBychr4incr.txt",[m mean(med3tr,dims=2)'])
      
    
# # for j=1:3
# #     Y=readdlm("../testdata/mouse_hs1940_simulated$(3j+3)traits.txt")
#         Y=convert(Array{Float64,2},Y');  m=size(Y,1)
#         mtime=runTime(Y)
#         writedlm("../Result/mouse_hs1940_flxqtl_64runtimes_$(m)traits_chr1-3-7-11-15-19.txt",mtime)
#     f= open("../Result/mouse_hs1940_flxqtl-mean_time_3-12traitsBychr4incr.txt","a")
#         writedlm(f,[m mean(mtime,dims=2)'])
#         close(f)
#       display([m mean(mtime,dims=2)'])
# end

#6traits: 304s
#9;683s
y=readdlm("../testdata/mouse_hs1940_simulated6traits.txt");y=Float64.(y');
@time lod6,b6,est06=gene1Scan(1,T,λ,y,X19);
@time lod6,b6,est06=gene1Scan(1,T,λ,y,X19;penalize=true);
y9=readdlm("../testdata/mouse_hs1940_simulated9traits.txt");y9=Float64.(y9');
@time lod9,b9,est09=gene1Scan(1,T,λ,y9,X19);
@time lod9,b9,est09=gene1Scan(1,T,λ,y9,X19;penalize=true);

#12; 1292s
y12=readdlm("../testdata/mouse_hs1940_simulated12traits.txt");y12=Float64.(y12');
@time lod12,b12,est12=gene1Scan(1,T,λ,y12,X19);
@time lod12,b12,est12=gene1Scan(1,T,λ,y12,X19;penalize=true);

 y15=readdlm("../testdata/mouse_hs1940_simulated15traits.txt"); y15=Float64.(y15');
 @time lod15,B,est15 =gene1Scan(1,T,λ,y15,X19); #2280 vs 1129s gemma
 @time lod15,B,est15 =gene1Scan(1,T,λ,y15,X19;penalize=true,df_prior=Int64(ceil(1.5*size(y15,1)))); #3778s
 
  y18=readdlm("../testdata/mouse_hs1940_simulated18traits.txt");y18=Float64.(y18')
@time lod18,B,est18 =gene1Scan(1,T,λ,y18,X19); # 3126,3153s vs 6404s gemma
@time lod18,B,est18 =gene1Scan(1,T,λ,y18,X19;penalize=true,df_prior=Int64(ceil(1.5*size(y18,1)))); #6529s
 y21=readdlm("../testdata/mouse_hs1940_simulated21traits.txt");y21=Float64.(y21');rank(y21)
@time lod21,B,est0 =gene1Scan(1,T,λ,y21,X19); #4595,4405s vs 5710s gemma
sum(lod21.<0)
# @time lod21,B,est21 =gene1Scan(1,T,λ,y21,X19;penalize=true,df_prior=Int64(ceil(1.9*size(y21,1)))); 


y24=readdlm("../testdata/mouse_hs1940_simulated24traits.txt");y24=Float64.(y24');rank(y24) # m=24
@time lod24,B,est24 =gene1Scan(1,T,λ,y24,X19); sum(lod24.<0) #6203 vs 48,396s gemma
@time lod,B,est24 =gene1Scan(1,T,λ,y2,X19;penalize=true,df_prior=Int64(ceil(1.9*size(y24,1)))); #
y27=readdlm("../testdata/mouse_hs1940_simulated27traits.txt");y27=Float64.(y27') #m=27
@time lod27,B,est27 =gene1Scan(1,T,λ,y27,X19); #8239 vs 128,784s gemma
@time lod27,B,est27 =gene1Scan(1,T,λ,y27,X19;penalize=true,df_prior=Int64(ceil(1.9*size(y27,1)))); #
y30=readdlm("../testdata/mouse_hs1940_simulated30traits.txt");y30=Float64.(y30');rank(y30) #m=30
@time lod30,B,est30 =gene1Scan(1,T,λ,y30,X19); sum(lod30.<0) #10519s vs 35624s gemma
@time lod30,B,est30 =gene1Scan(1,T,λ,y30,X19;penalize=true,df_prior=Int64(ceil(1.9*size(y30,1)))); #

#gemma test
(time gemma -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_simulated15traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 -k ./output/kinship.cXX.txt -lmm 2 -o lrt15_ch19-v3 -no-check)2>> gemma_runtime15_30trait.txt
(time gemma -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_simulated18traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 -k ./output/kinship.cXX.txt -lmm 2 -o lrt18_ch19-v3 -no-check)2>> gemma_runtime15_30trait.txt
(time gemma -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_simulated21traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 -k ./output/kinship.cXX.txt -lmm 2 -o lrt21_ch19-v3 -no-check)2>> gemma_runtime21_30trait1.txt
(time gemma -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_simulated24traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 -k ./output/kinship.cXX.txt -lmm 2 -o lrt24_ch19-v3 -no-check)2>> gemma_runtime15_30trait.txt
(time gemma -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_simulated27traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 -k ./output/kinship.cXX.txt -lmm 2 -o lrt27_ch19-v3 -no-check)2>> gemma_runtime21_30trait.txt
(time gemma -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_simulated30traits.txt -n 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 -k ./output/kinship.cXX.txt -lmm 2 -o lrt30_ch19-v3 -no-check)2>> gemma_runtime21_30trait1.txt


# #Chr1 only
# function fqtl_time(Y)
#     fqtl=zeros(64);
#     for j=1:64
#     fqtl[j]= @elapsed lod,B,est0 =geneScan(1,T,λ,Y,X1);
#     end
#     return fqtl
# end


#    f= open("mouse_hs1940_chr1runtime_6_9_12traits.txt","w")
# for j=1:3
#     Y=readdlm("../testdata/mouse_hs1940_simulated$(3j+3)traits.txt")
#         Y=convert(Array{Float64,2},Y');  m=size(Y,1)
#         mtime=fqtl_time(Y)
#         display([m median(mtime)])
#        f= open("mouse_hs1940_chr1runtime_6_9_12traits.txt","a")
#         writedlm(f,mtime')
#         close(f)
# end



# obtained by running gemma-time.jl (averages)
gem3=[5.505;9.162;15.882;22.883;28.806;39.414]
gem6=[46.096; 52.095; 56.170; 61.541; #1m1.541
65.014 ;#1m5.014
71.697;#1m11.697s
76.521;	#1m16.521s
81.597;	#1m21.597s
88.256;#1m28.256
91.733;  #1m31.733s
97.610; #1m37.610s
101.556;#1m41.556s
107.260;#1m47.260s
113.233; #1m53.233s
115.295; #1m55.295s
120.114;#2m0.114s
155.079; #2m35.079s
156.540; #2m36.540s
158.784; #2m38.784s
]
gem9=[119.616;#1m59.616s
132.672;#2m12.672s
142.125;#2m22.125s
165.681;#2m45.681s
175.495;#2m55.495s
188.866;#3m8.866s
199.973;#3m19.973s
213.962;#3m33.962s
223.197; #3m43.197s
228.511;  #3m48.511s
238.906;#3m58.906s
246.001;#4m6.001s
254.428;#4m14.428s
269.421;#4m49.421s
293.310;#4m53.310s
299.213;#4m59.213s
313.118;#5m13.118s
316.418;#5m16.418s √
320.700 #5m20.700s  √
 ]
#gemma 12 traits 1-5 accumulated chr's
gem12=  [
 489.096;#8m9.096s
 531.305;#8m51.305s
 558.899;#9m18.899s
 601.354;#10m1.354s
 625.388;#10m25.388s
 681.019; #11m21.019s
 854.534;#14m14.534s
 895.733; #14m55.733s
 926.200;#15m26.200s
 943.350;#15m43.350s
 986.028;#16m26.028s
 1007.041;#16m47.041s
 1030.773;#17m10.773s
 1050.310; #17m30.310s  √
 1089.820; #18m9.820s   √
 1104.356;#18m24.356s
1134.021; #18m54.021s √
1154.673;#19m14.673s √
1186.045;#19m46.045s √
]

#combine together

Tdata=[[3;gem3] [6;gem6[idx]] [9;gem9[idx]] [12;gem12[idx]]]
open("../test/mouse_hs1940_gemma-mean_time_chr1-3-7-11-15-19.txt","w")do io
    writedlm(io,Tdata)
end

--------
#### Fig4.eps
using Plots
pyplot()
using Plots.PlotMeasures


using FlxQTL
#scanned value comparison using -log10(P-value)
# flx=readdlm("../Result/HS1940mouse3traits-allscan_lnP.txt")
# flx= lod2logP(flx0[:,1] ,3)
flx=readdlm("./result/hs1940_w_no_locoByflxqtl_logP.txt")
gem0=readdlm("../Result/lrt3_ch19-v3.assoc.txt",'\t';header=true)

#qqplot
exptd=collect(1:size(flx,1))
lexp=-log10.(exptd/(size(flx,1)+1))
fobs=sort(flx[:,end];rev=true)
floco=sort(flx[:,1];rev=true);
gobs=-log10.(sort(gem0[1][:,end]))
# x=collect(1:8)
x=collect(1:11)


p1=scatter(lexp,fobs,label="FlxQTL",xlims=(0,11),ylims=(0,11),xlab="Expected -log10(P)",ylab="Observed -log10(P)",aspectratio=1,grid=false,guidefontsize=25,tickfontsize=20,legendfontsize=20, marker=8,
    bottom_margin=0mm,right_margin=0mm)
 scatter!(lexp,floco,label="FlxQTL_loco",marker=8)
scatter!(lexp,gobs,label="GEMMA",marker=8)
plot!(x,x,linecolor=:black, label="y=x")




#scan comparison

p2=plot(flx,label="FlxQTL",ylab="-log10(P)",xlab="SNPs",guidefontsize=25,tickfontsize=20,legendfontsize=20,line=3)
#right_margin=2.5mm, top_margin=0mm,left_margin=-2.0mm,line=3)
plot!(-log10.(gem0[1][:,end]),label="GEMMA",grid=false,line=3)
hline!([-log10(0.05/10783)],line=(3,:dash,:grey),label=false)
annotate!((10000,-log10(0.05/length(gem0[1][:,end]))+0.2 , Plots.text("α=0.05", 20, :grey, :right)))

#run time comparison
 idx=[1;3; 7; 11; 15; 19];
ptime=readdlm("../Result/runtime-flxqtl_vs_gemma_hsmouse_chr1-3-7-11-15-19.txt")

legendLabel =["FlxQTL 3 traits" "  6 traits" "  9 traits" "  12 traits" "GEMMA 3 traits" "  6 traits" "  9 traits" " 12 traits"]
markerLabel = [:rect :diamond :dot :X :rect :diamond :dot :X]
myLineStlye = [:dash :solid :dot :dashdot :dash :solid :dot :dashdot]
myColor = [RGB(0,155/255,250/255) RGB(227/255,111/255,71/225) ]

p3 = plot(idx, log2.(ptime[2:end,1]), color=myColor[Int(1>4)+1], line= (3, myLineStlye[1]),
    xticks=idx,yticks=collect(3:2:10),
    xlab="Cumulative chromosomes", ylab="log2(Time (sec))",guidefontsize=25,tickfontsize=20, legendfontsize=20,
    grid=false, label = legendLabel[1], left_margin=-2mm,right_margin=2.5mm,legend=:bottomright)
scatter!(idx, log2.(ptime[2:end,1]),color=myColor[Int(1>4)+1],  label="",marker=6)

for i in 2:8
    plot!(idx, log2.(ptime[2:end,i]), color=myColor[Int(i>4)+1], line = (3, myLineStlye[i]),
          xlab="Cumulative chromosomes", ylab="log2(Time (sec))",
          grid=false, label = legendLabel[i])
    scatter!(idx, log2.(ptime[2:end,i]), color=myColor[Int(i>4)+1], label="",marker=6)
end
p3

ll = @layout [ a{1.1h} b;
              c]
plot(p3,p1,p2,layout=ll)
