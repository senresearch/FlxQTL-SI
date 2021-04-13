#  implemented in linux machine & saved in linux & imac
@everywhere using FlxQTL
@time FlxQTL.setSeed(10,200)

# data:  unzipped and saved to a csv format
@time anno=readdlm("../testdata/mouse_hs1940_anno.csv",',');
@time gen=readdlm("../testdata/mouse_hs1940_geno.csv",',');
##using imputed phenotypes (1940x3)
@time phen=readdlm("../testdata/mouse_hs1940_imp_pheno.txt");

######  computation time comparison by increasing genotype & trait data size
#sort genotype by chr 

ch1=findall(anno[:,3].==1); #1044
ch2=findall(anno[:,3].==2); #948
ch3=findall(anno[:,3].==3); #857
ch4=findall(anno[:,3].==4); #778
ch5=findall(anno[:,3].==5); #770
ch6=findall(anno[:,3].==6); #709
ch7=findall(anno[:,3].==7); #658
ch8=findall(anno[:,3].==8); #615
ch9=findall(anno[:,3].==9); #630
ch10=findall(anno[:,3].==10);#481
ch11=findall(anno[:,3].==11);#706
ch12=findall(anno[:,3].==12);#550
ch13=findall(anno[:,3].==13);#573
ch14=findall(anno[:,3].==14);#590
ch15=findall(anno[:,3].==15);#527
ch16=findall(anno[:,3].==16);#497
ch17=findall(anno[:,3].==17);#535
ch18=findall(anno[:,3].==18);#456
ch19=findall(anno[:,3].==19);#302

#write files for accumulated chr's 
open("../testdata/mouse_hs1940_genochr18.txt","w") do io
    writedlm(io, gen[[ch1;ch2;ch3;ch4;ch5;ch6;ch7;ch8;ch9;ch10;ch11;ch12;ch13;ch14;ch15;ch16;ch17;ch18],:])
end


#simulated phenotype data
n,s=size(phen)
indx=zeros(Int64,n,3)
for j=1:3
 indx[:,j]=shuffle(1:n)
end

Y6=[phen phen[indx[:,1],:]]; rank(Y6)
Y9= [Y6 phen[indx[:,2],:]]; rank(Y9)
Y12=[Y9 phen[indx[:,3],:]]; rank(Y12)    

open("../testdata/mouse_hs1940_simulated6traits.txt","w")do io
        writedlm(io, Y6)
end

open("../testdata/mouse_hs1940_simulated9traits.txt","w")do io
        writedlm(io, Y9)
end
open("../testdata/mouse_hs1940_simulated12traits.txt","w")do io
        writedlm(io, Y12)
end
        



#FlxQTL

gen1=readdlm("../testdata/mouse_hs1940_genochr1.txt")[:,4:end]; @time rmid1=FlxQTL.getGenoidx(gen1 ,0.01); #950
gen2=readdlm("../testdata/mouse_hs1940_genochr2.txt")[:,4:end]; @time rmid2=FlxQTL.getGenoidx(gen2 ,0.01); #1776
gen3=readdlm("../testdata/mouse_hs1940_genochr3.txt")[:,4:end]; @time rmid3=FlxQTL.getGenoidx(gen3 ,0.01); #2570
gen4=readdlm("../testdata/mouse_hs1940_genochr4.txt")[:,4:end]; @time rmid4=FlxQTL.getGenoidx(gen4 ,0.01); #3311
gen5=readdlm("../testdata/mouse_hs1940_genochr5.txt")[:,4:end]; @time rmid5=FlxQTL.getGenoidx(gen5 ,0.01); #3909
gen6=readdlm("../testdata/mouse_hs1940_genochr6.txt")[:,4:end]; @time rmid6=FlxQTL.getGenoidx(gen6 ,0.01); #4564
gen7=readdlm("../testdata/mouse_hs1940_genochr7.txt")[:,4:end]; @time rmid7=FlxQTL.getGenoidx(gen7 ,0.01); #5151
gen8=readdlm("../testdata/mouse_hs1940_genochr8.txt")[:,4:end]; @time rmid8=FlxQTL.getGenoidx(gen8 ,0.01); #5696
gen9=readdlm("../testdata/mouse_hs1940_genochr9.txt")[:,4:end]; @time rmid9=FlxQTL.getGenoidx(gen9 ,0.01); #6259
gen10=readdlm("../testdata/mouse_hs1940_genochr10.txt")[:,4:end]; @time rmid10=FlxQTL.getGenoidx(gen10 ,0.01); #6634
gen11=readdlm("../testdata/mouse_hs1940_genochr11.txt")[:,4:end]; @time rmid11=FlxQTL.getGenoidx(gen11 ,0.01); #7304
gen12=readdlm("../testdata/mouse_hs1940_genochr12.txt")[:,4:end]; @time rmid12=FlxQTL.getGenoidx(gen12 ,0.01); #7826
gen13=readdlm("../testdata/mouse_hs1940_genochr13.txt")[:,4:end]; @time rmid13=FlxQTL.getGenoidx(gen13 ,0.01); #8308
gen14=readdlm("../testdata/mouse_hs1940_genochr14.txt")[:,4:end]; @time rmid14=FlxQTL.getGenoidx(gen14 ,0.01); #8784
gen15=readdlm("../testdata/mouse_hs1940_genochr15.txt")[:,4:end]; @time rmid15=FlxQTL.getGenoidx(gen15 ,0.01); #9245
gen16=readdlm("../testdata/mouse_hs1940_genochr16.txt")[:,4:end]; @time rmid16=FlxQTL.getGenoidx(gen16 ,0.01); #9700
gen17=readdlm("../testdata/mouse_hs1940_genochr17.txt")[:,4:end]; @time rmid17=FlxQTL.getGenoidx(gen17 ,0.01); #10138
gen18=readdlm("../testdata/mouse_hs1940_genochr18.txt")[:,4:end]; @time rmid18=FlxQTL.getGenoidx(gen18 ,0.01); #10522

X1=FlxQTL.Markers(anno[rmid1,1],anno[rmid1,3],anno[rmid1,4],gen1[rmid1,:])
X2=FlxQTL.Markers(anno[rmid2,1],anno[rmid2,3],anno[rmid2,4],gen2[rmid2,:])
X3=FlxQTL.Markers(anno[rmid3,1],anno[rmid3,3],anno[rmid3,4],gen3[rmid3,:])
X4=FlxQTL.Markers(anno[rmid4,1],anno[rmid4,3],anno[rmid4,4],gen4[rmid4,:])
X5=FlxQTL.Markers(anno[rmid5,1],anno[rmid5,3],anno[rmid5,4],gen5[rmid5,:])
X6=FlxQTL.Markers(anno[rmid6,1],anno[rmid6,3],anno[rmid6,4],gen6[rmid6,:])
X7=FlxQTL.Markers(anno[rmid7,1],anno[rmid7,3],anno[rmid7,4],gen7[rmid7,:])
X8=FlxQTL.Markers(anno[rmid8,1],anno[rmid8,3],anno[rmid8,4],gen8[rmid8,:])
X9=FlxQTL.Markers(anno[rmid9,1],anno[rmid9,3],anno[rmid9,4],gen9[rmid9,:])
X10=FlxQTL.Markers(anno[rmid10,1],anno[rmid10,3],anno[rmid10,4],gen10[rmid10,:])
X11=FlxQTL.Markers(anno[rmid11,1],anno[rmid11,3],anno[rmid11,4],gen11[rmid11,:])
X12=FlxQTL.Markers(anno[rmid12,1],anno[rmid12,3],anno[rmid12,4],gen12[rmid12,:])
X13=FlxQTL.Markers(anno[rmid13,1],anno[rmid13,3],anno[rmid13,4],gen13[rmid13,:])
X14=FlxQTL.Markers(anno[rmid14,1],anno[rmid14,3],anno[rmid14,4],gen14[rmid14,:])
X15=FlxQTL.Markers(anno[rmid15,1],anno[rmid15,3],anno[rmid15,4],gen15[rmid15,:])
X16=FlxQTL.Markers(anno[rmid16,1],anno[rmid16,3],anno[rmid16,4],gen16[rmid16,:])
X17=FlxQTL.Markers(anno[rmid17,1],anno[rmid17,3],anno[rmid17,4],gen17[rmid17,:])
X18=FlxQTL.Markers(anno[rmid18,1],anno[rmid18,3],anno[rmid18,4],gen18[rmid18,:])

@time rmid=FlxQTL.getGenoidx(gen[:,4:end] ,0.01); #10783
X19=FlxQTL.Markers(anno[rmid,1],anno[rmid,3],anno[rmid,4],gen[rmid,4:end]) 

Y=convert(Array{Float64,2},phen');
# m,n=size(Y)

# computing kinship (no loco) and save it to a file (using actual data from GEMMA)
time gemma.3 -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940.pheno.txt -a ../testdata/mouse_hs1940.anno.txt -gk -o kinship -no-check

K=readdlm("../runtest/output/kinship.cXX.txt");
@time T,λ =FlxQTL.K2eig(K+0.00001I); #2.15s

#####whole genescan for comparison w/o loco
@time lod,B,est0= FlxQTL.geneScan(1,T,Matrix(1.0I,m,m),λ,ones(m),Y,X19);
open("../test/HS1940mouse3traits-allscan.txt","w")do io
    writedlm(io, lod)
end

time gemma.3 -g ../testdata/mouse_hs1940.geno.txt.gz -p ../testdata/mouse_hs1940_imp_pheno.txt -n 1 2 3 -k ./output/kinship.cXX.txt -lmm 2 -o lrt3_ch19-v3 -no-check

###### computation times
# all but chr1
function runTime(m,Y)
    flxT=zeros(5,64);
   
    for j=1:64
flxT[1,j] = @elapsed lod,B,est0 =FlxQTL.geneScan(1,T,Matrix(1.0I,m,m),λ,ones(m),Y,X3);
flxT[2,j] = @elapsed lod,B,est0 =FlxQTL.geneScan(1,T,Matrix(1.0I,m,m),λ,ones(m),Y,X7);
flxT[3,j] = @elapsed lod,B,est0 =FlxQTL.geneScan(1,T,Matrix(1.0I,m,m),λ,ones(m),Y,X11);
flxT[4,j] = @elapsed lod,B,est0 =FlxQTL.geneScan(1,T,Matrix(1.0I,m,m),λ,ones(m),Y,X15);
flxT[5,j] = @elapsed lod,B,est0 =FlxQTL.geneScan(1,T,Matrix(1.0I,m,m),λ,ones(m),Y,X19);
    end
    return flxT
    end    


 mtim0= runTime(3,Y);
 median(mtim0,dims=2)
 open("mouse_hs1940_runtime4chr_3traits.txt","w")do io
        writedlm(io,mtim0)
        end

for j=1:3
    Y=readdlm("../testdata/mouse_hs1940_simulated$(3j+3)traits.txt")
        Y=convert(Array{Float64,2},Y');  m=size(Y,1)
        mtime=runTime(m,Y) 
    open("mouse_hs1940_runtime4chr_$(m)traits.txt","w")do io
        writedlm(io,mtime)
        end
     display([m;median(mtime,dims=2)])
end


#Chr1 only 
function fqtl_time(m,Y)
    fqtl=zeros(64);
    for j=1:64
    fqtl[j]= @elapsed lod,B,est0 =FlxQTL.geneScan(1,T,Matrix(1.0I,m,m),λ,ones(m),Y,X1); 
    end
    return fqtl
end


   f= open("mouse_hs1940_chr1runtime_6_9_12traits.txt","w")  
for j=1:3
    Y=readdlm("../testdata/mouse_hs1940_simulated$(3j+3)traits.txt")
        Y=convert(Array{Float64,2},Y');  m=size(Y,1)
        mtime=fqtl_time(m,Y)
        display([m median(mtime)])
       f= open("mouse_hs1940_chr1runtime_6_9_12traits.txt","a")    
        writedlm(f,mtime')
        close(f)
end


#concatenating outputs for file update
chr1=readdlm(homedir()*"/Mountpt/analysis/mouse_hs1940_chr1runtime_6_9_12traits.txt")
    
 idx=[1;3; 7; 11; 15; 19]
fq3=readdlm(homedir()*"/Mountpt/analysis/mouse_hs1940_runtime4chr_3traits.txt")
flx3=[3;median(fq3,dims=2)[idx,1].+2.15 ]
 
fq6=readdlm(homedir()*"/Mountpt/analysis/mouse_hs1940_runtime4chr_6traits.txt")
open("../test/mouse_hs1940_runtime4chr1-3-7-11-15-19_6traits.txt","w")do io
    writedlm(io,[chr1[[1],:];fq6])
end
f6=readdlm("../test/mouse_hs1940_runtime4chr1-3-7-11-15-19_6traits.txt")
flx6=[ 6.0;median(f6,dims=2)[:,1].+2.15]

fq9=readdlm(homedir()*"/Mountpt/analysis/mouse_hs1940_runtime4chr_9traits.txt")
open("../test/mouse_hs1940_runtime4chr1-3-7-11-15-19_9traits.txt","w")do io
    writedlm(io,[chr1[[2],:];fq9])
end
f9=readdlm("../test/mouse_hs1940_runtime4chr1-3-7-11-15-19_9traits.txt")
flx9=[9;median(f9,dims=2)[:,1].+2.15]

fq12=readdlm(homedir()*"/Mountpt/analysis/mouse_hs1940_runtime4chr_12traits.txt")
open("../test/mouse_hs1940_runtime4chr1-3-7-11-15-19_12traits.txt","w")do io
    writedlm(io,[chr1[[end],:];fq12])
end
f12=readdlm("../test/mouse_hs1940_runtime4chr1-3-7-11-15-19_12traits.txt")
flx12=[12;median(f12,dims=2)[:,1].+2.15]

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

Tdata=[flx3 flx6 flx9 flx12 [3;gem3] [6;gem6[idx]] [9;gem9[idx]] [12;gem12[idx]]]
open("../test/runtime-flxqtl_vs_gemma_hsmouse_chr1-3-7-11-15-19.txt","w")do io
    writedlm(io,Tdata)
end

--------
using Plots
pyplot()
using Plots.PlotMeasures


using FlxQTL
#scanned value comparison using -log10(P-value)
flx0=readdlm("../test/HS1940mouse3traits-allscan.txt")
flx= FlxQTL.lod2logP(flx0[:,1] ,3)
gem0=readdlm("../test/lrt3_ch19-v3.assoc.txt",'\t';header=true)

#qqplot
exptd=collect(1:length(flx))
lexp=-log10.(exptd/(length(flx)+1))
fobs=sort(flx;rev=true)
gobs=-log10.(sort(gem0[1][:,end]))
x=collect(1:8)

p1=scatter(lexp,fobs,label="FlxQTL",xlims=(0,8),ylims=(0,8),xlab="Expected -log10(P)",ylab="Observed -log10(P)",aspectratio=1,grid=false,guidefontsize=25,tickfontsize=20,legendfontsize=20, marker=8,
    bottom_margin=0mm,right_margin=0mm)
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
ptime=readdlm("../test/runtime-flxqtl_vs_gemma_hsmouse_chr1-3-7-11-15-19.txt")

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

