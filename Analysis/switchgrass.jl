using DelimitedFiles, StatsBase, Statistics, LinearAlgebra, Distributed
@everywhere using FlxQTL
FlxQTL.setSeed(3,1000)

#geno prob & trait data
@time geno=readdlm("../processedData/switchgrass36site_yr_geno_prob_imputed.csv",',';skipstart=1);
marker=readdlm("../processedData/switchgrass36site_yr_markers.csv",',';skipstart=1); #6118 markers x 4 genotype prob
XX=FlxQTL.Markers(marker[:,1],marker[:,2],marker[:,3],geno')

pheno=readdlm("../36_site_yr_fl50/FL50_36_Site_Yr.csv",',';skipstart=1)[:,2:end];
Y3=pheno.-mean(pheno,dims=1)
Y3=convert(Array{Float64,2},(Y3./std(Y3))')

m,n=size(Y3) #(36, 750)
# environment data: photo period only 
photo=readdlm("../processedData/DYLN_MONTH_FL50_36_Site_Yr.csv",',';skipstart=1)'

#standarization
Z_pht=(photo.-mean(photo))
Z_pht=Z_pht./std(Z_pht)

#genetic kinship
@time K=FlxQTL.kinshipLoco(FlxQTL.kinshipLin,4,XX);
@time T, Λ=K2eig(K,true);

#scan w/o Kc
@time  lod1,B1,est01 =FlxQTL.geneScan(4,T,Matrix(1.0I,m,m),Λ,ones(m),Y3,XX,true;ρ=0.01);

# the second outcome (lod) is for Kc=Gaussian (no part of the analysis in the manuscript)
open("../Result/switchgrass_fl50_36site_yr_Z=Kc=I-Kc=photo_Gauss-colcent-overstd.txt","w")do io
    writedlm(io, [lod1 lod])
end


# # permutation test
# Kg= QTL.kinshipLin(XX.X,4);
# @time maxLODs, cutoff = FlxQTL.permTest(1000,4,Kg,Matrix(1.0I,m,m),Y3,XX,ones(m,1);ρ=0.01);
# cutoff # 4.85103  5.50824


#plot for 1d scan: Fig6.eps
fl50=readdlm("../Result/switchgrass_fl50_36site_yr_Z=Kc=I-Kc=photo_Gauss-colcent-overstd.txt")

fl=layers(marker[:,2],marker[:,3],fl50[:,[1]])
# plot1d(fl;title="switchgrass 10 sites by 4 yrs (36 traits,Kc=I, Z=intercept) @ α = 0.1, 0.05, 0.01",yint=[4.4575 4.85103  5.50824])
plot1d(fl;title="FL50: switchgrass 10 sites by 4 yrs (36 traits,Z=Kc=I) ")

##select 3 major qtl
# q1=findall(fl50[:,1].==maximum(fl50)); # max:3398 

# idx5N=findall(marker[:,2].=="5N")#3369~3826 
# #2nd max
# findall(fl50[idx5N,:].==minimum(fl50[idx5N,:])) #207
# findall(fl50[idx5N[207:end],:].==maximum(fl50[idx5N[207:end],:])) #175

# #3rd max
# idx4K=findall(marker[:,2].=="4K") #2389~2647
# findall(fl50[idx4K,1].==maximum(fl50[idx4K,1])) #107

q1=3398;q2= 3749;q3=2495
Q=[q1;q2;q3]

#envirionmental scan
@time e_pht,B_p,est0_p = QTL.envScan(Q,4,T,Matrix(1.0I,m,m),Λ,ones(m),Y2,XX,Z_pht,true);e_pht
# @time e_max,B_x,est0_x = QTL.envScan(Q,4,T,Matrix(1.0I,m,m),Λ,ones(m),Y2,XX,Z_max,true); e_max
# @time e_min,B_m,est0_m = QTL.envScan(Q,4,T,Matrix(1.0I,m,m),Λ,ones(m),Y2,XX,Z_min,true); e_min
# @time e_ra,B_r,est0_r = QTL.envScan(Q,4,T,Matrix(1.0I,m,m),Λ,ones(m),Y2,XX,Z_ra,true); e_ra


open("../Result/switchgrass_fl50_centrl_36site_yr_Kc=I_dynl_tmax_min_rain.txt","w") do io
    writedlm(io,[e_pht e_max e_min e_ra])
end

using Plots
pyplot()
using Plots.PlotMeasures

x=collect(1:1:12)
myColor=["green" "orange" "blue"]

p0 = plot(x,Es[:,1], line= (3, :solid), 
    xticks=x,guidefontsize=18,tickfontsize=14, legendfontsize=16,titlefontsize=20,
    xlab="Month", ylab="LOD", color=myColor[1],right_margin=1cm,grid=false,
    #tickfontrotation=1.5,
     label = "Chr$(marker[Q[1],2])@$(round(marker[Q[1],3],digits=3))",title="Switchgrass FL50: monthly photoperiod",legend=:right )
scatter!(x, Es[:,1],color=myColor[1],  label="",marker=6)
for j=2:3
plot!(x,Es[:,j], color=myColor[j], line = (3, :solid),label = 
        "Chr$(marker[Q[j],2])@$(round(marker[Q[j],3],digits=3))")
scatter!(x, Es[:,j], color=myColor[j], label="",marker=6)
end

p0
savefig(p0,"~/GIT/manscripts/gxe/fig/Fig7.eps")

# using Gadfly

# #centralization -> standardized by the centered values
# Es=readdlm("../Result/switchgrass_fl50_centrl_36site_yr_Kc=I_dynl_tmax_min_rain.txt")
    
# #photoperiod only
# set_default_plot_size(15cm, 10cm)

# p=Gadfly.plot(layer(y=Es[:,1],Geom.point,Geom.line,Theme(default_color=color("blue"))),layer(y=Es[:,2],Geom.point,Geom.line,Theme(default_color=color("orange"))),layer(y=Es[:,3],Geom.point,Geom.line,Theme(default_color=color("purple"))),Guide.manual_color_key("Legend",["Chr$(marker[Q[1],2])@$(round(marker[Q[1],3],digits=3))","Chr$(marker[Q[2],2])@$(round(marker[Q[2],3],digits=3))","Chr$(marker[Q[3],2])@$(round(marker[Q[3],3],digits=3))" ], ["blue","orange","purple" ]),
#   Guide.title("Switchgrass FL50: Montly Photoperiod"),
#     Guide.YLabel("LOD",orientation=:vertical),Guide.XLabel("Month"),Coord.cartesian(xmin=0, xmax=12))

# # p |> PNG(homedir()*"/GIT/manscripts/gxe/fig/sw-fl50-environscan.png", 30cm, 18cm)


