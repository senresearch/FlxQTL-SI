using DelimitedFiles, LinearAlgebra
# cases for true parameters
τ2true=2.0.^-[10:-1:0.1;]

### Sigma (default) for pwr0
# a=0.4; b=0.07 
# A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);B=b*ones(3,3);
# Σtrue=[ A B;B A]
###


### power varying with B for each τ2
Btrue=zeros(2,2,10)
Btrue[:,:,1]= [0.25  -0.25 ;0.25   0.25]
for j=2:10
    Btrue[:,:,j].=round.(Btrue[:,:,1]*sqrt(2)^(j-1),digits=5)
end

pwr0=[τ2true zeros(10,10)]
for j=1:10
# pwr0[:,j+1]=readdlm(string(@__DIR__,"/../test/sim_Kc/agren_pwr_$(j)th_B$(Btrue[:,:,j]).txt"))[:,2]
pwr0[:,j+1]=readdlm("../Result/sim_Kc/agren_pwr_$(j)th_B$(Btrue[:,:,j]).txt")[:,2]
end


#Kc vs Kc=I for Σtrue (narrowed down B)
δ=sqrt(2).^[3.:0.5:5.;]
B1=zeros(2,2,length(δ));
for j=1:length(δ)
    B1[:,:,j]=round.(Btrue[:,:,1]*δ[j],digits=5)
end
    
pwr1=[τ2true zeros(10,10)]
for j=1:length(δ)
# pwr1[:,j+1]= readdlm(string(@__DIR__,"/../test/sim_Kc/agren1_pwr_$(j)th_B$(B1[:,:,j]).txt"))[:,2]
# pwr1[:,j+6]= readdlm(string(@__DIR__,"/../test/sim_Kc/agren2_I_pwr_$(j)th_B$(B1[:,:,j]).txt"))[:,2]
pwr1[:,j+1]= readdlm("../Result/sim_Kc/agren1_pwr_$(j)th_B$(B1[:,:,j]).txt")[:,2]
pwr1[:,j+6]= readdlm("../Result/sim_Kc/agren2_I_pwr_$(j)th_B$(B1[:,:,j]).txt")[:,2]
end

#Kc vs Kc=I for different Σtrue for B fixed

### Sigma for pwr2
# Sigma (default) for pwr2[:,2:3]

#  (a,b)=(0.3,0.07) for pwr2[:,4:5]
# A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);B=b*ones(3,3);
# Σtrue=[ A B;B A]

# a=0.3; b=0.0; c=0.1 for pwr2[:,6:7]
# A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);C=c*ones(3,3)+(1-c)*Matrix(1.0I,3,3);B=zeros(3,3)
# Σtrue=[A B;B C]

###### Aggregate simulated results into one file
Bt=round.([sqrt(2) -sqrt(2);sqrt(2) sqrt(2)],digits=5)
pwr2=[τ2true zeros(10,9)]

a=0.3;c=0.1
pwr2[:,2]= readdlm("../Result/sim_Kc/agren_pwr_B$(Bt).txt")[:,2]
pwr2[:,3]= readdlm("../Result/sim_Kc/agrenI_pwr_B$(Bt).txt")[:,2]
pwr2[:,4]= readdlm("../Result/sim_Kc/agren_pwr_B$(Bt)_Σa$(a).txt")[:,2]
pwr2[:,5]= readdlm("../Result/sim_Kc/agrenI_pwr_B$(Bt)_Σa$(a).txt")[:,2]    
pwr2[:,6]= readdlm("../Result/sim_Kc/agren_pwr_B$(Bt)_Σa$(a)bc$(c).txt")[:,2]
pwr2[:,7]= readdlm("../Result/sim_Kc/agrenI_pwr_B$(Bt)_Σa$(a)bc$(c).txt")[:,2]
pwr2[:,8]= readdlm("../Result/sim_Kc/agren_pwr_B$(Bt)_Σ=I.txt")[:,2]
pwr2[:,9]= readdlm("../Result/sim_Kc/agrenI_pwr_B$(Bt)_Σ=I.txt")[:,2]
pwr2[:,10]= readdlm("../Result/sim_Kc/agren_pwr_B$(Bt)_H(0.19)+Σ=I.txt")[:,2] # not included in fig

#generate a result file for pwr2
open("../Result/sim_Kc/agren_pwr4Sigmas.txt","w")do io
    writedlm(io,pwr2)
end


# using Gadfly,Colors

# set_default_plot_size(30cm, 30cm)

# # power by τ^2
# layer1=Array{Array{Layer,1}}(undef,5);layer2=Array{Array{Layer,1}}(undef,5)
# # color1=distinguishable_colors(5,[RGB24(0.8,0.5,0.2)]);color2=distinguishable_colors(5,[RGB24(0.8,0.5,0.2)])
# color1=["blue","purple","green","orange","grey"];
# x=convert(Array{Float64},Btrue[1,1,:]')
# for j=1:5
#     layer1[j]=layer(x=x,y=pwr0[j,2:end],Geom.point,Geom.line,Theme(default_color=color1[j]))
#     layer2[j]=layer(x=x,y=pwr0[j+5,2:end],Geom.point,Geom.line,Theme(default_color=color1[j]))
# end

# a1=plot(layer1...,Guide.XLabel("B"),Guide.YLabel("power"),Guide.title("Power varying with B for each τ^2"),Guide.manual_color_key("τ^2",string.(τ2true)[1:5],color1));

# a2=plot(layer2...,Guide.XLabel("B"),Guide.YLabel("power"),Guide.title("Power varying with B for each τ^2"),Guide.manual_color_key("τ^2",string.(τ2true)[6:end],color1));
# vstack(a1,a2)

##### # # power by τ^2
using Plots
#pyplot()
using Plots.PlotMeasures


B0=round.(Btrue[1,1,:],digits=2)[1:end-1]
tauLab=["1/1024" "1/512" "1/256" "1/128" "1/64" "1/32" "1/16" "1/8" "1/4" "1/2"]
myColor=[:orange :lightgreen :darkgreen :blue :purple]        

p0 = plot(B0,pwr0[1,2:end-1], line= (2, :solid), 
    xticks=B0,#guidefontsize=18,tickfontsize=8.6, legendfontsize=16,titlefontsize=20,xrotation=15,
    xrotation=55,
    xlab="B", ylab="Power", color=myColor[1],right_margin=1cm,grid=false,
    #tickfontrotation=1.5,
     label = tauLab[1],title="Power varying with B for each τ²" ,legend=:bottomright)
scatter!(B0, pwr0[1,2:end-1],color=myColor[1],  label="",marker=3)

for j in 2:5
    plot!(B0,pwr0[j,2:end-1], color=myColor[j], line = (2, :solid), 
          xlab="B", ylab="Power",label = tauLab[j])
    scatter!(B0, pwr0[j,2:end-1], color=myColor[j], label="",marker=3)
end
p0

p1 = plot(B0,pwr0[6,2:end-1], line= (2, :solid), 
    xticks=B0,#guidefontsize=9,tickfontsize=8.6, legendfontsize=10,titlefontsize=13,
    xrotation=55,
    xlab="B", ylab="Power", color=myColor[1],right_margin=1cm,grid=false,
     label = tauLab[6] ,legend=:bottomright)
scatter!(B0, pwr0[6,2:end-1],color=myColor[1],  label="",marker=3)

for j in 2:5
    plot!(B0,pwr0[j+5,2:end-1], color=myColor[j], line = (2, :solid), 
          xlab="B", ylab="Power",label = tauLab[j+5])
    scatter!(B0, pwr0[j+5,2:end-1], color=myColor[j], label="",marker=3)
end
p1

ll = @layout [ a;c]
#plot(p0,p1,layout=ll)
savefig(plot(p0,p1,layout=ll),"~/GIT/manscripts/gxe/fig/Fig3.eps")


### power by Kc=I vs Kc!=I (Fig_S1.png)
using Gadfly,Colors
Plots=Array{Plot,1}(undef,length(τ2true));ticks = [0.0,0.2, 0.4, 0.6,0.8,1.0]
for j=1:length(τ2true)
    Plots[j] =plot(layer(x=pwr1[j,2:6], y=pwr1[j,7:end],Geom.point,Geom.line,Geom.abline(color="red", style=:dash),
                Theme(default_color=colorant"blue")),Guide.XLabel("power (Kc(≠I))"),Guide.YLabel("power (Kc=I)",orientation=:vertical)
      ,Guide.yticks(ticks=ticks),Guide.xticks(ticks=ticks)
        ,Guide.title("Power varied by B for τ²= "*tauLab[j]));
    
end

set_default_plot_size(25cm, 30cm)
display(gridstack([Plots[1] Plots[2] ;Plots[3] Plots[4];Plots[5] Plots[6];Plots[7] Plots[8];Plots[9] Plots[10]]))

#power by Sigmas (Fig_S2.png)
pwr2=readdlm("../Result/sim_Kc/agren_pwr4Sigmas.txt")

layers1=Array{Array{Layer,1}}(undef,4);
colors1=["blue","purple","green","orange"];
ticks = [0.0,0.2, 0.4, 0.6,0.8,1.0]

for j=1:4    layers1[j]=layer(x=pwr2[:,2j],y=pwr2[:,2j+1],Geom.point,Geom.line,Theme(default_color=colors1[j]),Geom.abline(color="grey", style=:dash))
end
set_default_plot_size(13cm, 13cm)
Gadfly.plot(layers1...,Guide.XLabel("power for Kc(≠I)"),Guide.YLabel("power for Kc=I")
      ,Guide.yticks(ticks=ticks),Guide.xticks(ticks=ticks)
        ,Guide.title("Power varied by τ² for B=[√2 -√2; √2 √2])"),Guide.manual_color_key("Σ",["1","2","3","I"],colors1))


#Σ
# 1. a=0.4; b=0.07
# A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);B=b*ones(3,3);
# Σtrue=[ A B;B A]
# 2. a=0.3 the rest unchanged
# 3. a=0.3; b=0.0; c=0.1
# A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);C=c*ones(3,3)+(1-c)*Matrix(1.0I,3,3);B=zeros(3,3)
# Σtrue=[A B;B C]
#pwr2 : Kc!=I vs Kc=I for a pair of columns(1,2,3,I)
#  0.000976563  0.726  0.727  0.675  0.674  0.369  0.375  0.481  0.485
#  0.00195313   0.62   0.617  0.537  0.538  0.492  0.5    0.676  0.701
#  0.00390625   0.691  0.712  0.633  0.628  0.431  0.444  0.59   0.597
#  0.0078125    0.564  0.558  0.488  0.495  0.181  0.197  0.522  0.545
#  0.015625     0.523  0.524  0.745  0.739  0.633  0.659  0.66   0.653
#  0.03125      0.603  0.62   0.51   0.519  0.468  0.475  0.644  0.651
#  0.0625       0.657  0.643  0.557  0.554  0.382  0.392  0.63   0.644
#  0.125        0.764  0.75   0.599  0.605  0.514  0.526  0.55   0.562
#  0.25         0.754  0.763  0.615  0.656  0.606  0.613  0.49   0.488
#  0.5          0.331  0.317  0.537  0.545  0.483  0.495  0.579  0.609



