using DelimitedFiles, LinearAlgebra
# true parameter values
τ2true=2.0.^-[10:-1:0.1;]

a=0.4; b=0.07
A=a*ones(3,3)+(1-a)*Matrix(1.0I,3,3);B=b*ones(3,3);
Σtrue=[ A B;B A]

#narrow down the range of B
Bfix=zeros(2,2,4);
Bfix[:,:,1]=[-1.0 -1.0;1.0 1.0]

for j=2:4
    Bfix[:,:,j]=Bfix[:,:,1].*sqrt(2)^(j-1)
end
#or
#digit correction
for j=2:4
    Bfix[:,:,j]=round.(Bfix[:,:,1].*sqrt(2)^(j-1),digits=5)
end

Z=[τ2true zeros(10,12)]
for j=1:4
Z[:,j+1]=readdlm("../Result/sim_Z/agrenZ_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt")[:,2]
Z[:,j+5]=readdlm("../Result/sim_Z/agrenZi_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt")[:,2]
end
for j=1:3
Z[:,j+9]=readdlm("../Result/sim_Z/agren_mvlmm_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt")[:,2]   
end
Z[:,end]=[readdlm("../Result/sim_Z/agren_mvlmm_pwr_B$(Bfix[:,:,end])_Σa$(a)_b$(b).txt")[:,2] ;1.0]  
#redraw for τ2true[3]
z0=readdlm("../Result/sim_Z/agrenZ_pwr_τ2true$(τ2true[3])_Σa$(a)_b$(b).txt")[:,2]
z1=readdlm("../Result/sim_Z/agrenZi_pwr_τ2true$(τ2true[3])_Σa$(a)_b$(b).txt")[:,2]
Z[3,2:9]=[z0;z1]


using Plots
pyplot()
using Plots.PlotMeasures

x=round.(Bfix[2,1,:],digits=2)
myColor=["green" "orange" "blue"]

p0 = plot(x,Z[2,2:5], line= (3, :solid), 
    xticks=x,guidefontsize=18,tickfontsize=14, legendfontsize=16,titlefontsize=20,
    xlab="B", ylab="Power", color=myColor[1],right_margin=1cm,grid=false,
    #tickfontrotation=1.5,
     label = "FlxQTL (Z≠I)",title="Power varying with B for τ²= 1/512" ,legend=:bottomright)
scatter!(x, Z[2,2:5],color=myColor[1],  label="",marker=6)

plot!(x,Z[2,6:9], color=myColor[2], line = (3, :solid),label = "FlxQTL (Z=I)")
scatter!(x, Z[2,6:9], color=myColor[2], label="",marker=6)
plot!(x,Z[2,10:end], color=myColor[3], line = (3, :solid),label = "MLMM")
scatter!(x, Z[2,10:end], color=myColor[3], label="",marker=6)



p1 = plot(x,Z[5,2:5], line= (3, :solid), 
    xticks=x,guidefontsize=18,tickfontsize=14, legendfontsize=16,titlefontsize=20,
    xlab="B", ylab="Power", color=myColor[1],right_margin=1cm,grid=false,
    #tickfontrotation=1.5,
     label = "FlxQTL (Z≠I)",title="Power varying with B for τ²= 1/64" ,legend=:bottomright)
scatter!(x, Z[5,2:5],color=myColor[1],  label="",marker=6)

plot!(x,Z[5,6:9], color=myColor[2], line = (3, :solid),label = "FlxQTL (Z=I)")
scatter!(x, Z[5,6:9], color=myColor[2], label="",marker=6)
plot!(x,Z[5,10:end], color=myColor[3], line = (3, :solid),label = "MLMM")
scatter!(x, Z[5,10:end], color=myColor[3], label="",marker=6)

ll = @layout [ a;c]
savefig(plot(p0,p1,layout=ll),"~/GIT/manscripts/gxe/fig/Fig2.eps")




#### Fig_S2
using Gadfly

#generate effect plots from raw data
set_default_plot_size(25cm, 15cm)

tauLab=["1/1024" "1/512" "1/256" "1/128" "1/64" "1/32" "1/16" "1/8" "1/4" "1/2"]
Plots=Array{Plot,1}(undef,length(τ2true));
for j=1:length(τ2true)
    Plots[j] =Gadfly.plot(layer(x=[1.0 1.41421 2.0 2.82843], y=Z[j,2:5],Geom.point,Geom.line,
                Theme(default_color=colorant"green")),
        layer(x=[1.0 1.41421 2.0 2.82843],y=Z[j,6:9],Theme(default_color=colorant"orange"), Geom.point,Geom.line),layer(x=[1.0 1.41421 2.0 2.82843],y=Z[j,10:end],
           Theme(default_color=colorant"blue"),Geom.point,Geom.line),Guide.XLabel("B"),Guide.YLabel("power")
        ,Guide.manual_color_key("τ²= "*tauLab[j],["FlxQTL (Z≠I)","FlxQTL (Z=I)","MLMM" ], ["green","orange","blue"]) );
    
end



set_default_plot_size(25cm, 30cm)
display(gridstack([Plots[1] Plots[2] ;Plots[3] Plots[4];Plots[5] Plots[6];Plots[7] Plots[8];Plots[9] Plots[10]]))

