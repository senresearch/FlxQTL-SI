using DelimitedFiles
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
Z[:,j+1]=readdlm(string(@__DIR__,"/../test/sim_Z/agrenZ_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt"))[:,2]
Z[:,j+5]=readdlm(string(@__DIR__,"/../test/sim_Z/agrenZi_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt"))[:,2]
end
for j=1:3
Z[:,j+9]=readdlm(string(@__DIR__,"/../test/sim_Z/agren_mvlmm_pwr_B$(Bfix[:,:,j])_Σa$(a)_b$(b).txt"))[:,2]   
end
Z[:,end]=[readdlm(string(@__DIR__,"/../test/sim_Z/agren_mvlmm_pwr_B$(Bfix[:,:,end])_Σa$(a)_b$(b).txt"))[:,2] ;1.0]  
#redraw for τ2true[3]
z0=readdlm(string(@__DIR__,"/../test/sim_Z/agrenZ_pwr_τ2true$(τ2true[3])_Σa$(a)_b$(b).txt"))[:,2]
z1=readdlm(string(@__DIR__,"/../test/sim_Z/agrenZi_pwr_τ2true$(τ2true[3])_Σa$(a)_b$(b).txt"))[:,2]
Z[3,2:9]=[z0;z1]
using Gadfly, Plots

#generate effect plots from raw data
set_default_plot_size(25cm, 15cm)

Plots=Array{Plot,1}(undef,length(τ2true));
for j=1:length(τ2true)
    Plots[j] =Gadfly.plot(layer(x=[1.0 1.41421 2.0 2.82843], y=Z[j,2:5],Geom.point,Geom.line,
                Theme(default_color=colorant"purple")),
        layer(x=[1.0 1.41421 2.0 2.82843],y=Z[j,6:9],Theme(default_color=colorant"orange"), Geom.point,Geom.line),layer(x=[1.0 1.41421 2.0 2.82843],y=Z[j,10:end],
           Theme(default_color=colorant"blue"),Geom.point,Geom.line),Guide.XLabel("B"),Guide.YLabel("power")
        ,Guide.manual_color_key("τ^2=$(τ2true[j])",["Z(≠I): contrasts","Z=I","MLMM" ], ["purple","orange","blue"]) );
    
end



set_default_plot_size(25cm, 30cm)
display(gridstack([Plots[1] Plots[2] ;Plots[3] Plots[4];Plots[5] Plots[6];Plots[7] Plots[8];Plots[9] Plots[10]]))

