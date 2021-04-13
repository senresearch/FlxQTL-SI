#using ColorBrewer
#using DataFrames
using Colors
using Gadfly

struct Effect
  chr::Array{Any,1}
  pos::Array{Float64,1}
  eff::Array{Float64,1}
end  

function chrlabel(i::Int64)
    return chrs[i]
end

function effectplot(x::Effect;xlabel="Chromosome",ylabel="LOD scores",title=" ",ymin=minimum(x.eff)-1.5, ymax=maximum(x.eff)+1.5,yint=[], 
    color="grey",style=:dot)


    # chromosomes
    chrs = unique( x.chr )
    # number of chromosomes
    nChr = length( chrs )
    # chromosome lengths
    chrlens = zeros(nChr)
    # number of markers
    nMar = length(x.pos)

    # max/min of effects
    effectMax = maximum(x.eff)
    effectMin = minimum(x.eff)

    # find chromosome lengths
    cumpos = zeros(nMar)
    cumlen = 0.0
    pl = Array{Array{Layer,1}}(undef,nChr)
    for i=1:nChr
        chridx = (x.chr.==chrs[i])
        cumpos[chridx] = cumlen .+ x.pos[chridx]
        chrlens[i] = maximum(x.pos[chridx])
        cumlen += chrlens[i]
        pl[i] = layer(x=cumpos[chridx],y=x.eff[chridx],Geom.line,Theme(default_color="green"))
    end
 
    plot(pl...,Guide.xlabel(xlabel), yintercept=yint,Geom.hline(color=color,style=style),
       Guide.title(title),Guide.YLabel(ylabel),Coord.cartesian(ymin=ymin, ymax=ymax),Guide.xrug,
         Guide.xticks(ticks=cumsum(chrlens)))
  
end


# set_default_plot_size(25cm, 20cm)    
# x = Effect([1, 1, 1, 2, 2, 3, 3], [0.0, 10.0, 30.0, 0.0, 30.0, 2.0, 20.0], randn(7,1));

# https://timothyrenner.github.io/programming/2015/02/01/gadfly-splats-and-layers.html


### Generate multiple layers of effectplots (all effect plots have different colors)
##Input :
## x : a type of Mlayers including chromosome, position, and effect or lod score information.
## effect : default is true, i.e. create an effect plot. If effect=false, create a genenome scan plot.
struct Mlayers
    chr::Array{Any,1}
    pos::Array{Float64,1}
    effs::Array{Float64,2}
end


function meffectplot(x::Mlayers;xlabel="Chromosome",ylabel::String="y", title::String=" ",yint=[]
        ,color="grey",style=:dot,ymin=minimum(x.effs), ymax=maximum(x.effs),legend=["col1","col"])
    #chromosomes
    chrs=unique(x.chr)
    # #of chromosomes (for xticks)
    nChr=length(chrs)
    
    nper=size(x.effs,2)
   clength=zeros(nChr); cumpos=zeros(length(x.pos)); cumlen=zeros(nper);
     ply=Array{Array{Layer,1}}(undef,nChr*nper)
  #color seeds     
    colors=distinguishable_colors(nper,[RGB24(0.8,0.5,0.2)])
    for j=1:nper
        for i=1:nChr
        cidx=(x.chr.==chrs[i])
        cumpos[cidx]=cumlen[j] .+x.pos[cidx]
        clength[i]=maximum(x.pos[cidx])
        cumlen[j] +=clength[i]
        ply[nChr*(j-1)+i]= layer(x=cumpos[cidx],y=x.effs[cidx,j],Geom.line,Theme(default_color=colors[j]))
        end    
    end
        ticks=cumsum(clength)
    
    plot(ply...,Guide.XLabel(xlabel),Guide.YLabel(ylabel),Guide.xticks(ticks=ticks),Guide.title(title),Guide.xrug,
        Coord.cartesian(ymin=ymin, ymax=ymax),yintercept=yint,Geom.hline(color=color,style=style)
            ,Guide.manual_color_key("Legend",legend , colors))
      
  
   
end
 
## two color option with legend: x.effs[:,1] is colored differently from the rest.
function effect2plot(x::Mlayers;xlabel="Chromosome",ylabel::String="y", title::String=" ",yint=[]
        ,color="grey",style=:dot,ymin=minimum(x.effs), ymax=maximum(x.effs),colors=["black","skyblue"]
        ,legend=["black","skyblue"])
    #chromosomes
    chrs=unique(x.chr)
    # #of chromosomes (for xticks)
    nChr=length(chrs)
    
    nper=size(x.effs,2)
   clength=zeros(nChr); cumpos=zeros(length(x.pos)); cumlen=zeros(nper);
     ply=Array{Array{Layer,1}}(undef,nChr*nper)
  #color seeds     
    colors=colors
    for i=1:nChr
        cidx=(x.chr.==chrs[i])
        cumpos[cidx]=cumlen[1] .+x.pos[cidx]
        clength[i]=maximum(x.pos[cidx])
        cumlen[1] +=clength[i]
        ply[i]= layer(x=cumpos[cidx],y=x.effs[cidx,1],Geom.line,Theme(default_color=colors[1]))
     end    
    
    for j=2:nper
        for i=1:nChr
        cidx=(x.chr.==chrs[i])
        cumpos[cidx]=cumlen[j] .+x.pos[cidx]
        clength[i]=maximum(x.pos[cidx])
        cumlen[j] +=clength[i]
        ply[nChr*(j-1)+i]= layer(x=cumpos[cidx],y=x.effs[cidx,j],Geom.line,Theme(default_color=colors[2]))
        end    
    end
    
        ticks=cumsum(clength)
    
    plot(ply...,Guide.XLabel(xlabel),Guide.YLabel(ylabel),Guide.xticks(ticks=ticks),Guide.title(title),Guide.xrug,
        Coord.cartesian(ymin=ymin, ymax=ymax),yintercept=yint,Geom.hline(color=color,style=style)
            ,Guide.manual_color_key("Legend",legend , colors))
           
end
 
## two color option w/o legend:
function effectplot(x::Mlayers;xlabel="Chromosome",ylabel::String="y", title::String=" ",yint=[]
        ,color="grey",style=:dot,ymin=minimum(x.effs), ymax=maximum(x.effs),colors=["black","skyblue"])
    #chromosomes
    chrs=unique(x.chr)
    # #of chromosomes (for xticks)
    nChr=length(chrs)
    
    nper=size(x.effs,2)
   clength=zeros(nChr); cumpos=zeros(length(x.pos)); cumlen=zeros(nper);
     ply=Array{Array{Layer,1}}(undef,nChr*nper)
  #color seeds     
    colors=colors
    for i=1:nChr
        cidx=(x.chr.==chrs[i])
        cumpos[cidx]=cumlen[1] .+x.pos[cidx]
        clength[i]=maximum(x.pos[cidx])
        cumlen[1] +=clength[i]
        ply[i]= layer(x=cumpos[cidx],y=x.effs[cidx,1],Geom.line,Theme(default_color=colors[1]))
     end    
    
    for j=2:nper
        for i=1:nChr
        cidx=(x.chr.==chrs[i])
        cumpos[cidx]=cumlen[j] .+x.pos[cidx]
        clength[i]=maximum(x.pos[cidx])
        cumlen[j] +=clength[i]
        ply[nChr*(j-1)+i]= layer(x=cumpos[cidx],y=x.effs[cidx,j],Geom.line,Theme(default_color=colors[2]))
        end    
    end
    
        ticks=cumsum(clength)
    
    plot(ply...,Guide.XLabel(xlabel),Guide.YLabel(ylabel),Guide.xticks(ticks=ticks),Guide.title(title),Guide.xrug,
        Coord.cartesian(ymin=ymin, ymax=ymax),yintercept=yint,Geom.hline(color=color,style=style))
           
end

### generate multiple plot layers only to customize one's own plot.

function mplLayers(x::Mlayers)
    #chromosomes
    chrs=unique(x.chr)
    # #of chromosomes (for xticks)
    nChr=length(chrs)
    
    nper=size(x.effs,2)
   clength=zeros(nChr); cumpos=zeros(length(pos)); cumlen=zeros(nper);
   # colors=palette("Paired",nrep)  in colorbrewer
    #color in colors
    colors=distinguishable_colors(nper,[RGB24(0.8,0.5,0.2)])
    ply=Array{Array{Layer,1}}(undef,nChr*nper)
    for j=1:nper
        for i=1:nChr
            
        cidx=(x.chr.==chrs[i])
        cumpos[cidx]=cumlen[j] .+x.pos[cidx]
        clength[i]=maximum(x.pos[cidx])
        cumlen[j] +=clength[i]
        ply[nChr*(j-1)+i]= layer(x=cumpos[cidx],y=x.effs[cidx,j],Geom.line,Theme(default_color=colors[j]))
        end
    end
  return  ticks=cumsum(clength), ply
 #   plot(ply...,Guide.XLabel("Chromosome"),Guide.YLabel("Effect size"),Guide.xticks(ticks=ticks),Guide.title("Effect Plots"))
       # ,Coord.cartesian(xmin=0, xmax=ticks[end]))
end
     
    
    
    
    
    
    
    
    
    
    
    
    