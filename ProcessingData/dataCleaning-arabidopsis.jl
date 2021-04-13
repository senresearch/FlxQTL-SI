
##detect idential markers and drop one of them to create files for 2d-scan
using DelimitedFiles

labels=readdlm("../processedData/agren_markers_multipleqtl.csv",',';skipstart=1)
X=readdlm("../processedData/agren_gene_imp.csv",',';skipstart=1)'
Chr=unique(labels[:,2])
nChr=length(Chr)

#detect identical markers
for i=1:nChr
    maridx=findall(labels[:,2].==Chr[i]); 
    for j=maridx[1]:last(maridx)-1
        for l=j+1:last(maridx)
            if(X[j,:]==X[l,:])
                println([i j l])
            end
        end
    end
end


# filter identical markers and produce new data 
# (Note that 2 extra bad markers are dropped when genome scan gives any error )
Xnew=vcat(X[1:16,:],X[18:51,:],X[[53],:],X[55:118,:],X[120:142,:],X[144:152,:],X[154:168,:],X[170:205,:],X[207:226,:],
    X[228:247,:],X[249:263,:],X[265:328,:],X[330:336,:],X[338:339,:],X[341:347,:],X[349:358,:],X[360:364,:],X[366:368,:],
    X[370:376,:],X[378:405,:],X[407:end,:]);
newlabel=vcat(labels[1:16,:],labels[18:51,:],labels[[53],:],labels[55:118,:],labels[120:142,:],labels[144:152,:],
    labels[154:168,:],labels[170:205,:],labels[207:226,:],labels[228:247,:],labels[249:263,:],
    labels[265:328,:],labels[330:336,:],labels[338:339,:],labels[341:347,:],labels[349:358,:],labels[360:364,:],
    labels[366:368,:],labels[370:376,:],labels[378:405,:],labels[407:end,:]);

writedlm("../processedData/agren_newmarkers_multipleqtl.csv",newlabel,',');
writedlm("../processedData/agren_gene_imp_revision.csv",Xnew,',');

### climate data processing
#weather data from 7/1/2009 to 6/30/2012
weather_it=readdlm("../rawData/agren_climate_itCastelnuov2009_2012.csv",',';header=true); 
# variable order: min mean max
weather_sw=readdlm("../rawData/agren_climate_swRodasen2009_2012.csv",',';header=true);
# variable order: min max mean

#remove information on 2/29/2012 : 1095 days by 6 environments
soil_it=weather_it[1][1:end.!=974,8]-weather_it[1][1:end.!=974,6];
soil_sw=weather_sw[1][1:end.!=974,7]-weather_sw[1][1:end.!=974,6];
SDR= [soil_it[1:365,:] soil_it[366:730,:] soil_it[731:end,:] soil_sw[31:395,:] soil_sw[396:760,:] soil_sw[761:end,:]]
writedlm("../processedData/agren_soilDailyRng_it_sw_09-12.txt",SDR,',')


