clear all
close all
%clc
rng('shuffle');
addpath(['~/My Drive (john.c.rollins@gmail.com)//utils/'])
addpath(['~/My Drive (john.c.rollins@gmail.com)/utils/GeodeticToolbox/'])
addpath(['~/My Drive (john.c.rollins@gmail.com)/utils/surfacevel2strain/matlab/util_basic/'])
addpath(['~/My Drive (john.c.rollins@gmail.com)/utils/surfacevel2strain/matlab/util_est/'])
addpath(['~/My Drive (john.c.rollins@gmail.com)/utils/surfacevel2strain/matlab/util_euler/'])
addpath(['~/My Drive (john.c.rollins@gmail.com)/utils/surfacevel2strain/matlab/misc/'])

load colormaps.mat
coastlines = readtable(['~/My Drive (john.c.rollins@gmail.com)//utils/NZcoastlines/nz-coastline-mean-high-water/nz-coastline-mean-high-water_180.dat']);
lakes = readtable(['~/My Drive (john.c.rollins@gmail.com)//utils/NZcoastlines/nz-lake-polygons-topo-1500k/nz-lake-polygons-topo-1500k.dat']);

load(['~/My Drive (john.c.rollins@gmail.com)/NZ/Material/faults/CFM_2022/cfm_1pt0_vertices_NaNs.txt']);
cfm_1pt0_vertices_NaNs(cfm_1pt0_vertices_NaNs(:,1)<0,1) = cfm_1pt0_vertices_NaNs(cfm_1pt0_vertices_NaNs(:,1)<0,1) + 360;

addpath(['mbin/gps_ec'])
addpath(['mbin/stats'])

tol = 20.0;
minocctime = 2.49;
combdist = 1;

lonbounds = [90 270];
latbounds = [-55 -10];

E_b16 = [37.66	32.34  -0.635];
OM_b16 = eul2rot(E_b16)';

%%
beavan16 = readtable('./horiz/Beavan_2016/beavan16_flipped.dat');
for sitename = 1:length(beavan16.x_E)
    beavan16.sitecode{sitename} = beavan16.SITE{sitename}(1:4);
end
%return
%%
beavan12 = readtable('./Verticals/Beavan_2012/beavan12_awk.dat');
beavan12.lon(beavan12.lon<0) = beavan12.lon(beavan12.lon<0) + 360;
beavan12.occtime = load('./Verticals/Beavan_2012/beavan12_occtime.dat');

%%
houlie17 = readtable('./Verticals/Houlie_Stern_2017/houlie17.dat');
houlie17.lon(houlie17.lon<0) = houlie17.lon(houlie17.lon<0) + 360;
houlie17.occtime = houlie17.totalextent/365.243;
houlie17.nocc = houlie17.numdays;
houlie17.prcofdata = houlie17.numdays./houlie17.totalextent;

srwn = 1./sqrt(houlie17.occtime);
swn = 2.*sqrt(3).*1.1./(sqrt(houlie17.nocc).*houlie17.occtime);
%return
houlie17.su = max([houlie17.su, sqrt(srwn.^2 + swn.^2)/2 + houlie17.su/2],[],2);
houlie17.su = sqrt(houlie17.su.^2 + 0.3^2);

houlie17.vu_unadj = houlie17.vu - 2;

%%
ffile = fopen(['./horiz/midas/midas.IGS14.txt']);
midasinput = textscan(ffile,'%s %*s %f %f %f %*f %f %*f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %*f\n');
fclose(ffile);

midas.lons = midasinput{13} + 360;
midas.lats = midasinput{12};
midas.names = midasinput{1};
midas.starttime = midasinput{2};
midas.endtime = midasinput{3};
midas.occtime = midasinput{4};
midas.nocc = midasinput{5};
midas.veast = midasinput{6}*1e3;
midas.vnorth = midasinput{7}*1e3;
midas.vup = midasinput{8}*1e3;
midas.seast = midasinput{9}*1e3;
midas.snorth = midasinput{10}*1e3;
midas.su = midasinput{11}*1e3;
midas.coven = 0*midas.veast;
midas = struct2table(midas);

midas.inbounds = and(and(midas.lons>=lonbounds(1),midas.lons<=lonbounds(2)),and(midas.lats>=latbounds(1),midas.lats<=latbounds(2)));
midas.inbounds(or(or(midas.seast>tol,midas.snorth>tol),midas.su>tol)) = 0;
midas(midas.inbounds==0,:) = [];
midas.prcofdata = midas.nocc./(midas.occtime.*365.243);

srwn = 1./sqrt(midas.occtime);
swn = 2.*sqrt(3).*1.1./(sqrt(midas.nocc).*midas.occtime);
%return
midas.su = max([midas.su, midas.su/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
midas.seast = max([midas.seast, midas.seast/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
midas.snorth = max([midas.snorth, midas.snorth/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
%return
O = [midas.lats,midas.lons,zeros(length(midas.lons),1)];
V = [midas.vnorth,midas.veast,midas.vup]*1e-3;
SV = [midas.snorth,midas.seast,midas.su]*1e-3;

[x,y,z] = wgs2xyz(midas.lons,midas.lats,zeros(length(midas.lats),1));
[vxyz,~] = neu2xyz(O,V,SV,0*SV);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2014',[],'ITRF2008',[],[]);
[itrfout] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
midas.vnorth_itrf08 = itrfout(:,1)*1e3;
midas.veast_itrf08 = itrfout(:,2)*1e3;
midas.vup_itrf08 = itrfout(:,3)*1e3;

[x,y,z] = wgs2xyz(midas.lons,midas.lats,zeros(length(midas.lats),1));
[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
midas.veast_b16 = midas.veast_itrf08 + Venu(:,1);
midas.vnorth_b16 = midas.vnorth_itrf08 + Venu(:,2);
%return
%%
ulr5 = readtable('./Verticals/ULR5/ulr5_mod.txt');
ulr5.Lon(ulr5.Lon<0) = ulr5.Lon(ulr5.Lon<0) + 360;
ulr5.inbounds = and(and(ulr5.Lon>=lonbounds(1),ulr5.Lon<=lonbounds(2)),and(ulr5.Lat>=latbounds(1),ulr5.Lat<=latbounds(2)));
ulr5(ulr5.inbounds==0,:) = [];
ulr5.prcofdata = ulr5.Data/100;
ulr5.nocc = ulr5.T_GPS.*365.243.*ulr5.prcofdata;

%%
holden15 = readtable('./Verticals/Holden_2015/holden15_sed4.dat');
holden15.Site = upper(holden15.Site);
%holden15(or(holden15.sVUP>tol,or(holden15.sVE>tol,holden15.sVN>tol)),:) = [];
%return
holden15.occtime(strcmp(holden15.Site,'1291')) = 11;
holden15.nocc(strcmp(holden15.Site,'1291')) = 8;

holden15.occtime(strcmp(holden15.Site,'2308')) = 9;
holden15.nocc(strcmp(holden15.Site,'2308')) = 9;

holden15.occtime(strcmp(holden15.Site,'2309')) = 11;
holden15.nocc(strcmp(holden15.Site,'2309')) = 8;

holden15.occtime(strcmp(holden15.Site,'2313')) = 11;
holden15.nocc(strcmp(holden15.Site,'2313')) = 11;

holden15.occtime(strcmp(holden15.Site,'2318')) = 11;
holden15.nocc(strcmp(holden15.Site,'2318')) = 8;

holden15.occtime(strcmp(holden15.Site,'2344')) = 11;
holden15.nocc(strcmp(holden15.Site,'2344')) = 9;

holden15.occtime(strcmp(holden15.Site,'2348')) = 11;
holden15.nocc(strcmp(holden15.Site,'2348')) = 8;

holden15.occtime(strcmp(holden15.Site,'2350')) = 11;
holden15.nocc(strcmp(holden15.Site,'2350')) = 7;

holden15.occtime(strcmp(holden15.Site,'2351')) = 11;
holden15.nocc(strcmp(holden15.Site,'2351')) = 7;

holden15.occtime(strcmp(holden15.Site,'2382')) = 6;
holden15.nocc(strcmp(holden15.Site,'2382')) = 7;

holden15.occtime(strcmp(holden15.Site,'2383')) = 11;
holden15.nocc(strcmp(holden15.Site,'2383')) = 8;

holden15.occtime(strcmp(holden15.Site,'2401')) = 6;
holden15.nocc(strcmp(holden15.Site,'2401')) = 6;

holden15.occtime(strcmp(holden15.Site,'2405')) = 5;
holden15.nocc(strcmp(holden15.Site,'2405')) = 3;

holden15.occtime(strcmp(holden15.Site,'3002')) = 11;
holden15.nocc(strcmp(holden15.Site,'3002')) = 9;

holden15.occtime(strcmp(holden15.Site,'3006')) = 11;
holden15.nocc(strcmp(holden15.Site,'3006')) = 8;

holden15.occtime(strcmp(holden15.Site,'3051')) = 13;
holden15.nocc(strcmp(holden15.Site,'3051')) = 24;

holden15.occtime(strcmp(holden15.Site,'3052')) = 11;
holden15.nocc(strcmp(holden15.Site,'3052')) = 8;

holden15.occtime(strcmp(holden15.Site,'3056')) = 11;
holden15.nocc(strcmp(holden15.Site,'3056')) = 8;

holden15.occtime(strcmp(holden15.Site,'3061')) = 11;
holden15.nocc(strcmp(holden15.Site,'3061')) = 8;

holden15.occtime(strcmp(holden15.Site,'3063')) = 11;
holden15.nocc(strcmp(holden15.Site,'3063')) = 8;

holden15.occtime(strcmp(holden15.Site,'3066')) = 11;
holden15.nocc(strcmp(holden15.Site,'3066')) = 8;

holden15.occtime(strcmp(holden15.Site,'3067')) = 11;
holden15.nocc(strcmp(holden15.Site,'3067')) = 8;

holden15.occtime(strcmp(holden15.Site,'3077')) = 11;
holden15.nocc(strcmp(holden15.Site,'3077')) = 8;

holden15.occtime(strcmp(holden15.Site,'3078')) = 11;
holden15.nocc(strcmp(holden15.Site,'3078')) = 9;

holden15.occtime(strcmp(holden15.Site,'3079')) = 11;
holden15.nocc(strcmp(holden15.Site,'3079')) = 8;

holden15.occtime(strcmp(holden15.Site,'3080')) = 13;
holden15.nocc(strcmp(holden15.Site,'3080')) = 24;

holden15.occtime(strcmp(holden15.Site,'3081')) = 11;
holden15.nocc(strcmp(holden15.Site,'3081')) = 8;

holden15.occtime(strcmp(holden15.Site,'3093')) = 11;
holden15.nocc(strcmp(holden15.Site,'3093')) = 8;

holden15.occtime(strcmp(holden15.Site,'3094')) = 11;
holden15.nocc(strcmp(holden15.Site,'3094')) = 8;

holden15.occtime(strcmp(holden15.Site,'6881')) = 4;
holden15.nocc(strcmp(holden15.Site,'6881')) = 6;

holden15.occtime(strcmp(holden15.Site,'A4K6')) = 4;
holden15.nocc(strcmp(holden15.Site,'A4K6')) = 6;

holden15.occtime(strcmp(holden15.Site,'A7LF')) = 11;
holden15.nocc(strcmp(holden15.Site,'A7LF')) = 8;

holden15.occtime(strcmp(holden15.Site,'A7LJ')) = 11;
holden15.nocc(strcmp(holden15.Site,'A7LJ')) = 8;

holden15.occtime(strcmp(holden15.Site,'A86L')) = 11;
holden15.nocc(strcmp(holden15.Site,'A86L')) = 8;

holden15.occtime(strcmp(holden15.Site,'A86T')) = 4;
holden15.nocc(strcmp(holden15.Site,'A86T')) = 7;

holden15.occtime(strcmp(holden15.Site,'A86W')) = 4;
holden15.nocc(strcmp(holden15.Site,'A86W')) = 6;

holden15.occtime(strcmp(holden15.Site,'A87P')) = 4;
holden15.nocc(strcmp(holden15.Site,'A87P')) = 6;

holden15.occtime(strcmp(holden15.Site,'A88G')) = 11;
holden15.nocc(strcmp(holden15.Site,'A88G')) = 8;

holden15.occtime(strcmp(holden15.Site,'A88H')) = 11;
holden15.nocc(strcmp(holden15.Site,'A88H')) = 8;

holden15.occtime(strcmp(holden15.Site,'A886')) = 4;
holden15.nocc(strcmp(holden15.Site,'A886')) = 7;

holden15.occtime(strcmp(holden15.Site,'A89V')) = 5;
holden15.nocc(strcmp(holden15.Site,'A89V')) = 4;

holden15.occtime(strcmp(holden15.Site,'AGB9')) = 4;
holden15.nocc(strcmp(holden15.Site,'AGB9')) = 6;

holden15.occtime(strcmp(holden15.Site,'A9BM')) = 10;
holden15.nocc(strcmp(holden15.Site,'A9BM')) = 8;

holden15.occtime(strcmp(holden15.Site,'AGUG')) = 6;
holden15.nocc(strcmp(holden15.Site,'AGUG')) = 5;

holden15.occtime(strcmp(holden15.Site,'AH9J')) = 4;
holden15.nocc(strcmp(holden15.Site,'AH9J')) = 6;

holden15.occtime(strcmp(holden15.Site,'B26Q')) = 11;
holden15.nocc(strcmp(holden15.Site,'B26Q')) = 8;

holden15.occtime(strcmp(holden15.Site,'B4DP')) = 6;
holden15.nocc(strcmp(holden15.Site,'B4DP')) = 6;

holden15.occtime(strcmp(holden15.Site,'B4HJ')) = 6;
holden15.nocc(strcmp(holden15.Site,'B4HJ')) = 5;

holden15.occtime(strcmp(holden15.Site,'BE2T')) = 11;
holden15.nocc(strcmp(holden15.Site,'BE2T')) = 8;

holden15.occtime(strcmp(holden15.Site,'BEVJ')) = 13;
holden15.nocc(strcmp(holden15.Site,'BEVJ')) = 24;

holden15.occtime(strcmp(holden15.Site,'CBWR')) = 11;
holden15.nocc(strcmp(holden15.Site,'CBWR')) = 8;

holden15.occtime(strcmp(holden15.Site,'CBWU')) = 11;
holden15.nocc(strcmp(holden15.Site,'CBWU')) = 10;

holden15.occtime(strcmp(holden15.Site,'CBWV')) = 11;
holden15.nocc(strcmp(holden15.Site,'CBWV')) = 11;

holden15.occtime(strcmp(holden15.Site,'CBWW')) = 11;
holden15.nocc(strcmp(holden15.Site,'CBWW')) = 8;

holden15.occtime(strcmp(holden15.Site,'EAXX')) = 4;
holden15.nocc(strcmp(holden15.Site,'EAXX')) = 6;

holden15.occtime(strcmp(holden15.Site,'EAXY')) = 4;
holden15.nocc(strcmp(holden15.Site,'EAXY')) = 6;

holden15.occtime(strcmp(holden15.Site,'E4MH')) = 4;
holden15.nocc(strcmp(holden15.Site,'E4MH')) = 7;

for dummyindex = 1:length(holden15.Lon)
    if holden15.occtime(dummyindex)==0
        sameinmidas = strcmp(holden15.Site{dummyindex},midas.names);
        if sum(sameinmidas)>0
            holden15.occtime(dummyindex) = 2011.5 - midas.starttime(sameinmidas);
            holden15.nocc(dummyindex) = holden15.occtime(dummyindex)*365.243*midas.prcofdata(sameinmidas);
        end
    end
end
%return
%%
hamling22 = readtable('./Verticals/Ian_2022/NZ_Vertical_GPS_2003-2011_GRL.txt');
hamling22.occtime = NaN*hamling22.su;
hamling22.nocc = NaN*hamling22.su;

hamling22.disttobeavan = NaN*hamling22.lon;
for staindex = 1:length(hamling22.lon)
    [hamling22.disttobeavan(staindex),minindex] = min(sqrt(((hamling22.lon(staindex) - beavan16.Longitude).*111.1.*cosd(beavan16.Latitude)).^2 + ((hamling22.lat(staindex) - beavan16.Latitude).*111.1).^2));
    if hamling22.disttobeavan(staindex)<0.07
        hamling22.occtime(staindex) = beavan16.NYEARS(minindex);
        hamling22.nocc(staindex) = beavan16.NDAYS(minindex);
    end
    if isnan(hamling22.occtime(staindex))
        [hamling22.disttomidas(staindex),minindex] = min(sqrt(((hamling22.lon(staindex) - midas.lons).*111.1.*cosd(midas.lats)).^2 + ((hamling22.lat(staindex) - midas.lats).*111.1).^2));
        if hamling22.disttomidas(staindex)<0.07
            hamling22.occtime(staindex) = 2016.5 - midas.starttime(minindex);
            hamling22.nocc(staindex) = hamling22.occtime(staindex)*365.243*midas.prcofdata(sameinmidas);
        end
    end
    if isnan(hamling22.occtime(staindex))
        [hamling22.disttoholden(staindex),minindex] = min(sqrt(((hamling22.lon(staindex) - holden15.Lon).*111.1.*cosd(holden15.Lat)).^2 + ((hamling22.lat(staindex) - holden15.Lat).*111.1).^2));
        if hamling22.disttoholden(staindex)<0.07
            hamling22.occtime(staindex) = holden15.occtime(minindex);
            hamling22.nocc(staindex) = holden15.nocc(minindex);
        end
    end
end

hamling22.nocc(isnan(hamling22.occtime)) = 2;
hamling22.occtime(isnan(hamling22.occtime)) = 2.5;

srwn = 1./sqrt(hamling22.occtime);
swn = 2*sqrt(3)*1.1./(sqrt(hamling22.nocc).*hamling22.occtime);
%return
hamling22.su = max([hamling22.su, hamling22.su/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);

for dummyindex = 1:length(hamling22.lon)
    hamling22.names{dummyindex} = ['H22_' num2str(dummyindex)];
end

%return
%%
pickle19 = readtable('./horiz/Pickle_2019/pickle19.dat');
%return
pickle19.inbounds = and(and(pickle19.Longitude>=lonbounds(1),pickle19.Longitude<=lonbounds(2)),and(pickle19.Latitude>=latbounds(1),pickle19.Latitude<=latbounds(2)));
pickle19.inbounds(or(or(pickle19.Eunc>tol,pickle19.Nunc>tol),pickle19.Hunc>tol)) = 0;
pickle19(pickle19.inbounds==0,:) = [];

for sitename = 1:length(pickle19.Longitude)
    pickle19.sitecode{sitename} = pickle19.SITE{sitename}(1:4);
end

%pickle19.Eunc_max = pickle19.Eunc;
%pickle19.Nunc_max = pickle19.Nunc;
%pickle19.Hunc_max = pickle19.Hunc;

[~,uniquesiteindices] = unique(pickle19.sitecode);
for dummyindex = 1:length(uniquesiteindices)
    siteindex = uniquesiteindices(dummyindex);
    samesiteindices = find(strcmp(pickle19.sitecode{siteindex},pickle19.sitecode));
    numvels = length(samesiteindices);
    if numvels>1
        for component = 1:3
            if component==1
                local_vels = pickle19.E(samesiteindices);
                local_unc = pickle19.Eunc(samesiteindices);
                local_weights = 1./pickle19.Eunc(samesiteindices).^2;
            elseif component==2
                local_vels = pickle19.N(samesiteindices);
                local_unc = pickle19.Nunc(samesiteindices);
                local_weights = 1./pickle19.Nunc(samesiteindices).^2;
            elseif component==3
                local_vels = pickle19.H(samesiteindices);
                local_unc = pickle19.Hunc(samesiteindices);
                local_weights = 1./pickle19.Hunc(samesiteindices).^2;
            end

            wmeanvel = sum(local_vels.*local_weights)/sum(local_weights);

            velmat = local_vels + randn(numvels,1e6).*local_unc;
            weightmat = velmat*0 + local_weights;
            velmat = velmat(:);
            weightmat = weightmat(:);

            wmeanste1 = sqrt(sum(local_weights.*(local_vels - wmeanvel).^2)/sum(local_weights))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
            wmeanste = sqrt(sum(weightmat.*(velmat - wmeanvel).^2)/sum(weightmat))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
            wmeanste2 = sqrt(1./sum(local_weights));
            if or(wmeanste1>wmeanste,wmeanste2>wmeanste)
                pause
            end

            medste = std(local_vels)/(numvels^(1/4))*1.2533;

            wmeanmedvel = wmeanvel/2 + median(local_vels)/2;
            local_ste = sqrt((wmeanste/2).^2 + (medste/2).^2);

            if component==1
                pickle19.E(samesiteindices(1)) = wmeanmedvel;
                pickle19.Eunc(samesiteindices(1)) = local_ste;
%                pickle19.Eunc_max(samesiteindices(1)) = local_stemax;
                pickle19.E(samesiteindices(2:end)) = NaN;
                pickle19.Eunc(samesiteindices(2:end)) = NaN;
%                pickle19.Eunc_max(samesiteindices(2:end)) = NaN;
            elseif component==2
                pickle19.N(samesiteindices(1)) = wmeanmedvel;
                pickle19.Nunc(samesiteindices(1)) = local_ste;
%                pickle19.Nunc_max(samesiteindices(1)) = local_stemax;
                pickle19.N(samesiteindices(2:end)) = NaN;
                pickle19.Nunc(samesiteindices(2:end)) = NaN;
%                pickle19.Nunc_max(samesiteindices(2:end)) = NaN;
            elseif component==3
                pickle19.H(samesiteindices(1)) = wmeanmedvel;
                pickle19.Hunc(samesiteindices(1)) = local_ste;
%                pickle19.Hunc_max(samesiteindices(1)) = local_stemax;
                pickle19.H(samesiteindices(2:end)) = NaN;
                pickle19.Hunc(samesiteindices(2:end)) = NaN;
%                pickle19.Hunc_max(samesiteindices(2:end)) = NaN;
            end
        end
    end
end
pickle19(isnan(pickle19.E),:) = [];
for dummyindex = 1:length(pickle19.Longitude)
    stadists_beavan16 = sqrt(((beavan16.Longitude - pickle19.Longitude(dummyindex)).*111.1.*cosd(pickle19.Latitude(dummyindex))).^2 + ((beavan16.Latitude - pickle19.Latitude(dummyindex))*111.1).^2);
    sameinbeavan16 = and(stadists_beavan16<=combdist,strcmp(pickle19.sitecode{dummyindex},beavan16.sitecode));

    stadists_houlie = sqrt(((houlie17.lon - pickle19.Longitude(dummyindex)).*111.1.*cosd(pickle19.Latitude(dummyindex))).^2 + ((houlie17.lat - pickle19.Latitude(dummyindex))*111.1).^2);
    sameinhoulie = and(stadists_houlie<=combdist,strcmp(pickle19.sitecode{dummyindex},houlie17.name));

    stadists_midas = sqrt(((midas.lons - pickle19.Longitude(dummyindex)).*111.1.*cosd(pickle19.Latitude(dummyindex))).^2 + ((midas.lats - pickle19.Latitude(dummyindex))*111.1).^2);
    sameinmidas = and(stadists_midas<=combdist,strcmp(pickle19.sitecode{dummyindex},midas.names));

    stadists_beavan12 = sqrt(((beavan12.lon - pickle19.Longitude(dummyindex)).*111.1.*cosd(pickle19.Latitude(dummyindex))).^2 + ((beavan12.lat - pickle19.Latitude(dummyindex))*111.1).^2);
    sameinbeavan12 = and(stadists_beavan12<=combdist,strcmp(pickle19.sitecode{dummyindex},beavan12.site));

    stadists_ulr5 = sqrt(((ulr5.Lon - pickle19.Longitude(dummyindex)).*111.1.*cosd(pickle19.Latitude(dummyindex))).^2 + ((ulr5.Lat - pickle19.Latitude(dummyindex))*111.1).^2);
    sameinulr5 = and(stadists_ulr5<=combdist,strcmp(pickle19.sitecode{dummyindex},ulr5.Site));

    if sum(sameinmidas)>0
        pickle19.occtime(dummyindex) = 2017.5 - midas.starttime(sameinmidas);
        pickle19.nocc(dummyindex) = pickle19.occtime(dummyindex)*365.243*midas.prcofdata(sameinmidas);
    elseif sum(sameinbeavan16)>0
        if pickle19.Latitude(dummyindex)>-37.6
            pickle19.occtime(dummyindex) = beavan16.NYEARS(sameinbeavan16) + 4;
            pickle19.nocc(dummyindex) = beavan16.NDAYS(sameinbeavan16) + 1;
        else
            pickle19.occtime(dummyindex) = beavan16.NYEARS(sameinbeavan16);
            pickle19.nocc(dummyindex) = beavan16.NDAYS(sameinbeavan16);
        end
    elseif sum(sameinhoulie)>0
        pickle19.occtime(dummyindex) = houlie17.occtime(sameinhoulie) + 2;
        pickle19.nocc(dummyindex) = pickle19.occtime(dummyindex)*365.243*houlie17.prcofdata(sameinhoulie);
    elseif sum(sameinulr5)>0
        pickle19.occtime(dummyindex) = ulr5.T_GPS(sameinulr5) + 5;
        pickle19.nocc(dummyindex) = pickle19.occtime(dummyindex)*ulr5.prcofdata(sameinulr5);
    end
end
%return
pickle19.occtime(strcmp(pickle19.sitecode,'TONG')) = (2017.5 - 2002.1328);
pickle19.nocc(strcmp(pickle19.sitecode,'TONG')) = pickle19.occtime(strcmp(pickle19.sitecode,'TONG'))*365.243;

pickle19.occtime(strcmp(pickle19.sitecode,'VANU')) = (2017.5 - 2002.6968);
pickle19.nocc(strcmp(pickle19.sitecode,'VANU')) = pickle19.occtime(strcmp(pickle19.sitecode,'VANU'))*365.243;

pickle19.occtime(strcmp(pickle19.sitecode,'TAHI')) = (2017.5 - 2002.6968);
pickle19.nocc(strcmp(pickle19.sitecode,'TAHI')) = pickle19.occtime(strcmp(pickle19.sitecode,'TAHI'))*365.243;

pickle19.occtime(strcmp(pickle19.sitecode,'BUGK')) = 4;
pickle19.nocc(strcmp(pickle19.sitecode,'BUGK')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'A6LT')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'A6LT')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'AAJA')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'AAJA')) = 3;

pickle19.occtime(strcmp(pickle19.sitecode,'B36J')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'B36J')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'B42L')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'B42L')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'B15J')) = 4;
pickle19.nocc(strcmp(pickle19.sitecode,'B15J')) = 2;

pickle19.occtime(strcmp(pickle19.sitecode,'VTHM')) = 6;
pickle19.nocc(strcmp(pickle19.sitecode,'VTHM')) = 6*365.243;

pickle19.occtime(strcmp(pickle19.sitecode,'VWKU')) = 6;
pickle19.nocc(strcmp(pickle19.sitecode,'VWKU')) = 6*365.243;

pickle19.occtime(strcmp(pickle19.sitecode,'AAJ1')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'AAJ1')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'AAJ4')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'AAJ4')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'AAF6')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'AAF6')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'AAJ6')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'AAJ6')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'AAJ5')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'AAJ5')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'AAV5')) = 22;
%pickle19.nocc(strcmp(pickle19.sitecode,'AAV5')) = 2;

pickle19.occtime(strcmp(pickle19.sitecode,'A9PH')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'A9PH')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'B0PD')) = 4;
pickle19.nocc(strcmp(pickle19.sitecode,'B0PD')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'1320')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'1320')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'DFFL')) = 4;
pickle19.nocc(strcmp(pickle19.sitecode,'DFFL')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'1349')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'1349')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'A9FW')) = 4;
pickle19.nocc(strcmp(pickle19.sitecode,'A9FW')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'AAJ0')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'AAJ0')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'1331')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'1331')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'VWHK')) = 6.8;
pickle19.nocc(strcmp(pickle19.sitecode,'VWHK')) = 6.8*365.243;

pickle19.occtime(strcmp(pickle19.sitecode,'VPRK')) = 6.8;
pickle19.nocc(strcmp(pickle19.sitecode,'VPRK')) = 6.8*365.243;

pickle19.occtime(strcmp(pickle19.sitecode,'A9MK')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'A9MK')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'A6Q1')) = 4;
pickle19.nocc(strcmp(pickle19.sitecode,'A6Q1')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'AAEJ')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'AAEJ')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'AAEB')) = 22;
pickle19.nocc(strcmp(pickle19.sitecode,'AAEB')) = 4;

pickle19.occtime(strcmp(pickle19.sitecode,'A42F')) = 4;
pickle19.nocc(strcmp(pickle19.sitecode,'A42F')) = 2;

pickle19.occtime(strcmp(pickle19.sitecode,'A49Q')) = 4;
pickle19.nocc(strcmp(pickle19.sitecode,'A49Q')) = 2;

pickle19.occtime(strcmp(pickle19.sitecode,'A05E')) = 4;
pickle19.nocc(strcmp(pickle19.sitecode,'A05E')) = 2;

pickle19.occtime(strcmp(pickle19.sitecode,'1344')) = 22;

pickle19.occtime(strcmp(pickle19.sitecode,'A93W')) = 2;
pickle19.occtime(strcmp(pickle19.sitecode,'A8TX')) = 2;
pickle19.occtime(strcmp(pickle19.sitecode,'ER3B')) = 2;

pickle19.occtime(pickle19.occtime==0) = 2;
pickle19.nocc(and(pickle19.occtime==2,pickle19.nocc==0)) = 3;

srwn = 1./sqrt(pickle19.occtime);
swn = 2.*sqrt(3).*1.1./(sqrt(pickle19.nocc).*pickle19.occtime);
%return
pickle19.Eunc = max([pickle19.Eunc,pickle19.Eunc/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
pickle19.Nunc = max([pickle19.Nunc,pickle19.Nunc/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
pickle19.Hunc = max([pickle19.Hunc,pickle19.Hunc/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);

[x,y,z] = wgs2xyz(pickle19.Longitude,pickle19.Latitude,zeros(length(pickle19.Latitude),1));
OM_pickle = [0.411619    0.325891  0.343019]';
[Vxyz,Venu] = rotate(OM_pickle,[0*OM_pickle 0*OM_pickle 0*OM_pickle],[x y z]);
pickle19.E_IT08 = pickle19.E + Venu(:,1);
pickle19.N_IT08 = pickle19.N + Venu(:,2);

[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
pickle19.veast_b16 = pickle19.E_IT08 + Venu(:,1);
pickle19.vnorth_b16 = pickle19.N_IT08 + Venu(:,2);
%
%%
denys16_pre09 = readtable('./horiz/Denys_2016/denys16_pre2009.dat');

[x,y,z] = wgs2xyz(denys16_pre09.Lon,denys16_pre09.Lat,zeros(length(denys16_pre09.Lon),1));
[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
denys16_pre09.veast_b16 = denys16_pre09.VE + Venu(:,1);
denys16_pre09.vnorth_b16 = denys16_pre09.VN + Venu(:,2);
%return
%%
denys16_post09 = readtable('./horiz/Denys_2016/denys16_post2009.dat');
denys16_post09.Lon = denys16_post09.Var1;
denys16_post09.Lat = denys16_post09.Var2;
denys16_post09.VE = denys16_post09.Var3;
denys16_post09.VN = denys16_post09.Var4;
denys16_post09.SE = denys16_post09.Var5;
denys16_post09.SN = denys16_post09.Var6;
denys16_post09.ID = denys16_post09.Var7;

[x,y,z] = wgs2xyz(denys16_post09.Lon,denys16_post09.Lat,zeros(length(denys16_post09.Lon),1));
[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
denys16_post09.veast_b16 = denys16_post09.VE + Venu(:,1);
denys16_post09.vnorth_b16 = denys16_post09.VN + Venu(:,2);

%%
denys16_all = readtable('./horiz/Denys_2016/denys16_all.dat');

[x,y,z] = wgs2xyz(denys16_all.Lon,denys16_all.Lat,zeros(length(denys16_all.Lon),1));
[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
denys16_all.veast_b16 = denys16_all.VE + Venu(:,1);
denys16_all.vnorth_b16 = denys16_all.VN + Venu(:,2);

%%
page18_pre09 = readtable('./horiz/Page_2018/page18_pre09.dat');
[x,y,z] = wgs2xyz(page18_pre09.Longitude,page18_pre09.Latitude,zeros(length(page18_pre09.Latitude),1));

page18_pre09.occtime = 0*page18_pre09.Vn_ITRF14 + 4;
page18_pre09.nocc = 0*page18_pre09.Vn_ITRF14 + 5;

srwn = 1./sqrt(page18_pre09.occtime);
swn = 2.*sqrt(3).*1.1./(sqrt(page18_pre09.nocc).*page18_pre09.occtime);
page18_pre09.Se = max([page18_pre09.Se,page18_pre09.Se/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
page18_pre09.Sn = max([page18_pre09.Sn,page18_pre09.Sn/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);

O = [page18_pre09.Latitude,page18_pre09.Longitude,zeros(length(page18_pre09.Longitude),1)];
V = [page18_pre09.Vn_ITRF14,page18_pre09.Ve_ITRF14,0*page18_pre09.Ve_ITRF14]*1e-3;
SV = [page18_pre09.Sn,page18_pre09.Se,0*page18_pre09.Se]*1e-3;

[vxyz,cxyz] = neu2xyz(O,V,SV,0*SV);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2014',[],'ITRF2008',[],[]);
[itrfout] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
page18_pre09.vnorth_itrf08 = itrfout(:,1)*1e3;
page18_pre09.veast_itrf08 = itrfout(:,2)*1e3;

[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
page18_pre09.veast_b16 = page18_pre09.veast_itrf08 + Venu(:,1);
page18_pre09.vnorth_b16 = page18_pre09.vnorth_itrf08 + Venu(:,2);

%%
page18_post09 = readtable('./horiz/Page_2018/page18_post09.dat');
[x,y,z] = wgs2xyz(page18_post09.Longitude,page18_post09.Latitude,zeros(length(page18_post09.Latitude),1));

page18_post09.occtime = 0*page18_post09.Vn_ITRF14 + 4;
page18_post09.nocc = 0*page18_post09.Vn_ITRF14 + 5;

srwn = 1./sqrt(page18_post09.occtime);
swn = 2.*sqrt(3).*1.1./(sqrt(page18_post09.nocc).*page18_post09.occtime);
page18_post09.Se = max([page18_post09.Se,page18_post09.Se/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
page18_post09.Sn = max([page18_post09.Sn,page18_post09.Sn/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);

O = [page18_post09.Latitude,page18_post09.Longitude,zeros(length(page18_post09.Longitude),1)];
V = [page18_post09.Vn_ITRF14,page18_post09.Ve_ITRF14,0*page18_post09.Ve_ITRF14]*1e-3;
SV = [page18_post09.Sn,page18_post09.Se,0*page18_post09.Se]*1e-3;

[vxyz,cxyz] = neu2xyz(O,V,SV,0*SV);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2014',[],'ITRF2008',[],[]);
[itrfout] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
page18_post09.vnorth_itrf08 = itrfout(:,1)*1e3;
page18_post09.veast_itrf08 = itrfout(:,2)*1e3;

[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
page18_post09.veast_b16 = page18_post09.veast_itrf08 + Venu(:,1);
page18_post09.vnorth_b16 = page18_post09.vnorth_itrf08 + Venu(:,2);

%return
%%
denys20 = readtable('./Verticals/Denys_2020/denys20.dat');
%disp('denys20')
%return
%%
beavan10 = readtable('./Verticals/Beavan_2010/beavan10.dat');
beavan10.lon = beavan10.lon + 360;
%beavan10.su = sqrt(beavan10.su.^2 + 0.5^2);

O = [beavan10.lat,beavan10.lon,zeros(length(beavan10.lon),1)];
V = [0*beavan10.vu,0*beavan10.vu,beavan10.vu]*1e-3;
SV = [0*beavan10.vu,0*beavan10.vu,beavan10.su]*1e-3;

[x,y,z] = wgs2xyz(beavan10.lon,beavan10.lat,zeros(length(beavan10.lat),1));
[vxyz,~] = neu2xyz(O,V,SV,0*SV);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2000',[],'ITRF2008',[],[]);
[itrfout] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
beavan10.vup_itrf08 = itrfout(:,3)*1e3;
beavan10.vup_itrf08 = beavan10.vup_itrf08 - beavan10.vup_itrf08(strcmp(beavan10.name,'OUSD'));
%return
%%
ffile = fopen(['horiz/kreemer14/GPS_ITRF08.gmt']);
gsrminput = textscan(ffile,'%f %f %f %f %f %f %f %s %s\n');
fclose(ffile);

kreemer14.lons = gsrminput{1};
kreemer14.lats = gsrminput{2};
kreemer14.veast = gsrminput{3};
kreemer14.vnorth = gsrminput{4};
kreemer14.seast = gsrminput{5};
kreemer14.snorth = gsrminput{6};
kreemer14.names = gsrminput{8};
kreemer14.study = gsrminput{9};
kreemer14 = struct2table(kreemer14);

kreemer14.inbounds = and(and(kreemer14.lons>=lonbounds(1),kreemer14.lons<=lonbounds(2)),and(kreemer14.lats>=latbounds(1),kreemer14.lats<=latbounds(2)));
kreemer14(kreemer14.inbounds==0,:) = [];

UNR = kreemer14(strcmp(kreemer14.study,'UNR'),:);

[x,y,z] = wgs2xyz(UNR.lons,UNR.lats,zeros(length(UNR.lats),1));
[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
UNR.veast_b16 = UNR.veast + Venu(:,1);
UNR.vnorth_b16 = UNR.vnorth + Venu(:,2);

%%
ffile = fopen(['horiz/measures/measures_2025.dat']);
input = textscan(ffile,'%s %f %f %*f %f %f %f %f %f %f\n','CommentStyle','Site');
fclose(ffile);

measures.names = upper(input{1});
measures.lats = input{3};
measures.lons = input{2};
measures.vnorth = input{4};
measures.veast = input{5};
measures.vup = input{6};
measures.snorth = input{7};
measures.seast = input{8};
measures.su = input{9};
measures = struct2table(measures);

measures.inbounds = and(and(measures.lons>=lonbounds(1),measures.lons<=lonbounds(2)),and(measures.lats>=latbounds(1),measures.lats<=latbounds(2)));
measures(measures.inbounds==0,:) = [];

measures.occtime = 0*measures.vup + 2.5;
measures.nocc = measures.occtime*365.243;
for measuresindex = 1:length(measures.names)
    sameinmidas = strcmp(midas.names,measures.names(measuresindex));
    if sum(sameinmidas)>0
        measures.occtime(measuresindex) = midas.occtime(sameinmidas);
        measures.nocc(measuresindex) = midas.nocc(sameinmidas);
    end
end
%measures.occtime(measures.occtime==0) = NaN;

measures.occtime(strcmp(measures.names,'KUTA')) = (datenum(2025,02,27) - (datenum(2000,1,1) + 3749)/365.243);
measures.occtime(strcmp(measures.names,'TGWH')) = (datenum(2025,02,27) - (datenum(2000,1,1) + 2594)/365.243);
measures.occtime(strcmp(measures.names,'TONG')) = (datenum(2025,02,27) - 2002.1328)/365.243;
measures.occtime(strcmp(measures.names,'VANU')) = (datenum(2025,02,27) - 2002.6968)/365.243;
measures.nocc(strcmp(measures.names,'KUTA')) = measures.occtime(strcmp(measures.names,'KUTA'))*365.243;
measures.nocc(strcmp(measures.names,'TGWH')) = measures.occtime(strcmp(measures.names,'TGWH'))*365.243;
measures.nocc(strcmp(measures.names,'TONG')) = measures.occtime(strcmp(measures.names,'TONG'))*365.243;
measures.nocc(strcmp(measures.names,'VANU')) = measures.occtime(strcmp(measures.names,'VANU'))*365.243;
%measures.occtime = max(measures.occtime,2.5);
%return
srwn = 1./sqrt(measures.occtime);
swn = 2.*sqrt(3).*1.1./(sqrt(measures.nocc).*measures.occtime);
measures.su = max([measures.su, measures.su/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
measures.seast = max([measures.seast, measures.seast/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
measures.snorth = max([measures.snorth, measures.snorth/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);

O = [measures.lats,measures.lons,zeros(length(measures.lons),1)];
V = [measures.vnorth,measures.veast,measures.vup]*1e-3;
SV = [measures.snorth,measures.seast,measures.su]*1e-3;

[x,y,z] = wgs2xyz(measures.lons,measures.lats,zeros(length(measures.lats),1));
[vxyz,cxyz] = neu2xyz(O,V,SV,0*SV);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2014',[],'ITRF2008',[],[]);
[itrfout] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
measures.vnorth_itrf08 = itrfout(:,1)*1e3;
measures.veast_itrf08 = itrfout(:,2)*1e3;
measures.vup_itrf08 = itrfout(:,3)*1e3;

[x,y,z] = wgs2xyz(measures.lons,measures.lats,zeros(length(measures.lats),1));
[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
measures.veast_b16 = measures.veast_itrf08 + Venu(:,1);
measures.vnorth_b16 = measures.vnorth_itrf08 + Venu(:,2);

%%
input = readtable(['~/My Drive (john.c.rollins@gmail.com)/AHB/Material/gps/itrf14vlbi/ITRF2014_VLBI.SSC.txt']);
itrf14vlbi.x = input.Var5(5:2:end);
itrf14vlbi.y = input.Var6(5:2:end);
itrf14vlbi.z = input.Var7(5:2:end);
itrf14vlbi.vx = 0*itrf14vlbi.x;
itrf14vlbi.vy = 0*itrf14vlbi.y;
itrf14vlbi.vz = 0*itrf14vlbi.z;
itrf14vlbi.sx = input.Var5(6:2:end);
itrf14vlbi.sy = input.Var6(6:2:end);
itrf14vlbi.sz = input.Var7(6:2:end);

for index = 1:length(itrf14vlbi.x)
    itrf14vlbi.vx(index) = str2num(input.Var2{4+index*2});
    itrf14vlbi.vy(index) = str2num(input.Var3{4+index*2});
    itrf14vlbi.vz(index) = input.Var4(4+index*2);
end
itrf14vlbi.names = cell(length(itrf14vlbi.x),1);

for index = 1:length(itrf14vlbi.x)
    itrf14vlbi.names{index} = ['VLBI' num2str(index)];
end
%return

itrfll = xyz2wgs([0*itrf14vlbi.x itrf14vlbi.x itrf14vlbi.y itrf14vlbi.z]);
itrf14vlbi.lons = itrfll(:,2);
itrf14vlbi.lats = itrfll(:,3);
itrf14vlbi.lons(itrf14vlbi.lons<0) = itrf14vlbi.lons(itrf14vlbi.lons<0) + 360;

[itrfout,err] = xyz2neu([itrf14vlbi.x,itrf14vlbi.y,itrf14vlbi.z],[itrf14vlbi.vx,itrf14vlbi.vy,itrf14vlbi.vz],[itrf14vlbi.sx,itrf14vlbi.sy,itrf14vlbi.sz],0*[itrf14vlbi.sx,itrf14vlbi.sy,itrf14vlbi.sz]);
itrf14vlbi.vnorth_itrf14 = itrfout(:,1)*1e3;
itrf14vlbi.veast_itrf14 = itrfout(:,2)*1e3;
itrf14vlbi.vup_itrf14 = itrfout(:,3)*1e3;
itrf14vlbi.snorth = sqrt(err(:,1))*1e3;
itrf14vlbi.seast = sqrt(err(:,4))*1e3;
itrf14vlbi.sup = sqrt(err(:,6))*1e3;

% itrf14vlbi.occtime = 0*itrf14vlbi.lons + 10;
% itrf14vlbi.nocc = 0*itrf14vlbi.lons + 10;
itrf14vlbi = struct2table(itrf14vlbi);
%return
itrf14vlbi.inbounds = and(and(itrf14vlbi.lons>=lonbounds(1),itrf14vlbi.lons<=lonbounds(2)),and(itrf14vlbi.lats>=latbounds(1),itrf14vlbi.lats<=latbounds(2)));

itrf14vlbi(itrf14vlbi.inbounds==0,:) = [];

O = [itrf14vlbi.lats,itrf14vlbi.lons,zeros(length(itrf14vlbi.lons),1)];
V = [itrf14vlbi.vnorth_itrf14,itrf14vlbi.veast_itrf14,itrf14vlbi.vup_itrf14]*1e-3;
SV = [itrf14vlbi.snorth,itrf14vlbi.seast,itrf14vlbi.sup]*1e-3;

[x,y,z] = wgs2xyz(itrf14vlbi.lons,itrf14vlbi.lats,zeros(length(itrf14vlbi.lats),1));
[vxyz,cxyz] = neu2xyz(O,V,SV,0*SV);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2014',[],'ITRF2008',[],[]);
[itrfout] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
itrf14vlbi.vnorth_itrf08 = itrfout(:,1)*1e3;
itrf14vlbi.veast_itrf08 = itrfout(:,2)*1e3;
itrf14vlbi.vup_itrf08 = itrfout(:,3)*1e3;

[x,y,z] = wgs2xyz(itrf14vlbi.lons,itrf14vlbi.lats,zeros(length(itrf14vlbi.lats),1));
[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
itrf14vlbi.veast_b16 = itrf14vlbi.veast_itrf08 + Venu(:,1);
itrf14vlbi.vnorth_b16 = itrf14vlbi.vnorth_itrf08 + Venu(:,2);

%%
input = readtable(['~/My Drive (john.c.rollins@gmail.com)/AHB/Material/GPS/vardic21/vardic21.dat']);
vardic21.names = input.Var1;
vardic21.lats = input.Var2;
vardic21.lons = input.Var3;
vardic21.lons(vardic21.lons<0) = vardic21.lons(vardic21.lons<0) + 360;
vardic21.vnorth = input.Var4;
vardic21.snorth = input.Var5;
vardic21.veast = input.Var6;
vardic21.seast = input.Var7;
vardic21.vup = input.Var8;
vardic21.sup = input.Var9;
%vardic21.occtime = 0*vardic21.lons + 9;
%vardic21.nocc = vardic21.occtime*365.243;
vardic21 = struct2table(vardic21);

vardic21.inbounds = and(and(vardic21.lons>=lonbounds(1),vardic21.lons<=lonbounds(2)),and(vardic21.lats>=latbounds(1),vardic21.lats<=latbounds(2)));
vardic21(vardic21.inbounds==0,:) = [];

% srwn = 1./sqrt(vardic21.occtime);
% swn = 2.*sqrt(3).*1.1./(sqrt(vardic21.nocc).*vardic21.occtime);
% vardic21.sup = max([vardic21.sup, vardic21.sup/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
% vardic21.seast = max([vardic21.seast, vardic21.seast/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
% vardic21.snorth = max([vardic21.snorth, vardic21.snorth/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);

O = [vardic21.lats,vardic21.lons,zeros(length(vardic21.lons),1)];
V = [vardic21.vnorth,vardic21.veast,vardic21.vup]*1e-3;
SV = [vardic21.snorth,vardic21.seast,vardic21.sup]*1e-3;

[x,y,z] = wgs2xyz(vardic21.lons,vardic21.lats,zeros(length(vardic21.lats),1));
[vxyz,cxyz] = neu2xyz(O,V,SV,0*SV);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2014',[],'ITRF2008',[],[]);
[velout] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
vardic21.vnorth_itrf08 = velout(:,1)*1e3;
vardic21.veast_itrf08 = velout(:,2)*1e3;
vardic21.vup_itrf08 = velout(:,3)*1e3;

[x,y,z] = wgs2xyz(vardic21.lons,vardic21.lats,zeros(length(vardic21.lats),1));
[Vxyz,Venu] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
vardic21.veast_b16 = vardic21.veast_itrf08 + Venu(:,1);
vardic21.vnorth_b16 = vardic21.vnorth_itrf08 + Venu(:,2);

%return
%%
lamb13 = readtable('./horiz/Lamb_2013/lamb13_sed4.dat');
[x,y,z] = wgs2xyz(lamb13.lon,lamb13.lat,zeros(length(lamb13.lat),1));
E = [37.54	32.76  0.621];
OM_b02 = eul2rot(E)';
[~,Venu_b02] = rotate(OM_b02,[0*OM_b02 0*OM_b02 0*OM_b02],[x y z]);
lamb13.ve = lamb13.ve + Venu_b02(:,1);
lamb13.vn = lamb13.vn + Venu_b02(:,2);

lamb13.nocc = lamb13.endday - lamb13.startday;
lamb13.startdate = (datenum(2000,1,3) + lamb13.startday)/365.243;
lamb13.occtime = lamb13.nocc/365.243;
%return
srwn = 1./sqrt(lamb13.occtime);
swn = 2.*sqrt(3).*1.1./(sqrt(lamb13.nocc).*lamb13.occtime);
% lamb13.su = sqrt(lamb13.su.^2 + srwn.^2 + swn.^2);
% lamb13.se = sqrt(lamb13.se.^2 + srwn.^2 + swn.^2);
% lamb13.sn = sqrt(lamb13.sn.^2 + srwn.^2 + swn.^2);
lamb13.su = max([lamb13.su, lamb13.su/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
lamb13.se = max([lamb13.se, lamb13.se/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);
lamb13.sn = max([lamb13.sn, lamb13.sn/2 + sqrt(srwn.^2 + swn.^2)/2],[],2);

O = [lamb13.lat,lamb13.lon,zeros(length(lamb13.lon),1)];
V = [lamb13.vn,lamb13.ve,lamb13.vu]*1e-3;
SV = [lamb13.sn,lamb13.se,lamb13.su]*1e-3;

[x,y,z] = wgs2xyz(lamb13.lon,lamb13.lat,zeros(length(lamb13.lat),1));
[vxyz,~] = neu2xyz(O,V,SV,0*SV);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2000',[],'ITRF2008',[],[]);
[itrfout] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
lamb13.vn_IT08 = itrfout(:,1)*1e3;
lamb13.ve_IT08 = itrfout(:,2)*1e3;
lamb13.vu_IT08 = itrfout(:,3)*1e3;
% 
% [~,Venu_b16] = rotate(OM_b16,[0*OM_b16 0*OM_b16 0*OM_b16],[x y z]);
% lamb13.ve_b16 = lamb13.ve_IT08 + Venu_b16(:,1);
% lamb13.vn_b16 = lamb13.vn_IT08 + Venu_b16(:,2);
%lamb13.VE_AUS = lamb13.ve;
%lamb13.VN_AUS = lamb13.vn;

%return
%%
beavan16.outlier = or(beavan16.x_E>tol,beavan16.x_N>tol);
beavan16.outlier(strcmp(beavan16.sitecode,'AAV5')) = 1;
beavan16.outlier(strcmp(beavan16.sitecode,'1367')) = 1;
beavan16.outlier(strcmp(beavan16.sitecode,'DD1Y')) = 1;
beavan16.outlier(strcmp(beavan16.sitecode,'2404')) = 1;
beavan16.outlier(strcmp(beavan16.sitecode,'A934')) = 1;
beavan16.outlier(strcmp(beavan16.sitecode,'AB5P')) = 1;
beavan16.outlier(strcmp(beavan16.sitecode,'A8YC')) = 1;
%beavan16.outlier(strcmp(beavan16.sitecode,'2011')) = 1;
%beavan16.outlier(strcmp(beavan16.sitecode,'2353')) = 1;
beavan16.outlier(strcmp(beavan16.sitecode,'A1YC')) = 1;
beavan16.outlier(strcmp(beavan16.sitecode,'AH4Q')) = 1;
%beavan16.outlier(strcmp(beavan16.sitecode,'WARK')) = 1;
%beavan16.outlier(strcmp(beavan16.sitecode,'KTIA')) = 1;
%beavan16.outlier(strcmp(beavan16.sitecode,'AAJU')) = 1;

beavan12.outlier = beavan12.su>tol;

houlie17.outlier = houlie17.su>tol;
houlie17.outlier(strcmp(houlie17.name,'UTKU')) = 1;
houlie17.outlier(strcmp(houlie17.name,'UTK1')) = 1;
houlie17.outlier(strcmp(houlie17.name,'UTK2')) = 1;
houlie17.outlier(strcmp(houlie17.name,'UTK3')) = 1;
houlie17.outlier(strcmp(houlie17.name,'UTK4')) = 1;
houlie17.outlier(strcmp(houlie17.name,'BHST')) = 1;
houlie17.outlier(strcmp(houlie17.name,'VGMO')) = 1;
%houlie17.neardusky = and(and(houlie17.lon>165,houlie17.lon<168.5),and(houlie17.lat>-48,houlie17.lat<-43.5));
%houlie17.outlier(houlie17.neardusky==1) = 1;

midas.outlier = midas.occtime<minocctime;
midas.outlier(or(or(midas.seast>tol,midas.snorth>tol),midas.su>tol)) = 1;

midas.nearkaikourachch = and(and(midas.lons<174.7,midas.lats<-40.5),and(midas.lons>172.5,midas.lats>-43.6));
%midas.neardusky = and(and(midas.lons>165,midas.lons<168.5),and(midas.lats>-48,midas.lats<-43.5));
midas.outlier(midas.nearkaikourachch==1) = 1;
%midas.outlier(midas.neardusky==1) = 1;
%midas.nearchch = and(and(midas.lons>172.5,midas.lats<-43.5),midas.lons<173);
%midas.outlier(midas.nearchch==1) = 1;
midas.manawatu = and(and(and(midas.lons>170,midas.lons<176.25),midas.lats>-40.15),midas.lats<-39.6);
midas.outlier(midas.manawatu) = 1;
midas.outlier(strcmp(midas.names,'PNUI')) = 0;

midas.outlier(strcmp(midas.names,'MAKO')) = 1;
midas.outlier(strcmp(midas.names,'PARI')) = 1;
midas.outlier(strcmp(midas.names,'GUNR')) = 1;
midas.outlier(strcmp(midas.names,'UTKU')) = 1;
midas.outlier(strcmp(midas.names,'UTK1')) = 1;
midas.outlier(strcmp(midas.names,'UTK2')) = 1;
midas.outlier(strcmp(midas.names,'UTK3')) = 1;
midas.outlier(strcmp(midas.names,'UTK4')) = 1;
%midas.outlier(strcmp(midas.names,'RGKW')) = 1;
midas.outlier(strcmp(midas.names,'RGWI')) = 1;
midas.outlier(strcmp(midas.names,'RGWC')) = 1;
%midas.outlier(strcmp(midas.names,'BHST')) = 1;
midas.outlier(strcmp(midas.names,'VGMO')) = 1;
midas.outlier(strcmp(midas.names,'Z406')) = 1;
midas.outlier(strcmp(midas.names,'PTVL')) = 1;
midas.outlier(strcmp(midas.names,'4MCK')) = 1;
midas.outlier(strcmp(midas.names,'QWLP')) = 1;

ulr5.outlier = ulr5.S_GPS>tol;

holden15.outlier = or(holden15.sVUP>tol,or(holden15.sVE>tol,holden15.sVN>tol));
holden15.outlier(strcmp(holden15.Site,'WHNG')) = 1;
holden15.outlier(strcmp(holden15.Site,'KTIA')) = 1;
%holden15.outlier(strcmp(holden15.Site,'RGOP')) = 1;
%holden15.outlier(strcmp(holden15.Site,'RGMT')==1) = 1;
holden15.outlier(strcmp(holden15.Site,'2401')) = 1;
%holden15.outlier(strcmp(holden15.Site,'RGHL')) = 1;
%holden15.outlier(strcmp(holden15.Site,'RGKW')==1) = 1;
holden15.outlier(strcmp(holden15.Site,'HAMT')) = 1;
holden15.outlier(strcmp(holden15.Site,'MAHO')) = 1;

hamling22.outlier = hamling22.su>tol;
hamling22.outlier(and(hamling22.lat==-39.4892,hamling22.lon==176.0632)) = 1;

pickle19.outlier = or(or(pickle19.Eunc>tol,pickle19.Nunc>tol),pickle19.Hunc>tol);
pickle19.outlier(strcmp(pickle19.sitecode,'A6BK')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'A6X8')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'A8UK')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'A93Q')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'A9FJ')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'A9FL')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'B0T9')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'B1P3')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'EJTY')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'EC0T')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'A2J7')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'B0PC')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'CAJP')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'B0TB')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSWH')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSTN')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSSB')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSAA')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSAL')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSKK')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSDR')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSBL')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSTK')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSTH')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSCA')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSHT')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSWI')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSUT')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSHW')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'GSMH')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'1362')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'1354')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'1350')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'1373')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'A2L2')) = 1;
%
pickle19.outlier(strcmp(pickle19.sitecode,'GISB')) = 1;
%pickle19.outlier(strcmp(pickle19.sitecode,'HAST')) = 1;
%pickle19.outlier(strcmp(pickle19.sitecode,'CKID')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'BHST')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'VGMO')) = 1;
%pickle19.outlier(strcmp(pickle19.sitecode,'METH')) = 1;
%pickle19.outlier(strcmp(pickle19.sitecode,'RGKW')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'RGWI')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'RGWC')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'1367')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'HOKI')) = 1;
%pickle19.outlier(strcmp(pickle19.sitecode,'A1YC')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'2404')) = 1;
pickle19.outlier(strcmp(pickle19.sitecode,'B2X8')) = 1;
%pickle19.outlier(strcmp(pickle19.sitecode,'PYGR')) = 1;

%pickle19.nearkaikoura = and(and(pickle19.Longitude<174.7,pickle19.Latitude<=-40.5),and(pickle19.Longitude>168.5,pickle19.Latitude>-43.5));
%pickle19.outlier(pickle19.nearkaikoura==1) = 1;

denys16_pre09.outlier = or(denys16_pre09.SE>tol,denys16_pre09.SN>tol);
denys16_pre09.outlier(strcmp(denys16_pre09.ID,'DME2')) = 1;
denys16_pre09.outlier(strcmp(denys16_pre09.ID,'ECFG')) = 1;

denys16_post09.outlier = or(denys16_post09.SE>tol,denys16_post09.SN>tol);
denys16_post09.outlier(strcmp(denys16_post09.ID,'DME2')) = 1;

denys16_all.outlier = or(denys16_all.SE>tol,denys16_all.SN>tol);
denys16_all.outlier(strcmp(denys16_all.ID,'DME2')) = 1;

page18_pre09.outlier = or(page18_pre09.Se>tol,page18_pre09.Sn>tol);
page18_pre09.outlier(strcmp(page18_pre09.site,'ECFG')) = 1;
page18_pre09.outlier(strcmp(page18_pre09.site,'ECFK')) = 1;

page18_post09.outlier = or(page18_post09.Se>tol,page18_post09.Sn>tol);
page18_post09.outlier(strcmp(page18_post09.site,'AADY')) = 1;

denys20.outlier = denys20.su>tol;

beavan10.outlier = beavan10.su>tol;

UNR.outlier = or(or(UNR.seast>tol,UNR.snorth>tol),isnan(UNR.seast));
% UNR.outlier(and(UNR.lons>175.55,and(UNR.lats<=-38.2,UNR.lats>=-39.1))) = 1;
UNR.outlier(and(UNR.lons>174.5,and(UNR.lats<=-39,UNR.lats>=-40.2))) = 1;
% UNR.outlier(and(and(UNR.lons<177.6,UNR.lats>-38.65),and(UNR.seast<=tol,UNR.snorth<=tol))) = 0;
UNR.outlier(strcmp(UNR.names,'CLSK')) = 1;
UNR.outlier(strcmp(UNR.names,'PYGR')) = 1;
%UNR.outlier(strcmp(UNR.names,'AUKT')) = 1;
UNR.outlier(strcmp(UNR.names,'ANAU')) = 1;
UNR.outlier(strcmp(UNR.names,'MAKO')) = 1;
UNR.outlier(strcmp(UNR.names,'AHTI')) = 1;
UNR.outlier(strcmp(UNR.names,'GISB')) = 1;
UNR.outlier(strcmp(UNR.names,'HANA')) = 1;
UNR.outlier(strcmp(UNR.names,'PRTU')) = 1;
UNR.outlier(strcmp(UNR.names,'PRTU')) = 1;
%UNR.outlier(strcmp(UNR.names,'KOKO')) = 1;
%UNR.outlier(strcmp(UNR.names,'VAHU')) = 1;

measures.outlier = or(or(measures.seast>tol,measures.snorth>tol),measures.su>tol);
%measures.tooshort = measures.occtime<minocctime;
measures.outlier(strcmp(measures.names,'RGWI')) = 1;
measures.outlier(strcmp(measures.names,'RGWC')) = 1;
%measures.outlier(strcmp(measures.names,'RGKW')) = 1;
%measures.outlier(strcmp(measures.names,'BHST')) = 1;
measures.outlier(strcmp(measures.names,'MRBL')) = 1;
measures.outlier(strcmp(measures.names,'HANM')) = 1;
measures.outlier(strcmp(measures.names,'TENN')) = 1;
measures.outlier(strcmp(measures.names,'WRAU')) = 1;
measures.outlier(strcmp(measures.names,'CLRR')) = 1;
measures.outlier(strcmp(measures.names,'LOOK')) = 1;
measures.outlier(strcmp(measures.names,'GLOK')) = 1;
measures.outlier(strcmp(measures.names,'VGMO')) = 1;

itrf14vlbi.outlier = or(or(itrf14vlbi.seast>tol,itrf14vlbi.snorth>tol),itrf14vlbi.sup>tol);

vardic21.outlier = or(or(vardic21.seast>tol,vardic21.snorth>tol),vardic21.sup>tol);

lamb13.outlier = or(or(lamb13.occtime<minocctime,lamb13.se>tol),or(lamb13.sn>tol,lamb13.su>tol));
lamb13.outlier(strcmp(lamb13.name,'BHST')) = 1;
lamb13.outlier(strcmp(lamb13.name,'VGMO')) = 1;
lamb13.manawatu = and(and(lamb13.lon<176.55,lamb13.lat>-40.2),lamb13.lat<-39.6);
lamb13.outlier(lamb13.manawatu) = 1;

%%
% denys16.lons = [denys16_pre09.Lon; denys16_post09.Lon; denys16_all.Lon];
% denys16.lats = [denys16_pre09.Lat; denys16_post09.Lat; denys16_all.Lat];
% denys16.ID = [denys16_pre09.ID; denys16_post09.ID; denys16_all.ID];
% denys16.veast_b16 = [denys16_pre09.veast_b16; denys16_post09.veast_b16; denys16_all.veast_b16];
% denys16.vnorth_b16 = [denys16_pre09.vnorth_b16; denys16_post09.vnorth_b16; denys16_all.vnorth_b16];
% denys16.seast = [denys16_pre09.SE; denys16_post09.SE; denys16_all.SE];
% denys16.snorth = [denys16_pre09.SN; denys16_post09.SN; denys16_all.SN];
% denys16.pagepre09 = [denys16_pre09.Lon*0 + 1; denys16_post09.Lon*0; denys16_all.Lon*0 + 1];
% % denys16.nocc = [denys16_pre09.nocc; denys16_post09.nocc; denys16_all.nocc];
% denys16.outlier = [denys16_pre09.outlier; denys16_post09.outlier; denys16_all.outlier];

denys16.lons = [denys16_pre09.Lon; denys16_all.Lon];
denys16.lats = [denys16_pre09.Lat; denys16_all.Lat];
denys16.ID = [denys16_pre09.ID; denys16_all.ID];
denys16.veast_b16 = [denys16_pre09.veast_b16; denys16_all.veast_b16];
denys16.vnorth_b16 = [denys16_pre09.vnorth_b16; denys16_all.vnorth_b16];
denys16.seast = [denys16_pre09.SE; denys16_all.SE*2];
denys16.snorth = [denys16_pre09.SN; denys16_all.SN*2];
% denys16.nocc = [denys16_pre09.nocc; denys16_post09.nocc; denys16_all.nocc];
denys16.outlier = [denys16_pre09.outlier; denys16_all.outlier];

denys16 = struct2table(denys16);

page18.lons = [page18_pre09.Longitude; page18_post09.Longitude];
page18.lats = [page18_pre09.Latitude; page18_post09.Latitude];
page18.ID = [page18_pre09.site; page18_post09.site];
page18.veast_b16 = [page18_pre09.veast_b16; page18_post09.veast_b16];
page18.vnorth_b16 = [page18_pre09.vnorth_b16; page18_post09.vnorth_b16];
page18.seast = [page18_pre09.Se; page18_post09.Se];
page18.snorth = [page18_pre09.Sn; page18_post09.Sn];
page18.pagepre09 = [page18_pre09.Longitude*0 + 1; page18_post09.Longitude*0];
% page18.nocc = [page18_pre09.nocc; page18_post09.nocc];
page18.outlier = [page18_pre09.outlier; page18_post09.outlier];

% page18.lons = [page18_pre09.Longitude];
% page18.lats = [page18_pre09.Latitude];
% page18.ID = [page18_pre09.site];
% page18.veast_b16 = [page18_pre09.veast_b16];
% page18.vnorth_b16 = [page18_pre09.vnorth_b16];
% page18.seast = [page18_pre09.Se];
% page18.snorth = [page18_pre09.Sn];
% page18.pagepre09 = [page18_pre09.Longitude*0 + 1];
% % page18.nocc = [page18_pre09.nocc; page18_post09.nocc];
% page18.outlier = [page18_pre09.outlier];

page18 = struct2table(page18);
% denys16 = denys16_pre09;
% page18 = page18_pre09;

%%
for dummyindex = 1:length(beavan16.Longitude)
    beavan16.study{dummyindex} = 'beavan16';
end
for dummyindex = 1:length(denys16.lons)
    denys16.study{dummyindex} = 'denys16';
end
for dummyindex = 1:length(page18.lons)
    page18.study{dummyindex} = 'page18';
end
for dummyindex = 1:length(pickle19.Longitude)
    pickle19.study{dummyindex} = 'pickle19';
end
for dummyindex = 1:length(denys20.lon)
    denys20.study{dummyindex} = 'denys20';
end
for dummyindex = 1:length(hamling22.lon)
    hamling22.study{dummyindex} = 'hamling22';
end
for dummyindex = 1:length(midas.lons)
    midas.study{dummyindex} = 'midas';
end
for dummyindex = 1:length(measures.lons)
    measures.study{dummyindex} = 'measures';
end
for dummyindex = 1:length(beavan10.lon)
    beavan10.study{dummyindex} = 'beavan10';
end
for dummyindex = 1:length(beavan12.lon)
    beavan12.study{dummyindex} = 'beavan12';
end
for dummyindex = 1:length(houlie17.lon)
    houlie17.study{dummyindex} = 'houlie17';
    houlie17.study2{dummyindex} = 'houlie17_unadj';
end
% for dummyindex = 1:length(denys14.lons)
% denys14.study{dummyindex} = 'denys14';
% end
for dummyindex = 1:length(holden15.Lon)
    holden15.study{dummyindex} = 'holden15';
end
for dummyindex = 1:length(ulr5.Lon)
    ulr5.study{dummyindex} = 'ulr5';
end
for dummyindex = 1:length(UNR.lons)
    UNR.study{dummyindex} = 'UNR';
end
for dummyindex = 1:length(vardic21.lons)
    vardic21.study{dummyindex} = 'vardic21';
end
for dummyindex = 1:length(itrf14vlbi.lons)
    itrf14vlbi.study{dummyindex} = 'itrf14vlbi';
end
for dummyindex = 1:length(lamb13.lon)
    lamb13.study{dummyindex} = 'lamb13';
end


%%
gps.lons = [beavan16.Longitude; denys16.lons; page18.lons; pickle19.Longitude; denys20.lon; hamling22.lon; midas.lons; measures.lons; beavan10.lon; beavan12.lon; houlie17.lon; houlie17.lon; holden15.Lon; ulr5.Lon; UNR.lons; vardic21.lons; itrf14vlbi.lons; lamb13.lon];
gps.lats = [beavan16.Latitude; denys16.lats; page18.lats; pickle19.Latitude; denys20.lat; hamling22.lat; midas.lats; measures.lats; beavan10.lat; beavan12.lat; houlie17.lat; houlie17.lat; holden15.Lat; ulr5.Lat; UNR.lats; vardic21.lats; itrf14vlbi.lats; lamb13.lat];

gps.ve = [beavan16.VE_AUS; denys16.veast_b16; page18.veast_b16; pickle19.veast_b16; NaN*denys20.vu; NaN*hamling22.vu; midas.veast_b16; measures.veast_b16; NaN*beavan10.lon; NaN*beavan12.lon; NaN*houlie17.lon; NaN*houlie17.lon; holden15.VE; NaN*ulr5.Lon; UNR.veast_b16; vardic21.veast_b16; itrf14vlbi.veast_b16*NaN; NaN*lamb13.lon];
gps.vn = [beavan16.VN_AUS; denys16.vnorth_b16; page18.vnorth_b16; pickle19.vnorth_b16; NaN*denys20.vu; NaN*hamling22.vu; midas.vnorth_b16; measures.vnorth_b16; NaN*beavan10.lat; NaN*beavan12.lat; NaN*houlie17.lat; NaN*houlie17.lat; holden15.VN; NaN*ulr5.Lon; UNR.vnorth_b16; vardic21.vnorth_b16; itrf14vlbi.vnorth_b16*NaN; NaN*lamb13.lon];
gps.vu = [NaN*beavan16.VN_AUS; NaN*denys16.veast_b16; NaN*page18.vnorth_b16; pickle19.H; denys20.vu; hamling22.vu; midas.vup; measures.vup_itrf08; beavan10.vup_itrf08; beavan12.vu; houlie17.vu; houlie17.vu_unadj; holden15.VUP; ulr5.V_GPS; NaN*UNR.veast_b16; vardic21.vup_itrf08; itrf14vlbi.vup_itrf08; lamb13.vu_IT08];

gps.se = [beavan16.x_E; denys16.seast; page18.seast; pickle19.Eunc; NaN*denys20.vu; NaN*hamling22.vu; midas.seast; measures.seast; NaN*beavan10.lon; NaN*beavan12.lon; NaN*houlie17.vu; NaN*houlie17.vu; holden15.sVE; NaN*ulr5.Lon; UNR.seast; vardic21.seast; itrf14vlbi.seast*NaN; NaN*lamb13.se];
gps.sn = [beavan16.x_N; denys16.snorth; page18.snorth; pickle19.Nunc; NaN*denys20.vu; NaN*hamling22.vu; midas.snorth; measures.snorth; NaN*beavan10.lat; NaN*beavan12.lat; NaN*houlie17.vu; NaN*houlie17.vu; holden15.sVN; NaN*ulr5.Lon; UNR.snorth; vardic21.snorth; itrf14vlbi.snorth*NaN; NaN*lamb13.sn];
gps.su = [NaN*beavan16.VN_AUS; NaN*denys16.veast_b16; NaN*page18.vnorth_b16; pickle19.Hunc; denys20.su; hamling22.su; midas.su; measures.su; beavan10.su; beavan12.su; houlie17.su; houlie17.su; holden15.sVUP; ulr5.S_GPS; NaN*UNR.snorth; vardic21.sup; itrf14vlbi.sup; lamb13.su];

gps.outlier = [beavan16.outlier; denys16.outlier; page18.outlier; pickle19.outlier; denys20.outlier; hamling22.outlier; midas.outlier; measures.outlier; beavan10.outlier; beavan12.outlier; houlie17.outlier; houlie17.outlier; holden15.outlier; ulr5.outlier; UNR.outlier; vardic21.outlier; itrf14vlbi.outlier; lamb13.outlier];
gps.studies = [beavan16.study; denys16.study; page18.study; pickle19.study; denys20.study; hamling22.study; midas.study; measures.study; beavan10.study; beavan12.study; houlie17.study; houlie17.study2; holden15.study; ulr5.study; UNR.study; vardic21.study; itrf14vlbi.study; lamb13.study];
%return
%gps.nocc = [beavan16.NDAYS; ]
%hamling22.names = [];
gps.names = [beavan16.SITE; denys16.ID; page18.ID; pickle19.sitecode; denys20.name; hamling22.names; midas.names; measures.names; beavan10.name; beavan12.site; houlie17.name; houlie17.name; holden15.Site; ulr5.Site; UNR.names; vardic21.names; itrf14vlbi.names; lamb13.name];
gps.names = lower(gps.names);

gps.pagepre09 = [beavan16.outlier*0 + 1; denys16.outlier*0 + 1; page18.pagepre09; pickle19.outlier*0 + 1; denys20.outlier*0 + 1; hamling22.outlier*0 + 1; midas.outlier*0 + 1; measures.outlier*0 + 1; beavan10.outlier*0 + 1; beavan12.outlier*0 + 1; houlie17.outlier*0 + 1; houlie17.outlier*0 + 1; holden15.outlier*0 + 1; ulr5.outlier*0 + 1; UNR.outlier*0 + 1; vardic21.outlier*0 + 1; itrf14vlbi.outlier*0 + 1; lamb13.outlier*0 + 1];

gps.inhamling = strcmp(gps.studies,'hamling22');
gps.invlbi = strcmp(gps.studies,'itrf14vlbi');

%return
for gpsindex = 1:length(gps.lons)
    if and(gps.inhamling(gpsindex)==0,gps.invlbi(gpsindex)==0)
        gps.names{gpsindex} = gps.names{gpsindex}(1:4);
    end
end
gps = struct2table(gps);
%return
local_vels = gps.vu(and(strcmp(gps.names,'ousd'),and(strcmp(gps.studies,'beavan10')==0,and(~isnan(gps.vu),gps.outlier==0))));
local_unc = gps.su(and(strcmp(gps.names,'ousd'),and(strcmp(gps.studies,'beavan10')==0,and(~isnan(gps.vu),gps.outlier==0))));

local_weights = 1./local_unc.^2;

wmeanvel = sum(local_vels.*local_weights)/sum(local_weights);
velmat = local_vels + randn(length(local_weights),1e6).*local_unc;
weightmat = velmat*0 + local_weights;
velmat = velmat(:);
weightmat = weightmat(:);

wmeanste1 = sqrt(sum(local_weights.*(local_vels - wmeanvel).^2)/sum(local_weights))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
wmeanste = sqrt(sum(weightmat.*(velmat - wmeanvel).^2)/sum(weightmat))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
wmeanste2 = sqrt(1./sum(local_weights));
if or(wmeanste1>wmeanste,wmeanste2>wmeanste)
    pause
end

medste = std(local_vels)/(length(local_vels)^(1/4))*1.2533;

ousd_vu_avg = wmeanvel/2 + median(local_vels)/2;
ousd_su_avg = sqrt((wmeanste/2).^2 + (medste/2).^2);

gps.vu(strcmp(gps.studies,'beavan10')) = gps.vu(strcmp(gps.studies,'beavan10')) + ousd_vu_avg;
gps.su(strcmp(gps.studies,'beavan10')) = sqrt(gps.su(strcmp(gps.studies,'beavan10')).^2 + ousd_su_avg.^2);

gps_u.beavan16 = gps(strcmp(gps.studies,'beavan16'),:);
gps_u.pickle19 = gps(strcmp(gps.studies,'pickle19'),:);
gps_u.hamling22 = gps(strcmp(gps.studies,'hamling22'),:);
gps_u.vardic21 = gps(strcmp(gps.studies,'vardic21'),:);
gps_u.UNR = gps(strcmp(gps.studies,'UNR'),:);
gps_u.ulr5 = gps(strcmp(gps.studies,'ulr5'),:);
gps_u.houlie17 = gps(strcmp(gps.studies,'houlie17'),:);
gps_u.houlie17_unadj = gps(strcmp(gps.studies,'houlie17_unadj'),:);
gps_u.holden15 = gps(strcmp(gps.studies,'holden15'),:);
gps_u.beavan12 = gps(strcmp(gps.studies,'beavan12'),:);
gps_u.beavan10 = gps(strcmp(gps.studies,'beavan10'),:);
gps_u.measures = gps(strcmp(gps.studies,'measures'),:);
gps_u.midas = gps(strcmp(gps.studies,'midas'),:);
gps_u.denys20 = gps(strcmp(gps.studies,'denys20'),:);
gps_u.denys16 = gps(strcmp(gps.studies,'denys16'),:);
gps_u.page18 = gps(strcmp(gps.studies,'page18'),:);
gps_u.itrf14vlbi = gps(strcmp(gps.studies,'itrf14vlbi'),:);
gps_u.lamb13 = gps(strcmp(gps.studies,'lamb13'),:);
%return

%%
gps.touse = gps.invlbi==0;

gps.numcloseby = gps.lons*NaN;

for gpsindex = 1:length(gps.lons)
    if gps.touse(gpsindex)
        stadists = sqrt(((gps.lons - gps.lons(gpsindex)).*111.1.*cosd(gps.lats)).^2 + ((gps.lats - gps.lats(gpsindex))*111.1).^2);
        samestations = and(or(strcmp(gps.names{gpsindex},gps.names),gps.inhamling),and(stadists<=combdist,gps.touse));
        gps.numcloseby(gpsindex) = sum(samestations);
    end
end
%return
gps.numcloseby(gps.invlbi) = 1;
[~,sortorder] = sort(gps.numcloseby);
%return
gps = gps(flip(sortorder),:);

%%
gps.preciseloc_1e4 = or((round(gps.lons*1e3)/1e3)~=gps.lons,(round(gps.lats*1e3)/1e3)~=gps.lats);
gps.preciseloc_1e3 = or((round(gps.lons*1e2)/1e2)~=gps.lons,(round(gps.lats*1e2)/1e2)~=gps.lats);
gps.preciseloc_1e2 = or((round(gps.lons*1e1)/1e1)~=gps.lons,(round(gps.lats*1e1)/1e1)~=gps.lats);
for gpsindex = 1:length(gps.lons)
    if gps.touse(gpsindex)
        %gpsindex
        %pause
        stadists = sqrt(((gps.lons - gps.lons(gpsindex)).*111.1.*cosd(gps.lats)).^2 + ((gps.lats - gps.lats(gpsindex))*111.1).^2);
        samestations = and(or(strcmp(gps.names{gpsindex},gps.names),gps.inhamling),and(stadists<=combdist,gps.touse));
        if strcmp(gps.names{gpsindex},'noum')
            samestations = and(or(strcmp(gps.names{gpsindex},gps.names),gps.inhamling),and(stadists<=4,gps.touse));
        end
        samestations_1e4 = find(and(samestations,gps.preciseloc_1e4));
        samestations_1e3 = find(and(samestations,gps.preciseloc_1e3));
        samestations_1e2 = find(and(samestations,gps.preciseloc_1e2));
        samestations = find(samestations);

        if ~isempty(samestations_1e4)
            truelon = median(gps.lons(samestations_1e4));
            truelat = median(gps.lats(samestations_1e4));
        elseif ~isempty(samestations_1e3)
            truelon = median(gps.lons(samestations_1e3));
            truelat = median(gps.lats(samestations_1e3));
        elseif ~isempty(samestations_1e2)
            truelon = median(gps.lons(samestations_1e2));
            truelat = median(gps.lats(samestations_1e2));
        else
            truelon = median(gps.lons);
            truelat = median(gps.lats);
        end
    end
    gps.lons(samestations) = truelon;
    gps.lats(samestations) = truelat;
    gps.touse(samestations) = 0;
%    return
end
%return
%%
page18ind = find(and(and(strcmp(gps.studies,'page18'),gps.pagepre09==0),gps.outlier==0));
for fakeindex = 1:length(page18ind)
    gpsindex = page18ind(fakeindex);
    stadists = sqrt(((gps.lons - gps.lons(gpsindex)).*111.1.*cosd(gps.lats)).^2 + ((gps.lats - gps.lats(gpsindex))*111.1).^2);
    samestations = and(stadists<=combdist,gps.pagepre09==1);
    if sum(samestations)>0
        gps.outlier(gpsindex) = 1;
    end
end
sum(gps.outlier)

%return
%%
%outsidenz = (and(and(gps.lons>165,gps.lons<180),and(gps.lats>-48,gps.lats<-34))==0);
%outsidenz(strcmp(gps.names,'raul')) = 0;
%innorthland = and(gps.lats>-37.6,gps.lons<176);
%insouthland = or(and(gps.lats<-45.6,gps.lons>170.4),and(and(gps.lats<-44,gps.lons>170.4),strcmp(gps.names,'waim')));
%outsidenz_NS = or(outsidenz,or(innorthland,insouthland));
%outsidenz_NS = or(outsidenz,innorthland);

horizstudies = {'denys16';'page18';'pickle19';'holden15';'measures';'vardic21';'midas';'UNR'};
inbasis = and(strcmp(gps.studies,'beavan16'),gps.outlier==0);
%minsta = 12;

for studyindex = 1:length(horizstudies)

    horizstudies{studyindex}
    indicestocheckrange = and(strcmp(gps.studies,horizstudies{studyindex}),gps.outlier==0);
    localcheck = and(range(gps.lons(indicestocheckrange))<5,range(gps.lats(indicestocheckrange))<5);
    indicestorot = find(strcmp(gps.studies,horizstudies{studyindex}));
    localrot.lons = [];
    localrot.lats = [];
    localrot.dveast = [];
    localrot.dvnorth = [];
    localrot.sdveast = [];
    localrot.sdvnorth = [];
    localrot.gpsindex = [];
    localrot.remoteindex = [];

    for fakeindex = 1:length(indicestorot)
        gpsindex = indicestorot(fakeindex);
        if and(gps.outlier(gpsindex)==0,gps.pagepre09(gpsindex)==1)
        stadists = sqrt(((gps.lons - gps.lons(gpsindex)).*111.1.*cosd(gps.lats)).^2 + ((gps.lats - gps.lats(gpsindex))*111.1).^2);
        samestations = and(strcmp(gps.names{gpsindex},gps.names),stadists<=combdist);
        if strcmp(horizstudies{studyindex},'page18')
        sameinbasis = find(and(samestations,or(inbasis,and(strcmp(gps.studies,'denys16'),gps.outlier==0))));
        else
        sameinbasis = find(and(samestations,inbasis));
        end
        if ~isempty(sameinbasis)
            for dummyindex = 1:length(sameinbasis)
                localrot.lons = [localrot.lons; gps.lons(gpsindex)];
                localrot.lats = [localrot.lats; gps.lats(gpsindex)];
                localrot.dveast = [localrot.dveast; gps.ve(sameinbasis(dummyindex)) - gps.ve(gpsindex)];
                localrot.dvnorth = [localrot.dvnorth; gps.vn(sameinbasis(dummyindex)) - gps.vn(gpsindex)];
                localrot.sdveast = [localrot.sdveast; sqrt(gps.se(gpsindex).^2 + gps.se(sameinbasis(dummyindex)).^2)];
                localrot.sdvnorth = [localrot.sdvnorth; sqrt(gps.sn(gpsindex).^2 + gps.sn(sameinbasis(dummyindex)).^2)];
                localrot.gpsindex = [localrot.gpsindex; gpsindex];
                localrot.remoteindex = [localrot.remoteindex; sameinbasis(dummyindex)];
            end
        end
        end
    end
    localrot = struct2table(localrot);
    localrot_all = localrot;

    numcommonsta = length(localrot.dveast)

    d_alltouse = reshape([localrot.dveast'; localrot.dvnorth'],2*numcommonsta,1)*1e-3;
    sig_alltouse = reshape([localrot.sdveast'; localrot.sdvnorth'],2*numcommonsta,1)*1e-3;

    if localcheck==0
        G_alltouse = build_G_euler_m(localrot.lons,localrot.lats,numcommonsta);
        numtorot = length(indicestorot);
        G_alltorot = build_G_euler_m(gps.lons(indicestorot),gps.lats(indicestorot),numtorot);

        [x,y,z] = wgs2xyz(localrot.lons,localrot.lats,0*localrot.lons);
        Vxout = xyz2neu([x,y,z],[x*0 + 1 x*0 x*0],zeros(numcommonsta,3),zeros(numcommonsta,3));
        G_alltouse(:,4) = reshape([Vxout(:,2)'; Vxout(:,1)'],2*numcommonsta,1);
        Vyout = xyz2neu([x,y,z],[x*0 x*0 + 1 x*0],zeros(numcommonsta,3),zeros(numcommonsta,3));
        G_alltouse(:,5) = reshape([Vyout(:,2)'; Vyout(:,1)'],2*numcommonsta,1);
        Vzout = xyz2neu([x,y,z],[x*0 x*0 x*0 + 1],zeros(numcommonsta,3),zeros(numcommonsta,3));
        G_alltouse(:,6) = reshape([Vzout(:,2)'; Vzout(:,1)'],2*numcommonsta,1);

        [x,y,z] = wgs2xyz(gps.lons(indicestorot),gps.lats(indicestorot),0*indicestorot);
        Vxout = xyz2neu([x,y,z],[x*0 + 1 x*0 x*0],zeros(numtorot,3),zeros(numtorot,3));
        G_alltorot(:,4) = reshape([Vxout(:,2)'; Vxout(:,1)'],2*numtorot,1);
        Vyout = xyz2neu([x,y,z],[x*0 x*0 + 1 x*0],zeros(numtorot,3),zeros(numtorot,3));
        G_alltorot(:,5) = reshape([Vyout(:,2)'; Vyout(:,1)'],2*numtorot,1);
        Vzout = xyz2neu([x,y,z],[x*0 x*0 x*0 + 1],zeros(numtorot,3),zeros(numtorot,3));
        G_alltorot(:,6) = reshape([Vzout(:,2)'; Vzout(:,1)'],2*numtorot,1);
        cols1to3 = abs(G_alltouse(:,1:3));
        cols1to3 = mean(cols1to3(:));
        cols4to6 = abs(G_alltouse(:,4:6));
        cols4to6 = mean(cols4to6(:));
        G_alltouse(:,4:6) = G_alltouse(:,4:6)*cols1to3/cols4to6;
        G_alltorot(:,4:6) = G_alltorot(:,4:6)*cols1to3/cols4to6;

    else
        G_alltouse = [repmat([1;0],numcommonsta,1) repmat([0;1],numcommonsta,1)];
        numtorot = length(indicestorot);
        G_alltorot = [repmat([1;0],numtorot,1) repmat([0;1],numtorot,1)];
    end

    niter = 50000;

    indicesnottosave = [];
    for iteriter = 1:3

    dstar_local = NaN(numcommonsta*2,niter);
    for dummyindex = 1:niter
        lrnum = round(rand(1,1)*numcommonsta);
        if localcheck==0
        while lrnum<6
            lrnum = round(rand(1,1)*numcommonsta);
        end
        else
        while lrnum<2
            lrnum = round(rand(1,1)*numcommonsta);
        end
        end
        lrind = randperm(numcommonsta,lrnum)';
        lrind = sort([2*lrind-1; 2*lrind]);
        G_local = G_alltouse(lrind,:);
        d_local = d_alltouse(lrind);
        sig_local = sig_alltouse(lrind);

        mstar = (G_local'*diag(1./sig_local.^2)*G_local)\(G_local'*diag(1./sig_local.^2)*d_local);
        if rcond((G_local'*diag(1./sig_local.^2)*G_local))<eps
            disp('bad')
            dstar_local(:,dummyindex) = G_alltouse*mstar*NaN;
        else
            dstar_local(:,dummyindex) = G_alltouse*mstar;
        end
    end
    dstar_local = dstar_local(:,~isnan(dstar_local(1,:)));
    nrms = sqrt(sum(((dstar_local - d_alltouse)./sig_alltouse).^2)./(numcommonsta*2));
    %return
    dstar_local = dstar_local(:,nrms<=prctile(nrms,95));
    
    prc2resids_e = prctile(dstar_local(1:2:end,:) - d_alltouse(1:2:end),2.275,2)>(2*sig_alltouse(1:2:end));
    %return
    prc98resids_e = prctile(dstar_local(1:2:end,:) - d_alltouse(1:2:end),97.725,2)<(-2*sig_alltouse(1:2:end));
    prc2resids_n = prctile(dstar_local(2:2:end,:) - d_alltouse(2:2:end),2.275,2)>(2*sig_alltouse(2:2:end));    
    prc98resids_n = prctile(dstar_local(2:2:end,:) - d_alltouse(2:2:end),97.725,2)<(-2*sig_alltouse(2:2:end));
    
    localrot.tolose = or(or(prc2resids_e,prc98resids_e),or(prc2resids_n,prc98resids_n));

    indicesnottosave = [indicesnottosave; localrot.gpsindex(localrot.tolose)];
    indtolose = sort([2*find(localrot.tolose) - 1; 2*find(localrot.tolose)]);
    G_alltouse(indtolose,:) = [];
    d_alltouse(indtolose) = [];
    sig_alltouse(indtolose) = [];

    localrot(localrot.tolose,:) = [];

    numcommonsta = length(localrot.dveast)

    end

    dstar_local = NaN(numcommonsta*2,niter);
    dstar_rot = NaN(numtorot*2,niter);
    for dummyindex = 1:niter
        lrnum = round(rand(1,1)*numcommonsta);
        if localcheck==0
        while lrnum<6
            lrnum = round(rand(1,1)*numcommonsta);
        end
        else
        while lrnum<2
            lrnum = round(rand(1,1)*numcommonsta);
        end
        end
        lrind = randperm(numcommonsta,lrnum)';
        lrind = sort([2*lrind-1; 2*lrind]);
        G_local = G_alltouse(lrind,:);
        d_local = d_alltouse(lrind);
        sig_local = sig_alltouse(lrind);

        mstar = (G_local'*diag(1./sig_local.^2)*G_local)\(G_local'*diag(1./sig_local.^2)*d_local);
        if rcond((G_local'*diag(1./sig_local.^2)*G_local))<eps
            disp('bad')
            dstar_local(:,dummyindex) = G_alltouse*mstar*NaN;
            dstar_rot(:,dummyindex) = G_alltorot*mstar*NaN;
        else
            dstar_local(:,dummyindex) = G_alltouse*mstar;
            dstar_rot(:,dummyindex) = G_alltorot*mstar;
        end
    end

    dstar_local = dstar_local(:,~isnan(dstar_local(1,:)));
    dstar_rot = dstar_rot(:,~isnan(dstar_rot(1,:)));
    nrms = sqrt(sum(((dstar_local - d_alltouse)./sig_alltouse).^2)./(numcommonsta*2));
    clear dstar_local

    dstar_rot = dstar_rot(:,nrms<=prctile(nrms,95))*1e3;
    %dstar_local = dstar_local(:,nrms<=prctile(nrms,90))*1e3;

    dstar_rot_med = median(dstar_rot,2);
    dstar_rot_sig = mean(abs(dstar_rot - dstar_rot_med),2);
    % dstarprobs = exp(-sum(((dstar_local - (d_alltouse*1e3))./(sig_alltouse*1e3)).^2));
    % dstarprobs = dstarprobs/sum(dstarprobs);

%    dstar_local_med = sum(dstar_rot.*dstarprobs,2)/sum(dstarprobs);
%    dstar_local_sig = sum((abs(dstar_rot - dstar_local_med).*dstarprobs),2)/sum(dstarprobs);
    %pause
    % return
    gps.ve(indicestorot) = gps.ve(indicestorot) + dstar_rot_med(1:2:end);
    gps.vn(indicestorot) = gps.vn(indicestorot) + dstar_rot_med(2:2:end);
    %pause
    gps.se(indicestorot) = sqrt(gps.se(indicestorot).^2 + (dstar_rot_sig(1:2:end)).^2);
    gps.sn(indicestorot) = sqrt(gps.sn(indicestorot).^2 + (dstar_rot_sig(2:2:end)).^2);

    %sum(indicestosave)
    if strcmp(horizstudies{studyindex},'pickle19')
        inbasis(localrot.gpsindex) = 1;
    elseif strcmp(horizstudies{studyindex},'page18')
        gps.outlier(localrot.gpsindex) = 1;
%     elseif strcmp(horizstudies{studyindex},'measures')
% %        inbasis(and(strcmp(gps.studies,'measures'),gps.outlier==0)) = 1;
%         extratoadd = and(strcmp(gps.studies,horizstudies{studyindex}),gps.outlier==0);
%         extratoadd(indicesnottosave) = 0;
%         inbasis(extratoadd) = 1;
    end
    %return
    sum(inbasis)
    sc = 1;
    quiver(localrot_all.lons,localrot_all.lats,sc*localrot_all.dveast,sc*localrot_all.dvnorth,0,'linewidth',1.25)
hold on;
quiver(localrot.lons,localrot.lats,sc*localrot.dveast,sc*localrot.dvnorth,0,'linewidth',1.25)
quiver(gps.lons(indicestorot),gps.lats(indicestorot),sc*dstar_rot_med(1:2:end),sc*dstar_rot_med(2:2:end),0,'linewidth',1.25)
ellipse(sc*dstar_rot_sig(1:2:end),sc*dstar_rot_sig(2:2:end),0*gps.lons(indicestorot),gps.lons(indicestorot) + sc*dstar_rot_med(1:2:end),gps.lats(indicestorot) + sc*dstar_rot_med(2:2:end),[1 0.5 0],30);
    %[mean(dstar_rot_med(1:2:end)) mean(dstar_rot_med(2:2:end))]
    %pause


    close all
    clear localrot
end
%return
sum(gps.outlier)
% %%
% page18ind = find(strcmp(gps.studies,'page18'));
% for fakeindex = 1:length(page18ind)
%     gpsindex = page18ind(fakeindex);
%     stadists = sqrt(((gps.lons - gps.lons(gpsindex)).*111.1.*cosd(gps.lats)).^2 + ((gps.lats - gps.lats(gpsindex))*111.1).^2);
%     samestations = stadists<=combdist;
%     if sum(samestations)>0
%         gps.outlier(gpsindex) = 1;
%     end
% end
% sum(gps.outlier)

%%

%%
studies = {'beavan16','pickle19','hamling22','vardic21','UNR','ulr5','houlie17','houlie17_unadj','holden15','beavan12','beavan10','measures','midas','denys20','denys16','page18','itrf14vlbi','lamb13'};

clear gps_a
for studyindex = 1:length(studies)
     gps_a.(studies{studyindex}) = gps(strcmp(gps.studies,studies{studyindex}),:);
     export_indices = find(strcmp(gps.studies,studies{studyindex}));

     ffile = fopen(['NZ_GNSS_' studies{studyindex} '.dat'],'wt');
     fprintf(ffile,['#lon lat ve vn vu se sn su name\n']);
     for gpsindex = 1:length(export_indices)
          fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %s\n',...
               gps.lons(export_indices(gpsindex)),...
               gps.lats(export_indices(gpsindex)),...
               gps.ve(export_indices(gpsindex)),...
               gps.vn(export_indices(gpsindex)),...
               gps.vu(export_indices(gpsindex)),...
               gps.se(export_indices(gpsindex)),...
               gps.sn(export_indices(gpsindex)),...
               gps.su(export_indices(gpsindex)),...
               gps.names{export_indices(gpsindex)});
     end
     fclose(ffile);
    
end

%return
gps_a.(studies{studyindex})
%%
gps(gps.outlier==1,:) = [];

gps_c.beavan16 = gps(strcmp(gps.studies,'beavan16'),:);
gps_c.pickle19 = gps(strcmp(gps.studies,'pickle19'),:);
gps_c.hamling22 = gps(strcmp(gps.studies,'hamling22'),:);
gps_c.vardic21 = gps(strcmp(gps.studies,'vardic21'),:);
gps_c.UNR = gps(strcmp(gps.studies,'UNR'),:);
gps_c.ulr5 = gps(strcmp(gps.studies,'ulr5'),:);
gps_c.houlie17 = gps(strcmp(gps.studies,'houlie17'),:);
gps_c.houlie17_unadj = gps(strcmp(gps.studies,'houlie17_unadj'),:);
gps_c.holden15 = gps(strcmp(gps.studies,'holden15'),:);
gps_c.beavan12 = gps(strcmp(gps.studies,'beavan12'),:);
gps_c.beavan10 = gps(strcmp(gps.studies,'beavan10'),:);
gps_c.measures = gps(strcmp(gps.studies,'measures'),:);
gps_c.midas = gps(strcmp(gps.studies,'midas'),:);
gps_c.denys20 = gps(strcmp(gps.studies,'denys20'),:);
gps_c.denys16 = gps(strcmp(gps.studies,'denys16'),:);
gps_c.page18 = gps(strcmp(gps.studies,'page18'),:);
gps_c.itrf14vlbi = gps(strcmp(gps.studies,'itrf14vlbi'),:);
gps_c.lamb13 = gps(strcmp(gps.studies,'lamb13'),:);

%%
gps_m.names = [];
gps_m.lons = [];
gps_m.lats = [];
gps_m.ve = [];
gps_m.vn = [];
gps_m.vu = [];
gps_m.se = [];
gps_m.sn = [];
gps_m.su = [];

touse = ones(length(gps.lons),1);

for gpsindex = 1:length(gps.lons)
    if 1==touse(gpsindex)
        stadists = sqrt(((gps.lons - gps.lons(gpsindex)).*111.1.*cosd(gps.lats)).^2 + ((gps.lats - gps.lats(gpsindex))*111.1).^2);
        samestations = and(stadists<combdist,touse);
        if strcmp(gps.names(gpsindex),'h22_176')
            samestations(strcmp(gps.names,'a170')) = 1;
        end

        samestations = and(samestations,touse==1);
        samestations_withhoriz = and(and(samestations,and(and(~isnan(gps.ve),~isnan(gps.vn)),and(~isnan(gps.se),~isnan(gps.sn)))),touse==1);
        samestations_withvert = and(and(samestations,and(~isnan(gps.su),~isnan(gps.vu))),touse==1);

        samestations_1e4 = find(and(samestations,gps.preciseloc_1e4));
        samestations_1e3 = find(and(samestations,gps.preciseloc_1e3));
        samestations_1e2 = find(and(samestations,gps.preciseloc_1e2));
        samestations = find(samestations);

        if ~isempty(samestations_1e4)
            truelon = median(gps.lons(samestations_1e4));
            truelat = median(gps.lats(samestations_1e4));
        elseif ~isempty(samestations_1e3)
            truelon = median(gps.lons(samestations_1e3));
            truelat = median(gps.lats(samestations_1e3));
        elseif ~isempty(samestations_1e2)
            truelon = median(gps.lons(samestations_1e2));
            truelat = median(gps.lats(samestations_1e2));
        else
            truelon = median(gps.lons);
            truelat = median(gps.lats);
        end

        gps_m.lons(end+1) = truelon;
        gps_m.lats(end+1) = truelat;

        numcloseby_withhoriz = sum(samestations_withhoriz);
        numcloseby_withvert = sum(samestations_withvert);

        if numcloseby_withhoriz>1
            local_vels = gps.ve(samestations_withhoriz);
            local_unc = gps.se(samestations_withhoriz);
            local_weights = 1./local_unc.^2;

            wmeanvel = sum(local_vels.*local_weights)/sum(local_weights);
            velmat = local_vels + randn(length(local_weights),1e6).*local_unc;
            weightmat = velmat*0 + local_weights;
            velmat = velmat(:);
            weightmat = weightmat(:);

            wmeanste1 = sqrt(sum(local_weights.*(local_vels - wmeanvel).^2)/sum(local_weights))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
            wmeanste = sqrt(sum(weightmat.*(velmat - wmeanvel).^2)/sum(weightmat))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
            wmeanste2 = sqrt(1./sum(local_weights));
            if or(wmeanste1>wmeanste,wmeanste2>wmeanste)
                pause
            end

            medste = std(local_vels)/(numcloseby_withhoriz^(1/4))*1.2533;

            gps_m.ve(end+1) = wmeanvel/2 + median(local_vels)/2;
            gps_m.se(end+1) = sqrt((wmeanste/2).^2 + (medste/2).^2);

            local_vels = gps.vn(samestations_withhoriz);
            local_unc = gps.sn(samestations_withhoriz);
            local_weights = 1./local_unc.^2;

            wmeanvel = sum(local_vels.*local_weights)/sum(local_weights);
            velmat = local_vels + randn(length(local_weights),1e6).*local_unc;
            weightmat = velmat*0 + local_weights;
            velmat = velmat(:);
            weightmat = weightmat(:);

            wmeanste1 = sqrt(sum(local_weights.*(local_vels - wmeanvel).^2)/sum(local_weights))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
            wmeanste = sqrt(sum(weightmat.*(velmat - wmeanvel).^2)/sum(weightmat))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
            wmeanste2 = sqrt(1./sum(local_weights));
            if or(wmeanste1>wmeanste,wmeanste2>wmeanste)
                pause
            end

            medste = std(local_vels)/(numcloseby_withhoriz^(1/4))*1.2533;

            gps_m.vn(end+1) = wmeanvel/2 + median(local_vels)/2;
            gps_m.sn(end+1) = sqrt((wmeanste/2).^2 + (medste/2).^2);

        elseif numcloseby_withhoriz==1
            gps_m.ve(end+1) = gps.ve(samestations_withhoriz);
            gps_m.se(end+1) = gps.se(samestations_withhoriz);

            gps_m.vn(end+1) = gps.vn(samestations_withhoriz);
            gps_m.sn(end+1) = gps.sn(samestations_withhoriz);

        else
            gps_m.ve(end+1) = NaN;
            gps_m.se(end+1) = NaN;
            gps_m.vn(end+1) = NaN;
            gps_m.sn(end+1) = NaN;
        end

        if numcloseby_withvert>1
            local_vels = gps.vu(samestations_withvert);
            local_unc = gps.su(samestations_withvert);
            local_weights = 1./local_unc.^2;

            wmeanvel = sum(local_vels.*local_weights)/sum(local_weights);
            velmat = local_vels + randn(length(local_weights),1e6).*local_unc;
            weightmat = velmat*0 + local_weights;
            velmat = velmat(:);
            weightmat = weightmat(:);
            wmeanste1 = sqrt(sum(local_weights.*(local_vels - wmeanvel).^2)/sum(local_weights))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
            wmeanste = sqrt(sum(weightmat.*(velmat - wmeanvel).^2)/sum(weightmat))/((sum(local_weights)^2)/sum(local_weights.^2) - 1)^(1/4);
            wmeanste2 = sqrt(1./sum(local_weights));
            if or(wmeanste1>wmeanste,wmeanste2>wmeanste)
                pause
            end

            medste = std(local_vels)/(numcloseby_withvert^(1/4))*1.2533;

            gps_m.vu(end+1) = wmeanvel/2 + median(local_vels)/2;
            gps_m.su(end+1) = sqrt((wmeanste/2).^2 + (medste/2).^2);

        elseif numcloseby_withvert==1
            gps_m.vu(end+1) = gps.vu(samestations_withvert);
            gps_m.su(end+1) = gps.su(samestations_withvert);

 %           weights_u = 1/gps.su(samestations_withvert).^2;

        else
%            pause
             gps_m.vu(end+1) = NaN;
             gps_m.su(end+1) = NaN;
        end
        
        gps_m.names{end+1} = gps.names{gpsindex};
        touse(samestations) = 0;
    end
end
%return
%%
gps_m.names = gps_m.names';
gps_m.lons = gps_m.lons';
gps_m.lats = gps_m.lats';
gps_m.ve = gps_m.ve';
gps_m.vn = gps_m.vn';
gps_m.vu = gps_m.vu';
gps_m.se = gps_m.se';
gps_m.sn = gps_m.sn';
gps_m.su = gps_m.su';

gps_m = struct2table(gps_m);
%gps_temp = gps_m;
%return
%gps_m.innz = and(and(gps_m.lons>165,gps_m.lons<180),and(gps_m.lats>-48,gps_m.lats<-34));
%%
%gps_m = gps_m;
%return
%gps_m = gps_m_08;
gps_08 = gps_m;

O = [gps_08.lats,gps_08.lons,zeros(length(gps_08.lats),1)];
V = [gps_08.vn,gps_08.ve,0*gps_08.ve]*1e-3;
SV = [gps_08.sn,gps_08.se,0*gps_08.se]*1e-3;

[x,y,z] = wgs2xyz(gps_08.lons,gps_08.lats,zeros(length(gps_08.lats),1));
[vxyz,cxyz] = neu2xyz(O,V,SV,0*SV);
%[vxyz] = neu2xyz_old(O2,V);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2008',[],'ITRF2014',[],[]);
[itrf14out] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
gps_m.vn = itrf14out(:,1)*1e3;
gps_m.ve = itrf14out(:,2)*1e3;
gps_m.vu = itrf14out(:,3)*NaN;
%return
%%
is3D = ~isnan(gps_08.vu);
O = [gps_08.lats(is3D),gps_08.lons(is3D),zeros(length(gps_08.lats(is3D)),1)];
V = [gps_08.vn(is3D),gps_08.ve(is3D),gps_08.vu(is3D)]*1e-3;
SV = [gps_08.sn(is3D),gps_08.se(is3D),gps_08.su(is3D)]*1e-3;

[x,y,z] = wgs2xyz(gps_08.lons(is3D),gps_08.lats(is3D),zeros(length(gps_08.lats(is3D)),1));
[vxyz,cxyz] = neu2xyz(O,V,SV,0*SV);
%[vxyz] = neu2xyz_old(O2,V);
xyzout = itrstrafo([x,y,z,vxyz],'ITRF2008',[],'ITRF2014',[],[]);
[itrf14out] = xyz2neu([x,y,z],xyzout(:,4:6),SV,0*SV);
gps_m.vn(is3D) = itrf14out(:,1)*1e3;
gps_m.ve(is3D) = itrf14out(:,2)*1e3;
gps_m.vu(is3D) = itrf14out(:,3)*1e3;

%return
%%

%return
%%
ffile = fopen(['NZ_GNSS_comb_v10_sillleftin_' date '.dat'],'wt');
fprintf(ffile,['#lon lat ve vn vu se sn su name\n']);
for gpsindex = 1:length(gps_m.lons)
        fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %s\n',...
            gps_m.lons(gpsindex),...
            gps_m.lats(gpsindex),...
            gps_m.ve(gpsindex),...
            gps_m.vn(gpsindex),...
            gps_m.vu(gpsindex),...
            gps_m.se(gpsindex),...
            gps_m.sn(gpsindex),...
            gps_m.su(gpsindex),...
            gps_m.names{gpsindex});
end
fclose(ffile);

toexport_bfor = ~isnan(gps_m.ve);
ffile = fopen(['./NZ_GNSS_comb_v10_sillleftin_' date '_bfor.dat'],'wt');
for gpsindex = 1:length(gps_m.lons)
    if toexport_bfor(gpsindex)
        fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f\n',...
            gps_m.lons(gpsindex),...
            gps_m.lats(gpsindex),...
            gps_m.ve(gpsindex),...
            gps_m.vn(gpsindex),...
            gps_m.se(gpsindex),...
            gps_m.sn(gpsindex));
    end
end
fclose(ffile);

%%
gps_v = gps_m;

toexport_2D = isnan(gps_v.vu);
nohorizyesvert = and(isnan(gps_v.ve),~isnan(gps_v.vu));
%gps_v.vu(horizandnovert) = 0;
%gps_v.su(horizandnovert) = 1e6;
gps_v.ve(nohorizyesvert) = 0;
gps_v.se(nohorizyesvert) = 1e6;
gps_v.vn(nohorizyesvert) = 0;
gps_v.sn(nohorizyesvert) = 1e6;
toexport_3D = ~isnan(gps_v.vu);

ffile = fopen(['NZ_GNSS_comb_v10_sillleftin_' date '_velmap2D.dat'],'wt');
fprintf(ffile,['#lon lat ve vn se sn corr name\n']);
for gpsindex = 1:length(gps_v.lons)
    if toexport_2D(gpsindex)
        fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f %.0f %s\n',...
            gps_v.lons(gpsindex),...
            gps_v.lats(gpsindex),...
            gps_v.ve(gpsindex),...
            gps_v.vn(gpsindex),...
            gps_v.se(gpsindex),...
            gps_v.sn(gpsindex),...
            0,...
            gps_v.names{gpsindex});
    end
end
fclose(ffile);

ffile = fopen(['NZ_GNSS_comb_v10_sillleftin_' date '_velmap3D.dat'],'wt');
fprintf(ffile,['#lon lat ve vn vu se sn su corr corr corr name\n']);
for gpsindex = 1:length(gps_v.lons)
    if toexport_3D(gpsindex)
        fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.0f %.0f %.0f %s\n',...
            gps_v.lons(gpsindex),...
            gps_v.lats(gpsindex),...
            gps_v.ve(gpsindex),...
            gps_v.vn(gpsindex),...
            gps_v.vu(gpsindex),...
            gps_v.se(gpsindex),...
            gps_v.sn(gpsindex),...
            gps_v.su(gpsindex),...
            0,...
            0,...
            0,...
            gps_v.names{gpsindex});
    end
end
fclose(ffile);

%%
gps_s = gps_m;
load forward_model_sill.dat
%forward_model_sill(and(forward_model_sill(:,1)==170.22790,forward_model_sill(:,2)==-45.20290),:) = [];
disttosillsta = NaN(length(gps_s.ve),1);
sillcorr = NaN(length(gps_s.ve),1);
for staindex = 1:length(gps_s.ve)
    [disttosillsta(staindex),minindex] = min(sqrt(((gps_s.lons(staindex) - forward_model_sill(:,1)).*111.1.*cosd(forward_model_sill(:,2))).^2 + ((gps_s.lats(staindex) - forward_model_sill(:,2)).*111.1).^2));
    if disttosillsta(staindex)>1
    %     [gps_s.lons(staindex),gps_s.lats(staindex),disttosillsta(staindex)]
        %staindex
    end
    sillcorr(staindex) = forward_model_sill(minindex,5);
    gps_s.ve(staindex) = gps_s.ve(staindex) - forward_model_sill(minindex,3);
    gps_s.vn(staindex) = gps_s.vn(staindex) - forward_model_sill(minindex,4);
    gps_s.vu(staindex) = gps_s.vu(staindex) - forward_model_sill(minindex,5);
    %nearestsillsta = 
end

%%
ffile = fopen(['NZ_GNSS_comb_v10_sillremoved_' date '.dat'],'wt');
fprintf(ffile,['#lon lat ve vn vu se sn su name\n']);
for gpsindex = 1:length(gps_s.lons)
        fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %s\n',...
            gps_s.lons(gpsindex),...
            gps_s.lats(gpsindex),...
            gps_s.ve(gpsindex),...
            gps_s.vn(gpsindex),...
            gps_s.vu(gpsindex),...
            gps_s.se(gpsindex),...
            gps_s.sn(gpsindex),...
            gps_s.su(gpsindex),...
            gps_s.names{gpsindex});
end
fclose(ffile);

toexport_bfor = ~isnan(gps_s.ve);
ffile = fopen(['./NZ_GNSS_comb_v10_sillremoved_' date '_bfor.dat'],'wt');
for gpsindex = 1:length(gps_s.lons)
    if toexport_bfor(gpsindex)
        fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f\n',...
            gps_s.lons(gpsindex),...
            gps_s.lats(gpsindex),...
            gps_s.ve(gpsindex),...
            gps_s.vn(gpsindex),...
            gps_s.se(gpsindex),...
            gps_s.sn(gpsindex));
    end
end
fclose(ffile);

%save(['NZ_GNSS_comb_v10_sillremoved_' date '.mat'],'gps_s')
% else
% load(['NZ_GNSS_comb_v10_sillremoved_' date '.mat'])
% end
%return
%%
% sc = 1e-2;
% figure(52); clf; hold on; box on;
% set(gcf,'position',[0 720 1305 720])
% %set(gcf,'position',[0 1025 1450 1025])
% set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
%; plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
%; plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
% %plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')
% 
% % scatter(gps_a.hamling22.lons,gps_a.hamling22.lats,7^2,gps_a.hamling22.vu,'filled')
% % scatter(gps_a.pickle19.lons,gps_a.pickle19.lats,10.5^2,gps_a.pickle19.vu,'linewidth',1.25)
% % scatter(gps_a.midas.lons,gps_a.midas.lats,14^2,gps_a.midas.vu,'linewidth',1.25)
% % scatter(gps_a.measures.lons,gps_a.measures.lats,17.5^2,gps_a.measures.vu,'linewidth',1.25)
% % scatter(gps_a.holden15.lons,gps_a.holden15.lats,21^2,gps_a.holden15.vu,'linewidth',1.25)
% % scatter(gps_a.lamb13.lons,gps_a.lamb13.lats,24.5^2,gps_a.lamb13.vu,'linewidth',1.25)
% % scatter(gps_a.beavan12.lons,gps_a.beavan12.lats,28^2,gps_a.beavan12.vu,'linewidth',1.25)
% % scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,31.5^2,gps_a.houlie17.vu,'linewidth',1.25)
% % scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,35^2,gps_a.houlie17_unadj.vu,'linewidth',1.25)
% % scatter(gps_a.denys20.lons,gps_a.denys20.lats,38.5^2,gps_a.denys20.vu,'linewidth',1.25)
% % scatter(gps_a.ulr5.lons,gps_a.ulr5.lats,42^2,gps_a.ulr5.vu,'linewidth',1.25)
% 
% quiver(gps_a.midas.lons,gps_a.midas.lats,sc*lonsc*gps_a.midas.ve,sc*gps_a.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-2,'autoscale','off');
% thisellipse = ellipse(sc*lonsc*gps_a.midas.se,sc*gps_a.midas.sn,0*gps_a.midas.se,gps_a.midas.lons + sc*lonsc*gps_a.midas.ve,gps_a.midas.lats + sc*gps_a.midas.vn,color_midas,30);
% set(thisellipse,'linewidth',0.25)
% %text(gps_a.midas.lons(gps_a.midas.outlier),gps_a.midas.lats(gps_a.midas.outlier),upper(gps_a.midas.names(gps_a.midas.outlier)))
% % 
% % quiver(gps_a.UNR.lons,gps_a.UNR.lats,sc*lonsc*gps_a.UNR.ve,sc*gps_a.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-2,'autoscale','off','color',color_unr);
% % thisellipse = ellipse(sc*lonsc*gps_a.UNR.se,sc*gps_a.UNR.sn,0*gps_a.UNR.se,gps_a.UNR.lons + sc*lonsc*gps_a.UNR.ve,gps_a.UNR.lats + sc*gps_a.UNR.vn,color_unr,30);
% % set(thisellipse,'linewidth',0.25)
% % %text(gps_a.UNR.lons(gps_a.UNR.outlier),gps_a.UNR.lats(gps_a.UNR.outlier),upper(gps_a.UNR.names(gps_a.UNR.outlier)))
% 
% quiver(gps_a.measures.lons,gps_a.measures.lats,sc*lonsc*gps_a.measures.ve,sc*gps_a.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-2,'autoscale','off');
% thisellipse = ellipse(sc*lonsc*gps_a.measures.se,sc*gps_a.measures.sn,0*gps_a.measures.se,gps_a.measures.lons + sc*lonsc*gps_a.measures.ve,gps_a.measures.lats + sc*gps_a.measures.vn,color_measures,30);
% set(thisellipse,'linewidth',0.25)
% %text(gps_a.measures.lons(gps_a.measures.outlier),gps_a.measures.lats(gps_a.measures.outlier),upper(gps_a.measures.names(gps_a.measures.outlier)))
% 
% quiver(gps_a.pickle19.lons,gps_a.pickle19.lats,sc*lonsc*gps_a.pickle19.ve,sc*gps_a.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-2,'autoscale','off');
% thisellipse = ellipse(sc*lonsc*gps_a.pickle19.se,sc*gps_a.pickle19.sn,0*gps_a.pickle19.se,gps_a.pickle19.lons + sc*lonsc*gps_a.pickle19.ve,gps_a.pickle19.lats + sc*gps_a.pickle19.vn,color_pickle19,30);
% set(thisellipse,'linewidth',0.25)
% %text(gps_a.pickle19.lons(gps_a.pickle19.outlier),gps_a.pickle19.lats(gps_a.pickle19.outlier),upper(gps_a.pickle19.names(gps_a.pickle19.outlier)))
% 
% quiver(gps_a.beavan16.lons,gps_a.beavan16.lats,sc*lonsc*gps_a.beavan16.ve,sc*gps_a.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-2,'autoscale','off');
% thisellipse = ellipse(sc*lonsc*gps_a.beavan16.se,sc*gps_a.beavan16.sn,0*gps_a.beavan16.ve,gps_a.beavan16.lons + sc*lonsc*gps_a.beavan16.ve,gps_a.beavan16.lats + sc*gps_a.beavan16.vn,color_beavan16,30);
% set(thisellipse,'linewidth',0.25)
% %text(gps_a.beavan16.lons(gps_a.beavan16.outlier),gps_a.beavan16.lats(gps_a.beavan16.outlier),upper(gps_a.beavan16.names(gps_a.beavan16.outlier)))
% 
% colormap(jet6_1000);
% clim([-2 2])
% hcb = colorbar;
% hcb.Location = "southoutside";
% %hcb.Label.String = "Uplift (mm/yr)";
% hcb.TickLabels{5} = "Uplift (mm/yr)";
% 
% xlim([110 230])
% ylim([-55 -10])
% fontname(gcf,'Minion Pro')
% ax = gca;
% ax.FontSize = 12.5;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) - 0.02;
% ax_width = outerpos(3) - 0.02;
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% %title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
% title(['Compiled horizontal and vertical GNSS velocities, aligned'],'fontsize',25)
% %return
% print(52,['NZGNSS_northsouth_' date '_aligned.png'],'-dpng');
% return

% %%
% sc = 7.5e-3;
% figure(23); clf; hold on; box on;
% set(gcf,'position',[0 720 1082 720])
% %set(gcf,'position',[0 1025 1450 1025])
% set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
%; plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
%; plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
% %plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')
% 
% scatter(gps_c.hamling22.lons,gps_c.hamling22.lats,7^2,gps_c.hamling22.vu,'filled')
% scatter(gps_c.pickle19.lons,gps_c.pickle19.lats,10.5^2,gps_c.pickle19.vu,'linewidth',1.25)
% scatter(gps_c.midas.lons,gps_c.midas.lats,14^2,gps_c.midas.vu,'linewidth',1.25)
% scatter(gps_c.measures.lons,gps_c.measures.lats,17.5^2,gps_c.measures.vu,'linewidth',1.25)
% scatter(gps_c.holden15.lons,gps_c.holden15.lats,21^2,gps_c.holden15.vu,'linewidth',1.25)
% scatter(gps_c.lamb13.lons,gps_c.lamb13.lats,24.5^2,gps_c.lamb13.vu,'linewidth',1.25)
% scatter(gps_c.beavan12.lons,gps_c.beavan12.lats,28^2,gps_c.beavan12.vu,'linewidth',1.25)
% scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,31.5^2,gps_c.houlie17.vu,'linewidth',1.25)
% scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,35^2,gps_c.houlie17_unadj.vu,'linewidth',1.25)
% scatter(gps_c.denys20.lons,gps_c.denys20.lats,38.5^2,gps_c.denys20.vu,'linewidth',1.25)
% scatter(gps_c.ulr5.lons,gps_c.ulr5.lats,42^2,gps_c.ulr5.vu,'linewidth',1.25)
% 
% % quiver(gps_c.lamb13.lons,gps_c.lamb13.lats,sc*lonsc*gps_c.lamb13.ve,sc*gps_c.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% % thisellipse = ellipse(sc*lonsc*gps_c.lamb13.se,sc*gps_c.lamb13.sn,0*gps_c.lamb13.se,gps_c.lamb13.lons + sc*lonsc*gps_c.lamb13.ve,gps_c.lamb13.lats + sc*gps_c.lamb13.vn,color_lamb13,30);
% % set(thisellipse,'linewidth',0.25)
% 
% quiver(gps_c.midas.lons,gps_c.midas.lats,sc*lonsc*gps_c.midas.ve,sc*gps_c.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
% thisellipse = ellipse(sc*lonsc*gps_c.midas.se,sc*gps_c.midas.sn,0*gps_c.midas.se,gps_c.midas.lons + sc*lonsc*gps_c.midas.ve,gps_c.midas.lats + sc*gps_c.midas.vn,color_midas,30);
% set(thisellipse,'linewidth',0.25)
% 
% quiver(gps_c.UNR.lons,gps_c.UNR.lats,sc*lonsc*gps_c.UNR.ve,sc*gps_c.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
% thisellipse = ellipse(sc*lonsc*gps_c.UNR.se,sc*gps_c.UNR.sn,0*gps_c.UNR.se,gps_c.UNR.lons + sc*lonsc*gps_c.UNR.ve,gps_c.UNR.lats + sc*gps_c.UNR.vn,color_unr,30);
% set(thisellipse,'linewidth',0.25)
% 
% quiver(gps_c.measures.lons,gps_c.measures.lats,sc*lonsc*gps_c.measures.ve,sc*gps_c.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
% thisellipse = ellipse(sc*lonsc*gps_c.measures.se,sc*gps_c.measures.sn,0*gps_c.measures.se,gps_c.measures.lons + sc*lonsc*gps_c.measures.ve,gps_c.measures.lats + sc*gps_c.measures.vn,color_measures,30);
% set(thisellipse,'linewidth',0.25)
% 
% quiver(gps_c.pickle19.lons,gps_c.pickle19.lats,sc*lonsc*gps_c.pickle19.ve,sc*gps_c.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
% thisellipse = ellipse(sc*lonsc*gps_c.pickle19.se,sc*gps_c.pickle19.sn,0*gps_c.pickle19.se,gps_c.pickle19.lons + sc*lonsc*gps_c.pickle19.ve,gps_c.pickle19.lats + sc*gps_c.pickle19.vn,color_pickle19,30);
% set(thisellipse,'linewidth',0.25)
% %text(gps_u.measures.lons,gps_u.measures.lats,upper(gps_u.measures.names))
% 
% quiver(gps_c.beavan16.lons,gps_c.beavan16.lats,sc*lonsc*gps_c.beavan16.ve,sc*gps_c.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
% thisellipse = ellipse(sc*lonsc*gps_c.beavan16.se,sc*gps_c.beavan16.sn,0*gps_c.beavan16.ve,gps_c.beavan16.lons + sc*lonsc*gps_c.beavan16.ve,gps_c.beavan16.lats + sc*gps_c.beavan16.vn,color_beavan16,30);
% set(thisellipse,'linewidth',0.25)
% %text(gps_c.beavan16.lons(gps_c.beavan16.outlier),gps_c.beavan16.lats(gps_c.beavan16.outlier),upper(gps_c.beavan16.names(gps_c.beavan16.outlier)))
% 
% colormap(jet6_1000);
% clim([-2 2])
% hcb = colorbar;
% hcb.Location = "southoutside";
% %hcb.Label.String = "Uplift (mm/yr)";
% hcb.TickLabels{5} = "Uplift (mm/yr)";
% 
% xlim([171.25 177.25])
% ylim([-42 -39.25])
% fontname(gcf,'Minion Pro')
% ax = gca;
% ax.FontSize = 12.5;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) - 0.019;
% ax_width = outerpos(3) - 0.0325;
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% %title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
% title(['Compiled horizontal and vertical GNSS velocities, aligned, outliers removed'],'fontsize',25)
% %return
% print(23,['NZGNSS_northsouth_' date '_aligned_culled.png'],'-dpng');
% %return
% 
% %%
% sc = 7.5e-3;
% figure(24); clf; hold on; box on;
% set(gcf,'position',[0 720 1082 720])
% %set(gcf,'position',[0 1025 1450 1025])
% set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
%; plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
%; plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
% %plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')
% 
% scatter(gps_s.lons(~isnan(gps_s.vu)),gps_s.lats(~isnan(gps_s.vu)),9^2,gps_s.vu(~isnan(gps_s.vu)),'filled','markeredgecolor','none')
% scatter(gps_s.lons(~isnan(gps_s.vu)),gps_s.lats(~isnan(gps_s.vu)),14^2,gps_s.vu(~isnan(gps_s.vu)) - gps_s.su(~isnan(gps_s.vu)),'linewidth',1.25)
% scatter(gps_s.lons(~isnan(gps_s.vu)),gps_s.lats(~isnan(gps_s.vu)),19.5^2,gps_s.vu(~isnan(gps_s.vu)) + gps_s.su(~isnan(gps_s.vu)),'linewidth',1.25)
% gpsplot = quiver(gps_s.lons,gps_s.lats,sc/cosd(-41)*gps_s.ve,sc*gps_s.vn,0,'color','k','linewidth',1.25,'maxheadsize',1e-3);
% gpsellipse = ellipse(sc/cosd(-41)*gps_s.se,sc*gps_s.sn,0*gps_s.sn,gps_s.lons + sc/cosd(-41)*gps_s.ve,gps_s.lats + sc*gps_s.vn,'k',30);
% set(gpsellipse,'linewidth',0.25)
% 
% colormap(jet6_1000);
% clim([-2 2])
% hcb = colorbar;
% hcb.Location = "southoutside";
% %hcb.Label.String = "Uplift (mm/yr)";
% hcb.TickLabels{5} = "Uplift (mm/yr)";
% 
% xlim([171.25 177.25])
% ylim([-42 -39.25])
% fontname(gcf,'Minion Pro')
% ax = gca;
% ax.FontSize = 12.5;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) - 0.019;
% ax_width = outerpos(3) - 0.0325;
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% %title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
% title(['Compiled horizontal and vertical GNSS vels., aligned, outliers removed, averaged'],'fontsize',25)
% %return
% print(24,['NZGNSS_northsouth_' date '_sillremoved.png'],'-dpng');
% return

%%
%%
color_beavan16 = 'k';
color_pickle19 = [0.3125 0.6875 0.125];
color_measures = [0.5 0 0.75];
color_midas = [0.825 0.575 0.325];
%color_midas = [0.9375 0.75 0];
%color_unr = [0.825 0.5875 0.35];
%color_unr = [0.7875 0.54375 0.3];
color_unr = [0.575 0.575 0.575];
color_holden15 = [1 0.6875 0.375];
color_vardic21 = [1 0.5 0.5];
%color_vardic21 = [1 0.625 0.5];
color_denys16 = [0 0.25 1];
color_page18 = [1 0.5 1];
color_itrf14vlbi = [1 0.25 0.25];

lonsc = 1/cosd(-41);

%%
sc = 3.75e-2;
figure(1); clf; hold on; box on;
set(gcf,'position',[0 720 924 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_u.hamling22.lons,gps_u.hamling22.lats,7^2,gps_u.hamling22.vu,'filled')
scatter(gps_u.pickle19.lons,gps_u.pickle19.lats,10.5^2,gps_u.pickle19.vu,'linewidth',1.25)
scatter(gps_u.midas.lons,gps_u.midas.lats,14^2,gps_u.midas.vu,'linewidth',1.25)
scatter(gps_u.measures.lons,gps_u.measures.lats,17.5^2,gps_u.measures.vu,'linewidth',1.25)
scatter(gps_u.holden15.lons,gps_u.holden15.lats,21^2,gps_u.holden15.vu,'linewidth',1.25)
scatter(gps_u.lamb13.lons,gps_u.lamb13.lats,24.5^2,gps_u.lamb13.vu,'linewidth',1.25)
scatter(gps_u.beavan12.lons,gps_u.beavan12.lats,28^2,gps_u.beavan12.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,31.5^2,gps_u.houlie17.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,35^2,gps_u.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_u.denys20.lons,gps_u.denys20.lats,38.5^2,gps_u.denys20.vu,'linewidth',1.25)
scatter(gps_u.ulr5.lons,gps_u.ulr5.lats,42^2,gps_u.ulr5.vu,'linewidth',1.25)
scatter(gps_u.vardic21.lons,gps_u.vardic21.lats,45.5^2,gps_u.vardic21.vu,'linewidth',1.25)
scatter(gps_u.itrf14vlbi.lons,gps_u.itrf14vlbi.lats,49^2,gps_u.itrf14vlbi.vu,'linewidth',1.25)

% quiver(gps_u.lamb13.lons,gps_u.lamb13.lats,sc*lonsc*gps_u.lamb13.ve,sc*gps_u.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% thisellipse = ellipse(sc*lonsc*gps_u.lamb13.se,sc*gps_u.lamb13.sn,0*gps_u.lamb13.se,gps_u.lamb13.lons + sc*lonsc*gps_u.lamb13.ve,gps_u.lamb13.lats + sc*gps_u.lamb13.vn,color_lamb13,30);
% set(thisellipse,'linewidth',0.25)

quiver(gps_u.midas.lons,gps_u.midas.lats,sc*lonsc*gps_u.midas.ve,sc*gps_u.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.midas.se,sc*gps_u.midas.sn,0*gps_u.midas.se,gps_u.midas.lons + sc*lonsc*gps_u.midas.ve,gps_u.midas.lats + sc*gps_u.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.midas.lons,gps_u.midas.lats,upper(gps_u.midas.names))

quiver(gps_u.UNR.lons,gps_u.UNR.lats,sc*lonsc*gps_u.UNR.ve,sc*gps_u.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_u.UNR.se,sc*gps_u.UNR.sn,0*gps_u.UNR.se,gps_u.UNR.lons + sc*lonsc*gps_u.UNR.ve,gps_u.UNR.lats + sc*gps_u.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.holden15.lons,gps_u.holden15.lats,sc*lonsc*gps_u.holden15.ve,sc*gps_u.holden15.vn,0,'linewidth',1.25,'maxheadsize',2e-2,'autoscale','off','color',color_holden15);
thisellipse = ellipse(sc*lonsc*gps_u.holden15.se,sc*gps_u.holden15.sn,0*gps_u.holden15.se,gps_u.holden15.lons + sc*lonsc*gps_u.holden15.ve,gps_u.holden15.lats + sc*gps_u.holden15.vn,color_holden15,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.measures.lons,gps_u.measures.lats,sc*lonsc*gps_u.measures.ve,sc*gps_u.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.measures.se,sc*gps_u.measures.sn,0*gps_u.measures.se,gps_u.measures.lons + sc*lonsc*gps_u.measures.ve,gps_u.measures.lats + sc*gps_u.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.pickle19.lons,gps_u.pickle19.lats,sc*lonsc*gps_u.pickle19.ve,sc*gps_u.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.pickle19.se,sc*gps_u.pickle19.sn,0*gps_u.pickle19.se,gps_u.pickle19.lons + sc*lonsc*gps_u.pickle19.ve,gps_u.pickle19.lats + sc*gps_u.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.beavan16.lons,gps_u.beavan16.lats,sc*lonsc*gps_u.beavan16.ve,sc*gps_u.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.beavan16.se,sc*gps_u.beavan16.sn,0*gps_u.beavan16.ve,gps_u.beavan16.lons + sc*lonsc*gps_u.beavan16.ve,gps_u.beavan16.lats + sc*gps_u.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.beavan16.lons(gps_u.beavan16.outlier),gps_u.beavan16.lats(gps_u.beavan16.outlier),upper(gps_u.beavan16.names(gps_u.beavan16.outlier)))

quiver(gps_u.vardic21.lons,gps_u.vardic21.lats,sc*lonsc*gps_u.vardic21.ve,sc*gps_u.vardic21.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_vardic21);
thisellipse = ellipse(sc*lonsc*gps_u.vardic21.se,sc*gps_u.vardic21.sn,0*gps_u.vardic21.se,gps_u.vardic21.lons + sc*lonsc*gps_u.vardic21.ve,gps_u.vardic21.lats + sc*gps_u.vardic21.vn,color_vardic21,30);
set(thisellipse,'linewidth',0.25)

text(177.5,-34.5 - 0.006,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275,-34.68975 - 0.006,'8 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075,-34.625,sc*lonsc*8,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.625 - 0.006,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.075,-34.75,sc*lonsc*8,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.75 - 0.006,'Pickle [2019+]','fontsize',13.5)
quiver(177.075,-34.875,sc*lonsc*8,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.875 - 0.006,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.075,-35,sc*lonsc*8,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35 - 0.006,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.075,-35.125,sc*lonsc*8,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35.125 - 0.006,'UNR [Kreemer et al., 2014]','fontsize',13.5)
quiver(177.075,-35.25,sc*lonsc*8,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35.25 - 0.006,'Holden et al. [2015]','fontsize',13.5)
quiver(177.075,-35.375,sc*lonsc*8,sc*0,0,'color',color_vardic21,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35.375 - 0.006,'Vardic et al. [2021]','fontsize',13.5)

text(177.5,-34.75 - 0.006 - 0.8125,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(177.44,-34.875 - 0.8125,7^2,-1,'filled')
text(177.5,-34.875 - 0.006 - 0.8125,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(177.426,-35 - 0.8125,10.5^2,-1,'linewidth',1.25)
text(177.5,-35 - 0.006 - 0.8125,['Pickle [2019+]'],'fontsize',13.5)
scatter(177.412,-35.125 - 0.8125,14^2,-1,'linewidth',1.25)
text(177.5,-35.125 - 0.006 - 0.8125,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(177.398,-35.25 - 0.8125,17.5^2,-1,'linewidth',1.25)
text(177.5,-35.25 - 0.006 - 0.8125,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(177.384,-35.375 - 0.8125,21^2,-1,'linewidth',1.25)
text(177.5,-35.375 - 0.006 - 0.8125,'Holden et al. [2015]','fontsize',13.5)
scatter(177.37,-35.5 - 0.8125,24.5^2,-1,'linewidth',1.25)
text(177.5,-35.5 - 0.006 - 0.8125,'Lamb and Smith [2013]','fontsize',13.5)
scatter(177.356,-35.625 - 0.8125,28^2,-1,'linewidth',1.25)
text(177.5,-35.625 - 0.006 - 0.8125,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(177.342,-35.75 - 0.8125,31.5^2,-1,'linewidth',1.25)
text(177.5,-35.75 - 0.006 - 0.8125,'Houlie & Stern [2017]','fontsize',13.5)
scatter(177.328,-35.875 - 0.8125,35^2,-1,'linewidth',1.25)
text(177.5,-35.875 - 0.006 - 0.8125,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(177.314,-36 - 0.8125,38.5^2,-1,'linewidth',1.25)
text(177.5,-36 - 0.006 - 0.8125,['Denys et al. [2020]'],'fontsize',13.5)
scatter(177.3,-36.125 - 0.8125,42^2,-1,'linewidth',1.25)
text(177.5,-36.125 - 0.006 - 0.8125,'ULR5 [2012]','fontsize',13.5)
scatter(177.286,-36.25 - 0.8125,45.5^2,-1,'linewidth',1.25)
text(177.5,-36.25 - 0.006 - 0.8125,['Vardic et al. [2021]'],'fontsize',13.5)
scatter(177.272,-36.375 - 0.8125,49^2,-1,'linewidth',1.25)
text(177.5,-36.375 - 0.006 - 0.8125,['ITRF2014 VLBI'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([172.5 178.75])
ylim([-37.75 -34.375])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0175;
ax_width = outerpos(3) - 0.0375;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities'],'fontsize',25)
%return
%print(1,['NZGNSS_aucknorth_' date '.png'],'-dpng');
%return
% 
%%
sc = 3.75e-2;
figure(2); clf; hold on; box on;
set(gcf,'position',[0 720 924 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_a.hamling22.lons,gps_a.hamling22.lats,7^2,gps_a.hamling22.vu,'filled')
scatter(gps_a.pickle19.lons,gps_a.pickle19.lats,10.5^2,gps_a.pickle19.vu,'linewidth',1.25)
scatter(gps_a.midas.lons,gps_a.midas.lats,14^2,gps_a.midas.vu,'linewidth',1.25)
scatter(gps_a.measures.lons,gps_a.measures.lats,17.5^2,gps_a.measures.vu,'linewidth',1.25)
scatter(gps_a.holden15.lons,gps_a.holden15.lats,21^2,gps_a.holden15.vu,'linewidth',1.25)
scatter(gps_a.lamb13.lons,gps_a.lamb13.lats,24.5^2,gps_a.lamb13.vu,'linewidth',1.25)
scatter(gps_a.beavan12.lons,gps_a.beavan12.lats,28^2,gps_a.beavan12.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,31.5^2,gps_a.houlie17.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,35^2,gps_a.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_a.denys20.lons,gps_a.denys20.lats,38.5^2,gps_a.denys20.vu,'linewidth',1.25)
scatter(gps_a.ulr5.lons,gps_a.ulr5.lats,42^2,gps_a.ulr5.vu,'linewidth',1.25)
scatter(gps_a.vardic21.lons,gps_a.vardic21.lats,45.5^2,gps_a.vardic21.vu,'linewidth',1.25)
scatter(gps_a.itrf14vlbi.lons,gps_a.itrf14vlbi.lats,49^2,gps_a.itrf14vlbi.vu,'linewidth',1.25)

quiver(gps_a.midas.lons,gps_a.midas.lats,sc*lonsc*gps_a.midas.ve,sc*gps_a.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.midas.se,sc*gps_a.midas.sn,0*gps_a.midas.se,gps_a.midas.lons + sc*lonsc*gps_a.midas.ve,gps_a.midas.lats + sc*gps_a.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.midas.lons(gps_a.midas.outlier),gps_a.midas.lats(gps_a.midas.outlier),upper(gps_a.midas.names(gps_a.midas.outlier)))

quiver(gps_a.UNR.lons,gps_a.UNR.lats,sc*lonsc*gps_a.UNR.ve,sc*gps_a.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_a.UNR.se,sc*gps_a.UNR.sn,0*gps_a.UNR.se,gps_a.UNR.lons + sc*lonsc*gps_a.UNR.ve,gps_a.UNR.lats + sc*gps_a.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.UNR.lons(gps_a.UNR.outlier),gps_a.UNR.lats(gps_a.UNR.outlier),upper(gps_a.UNR.names(gps_a.UNR.outlier)))

quiver(gps_a.holden15.lons,gps_a.holden15.lats,sc*lonsc*gps_a.holden15.ve,sc*gps_a.holden15.vn,0,'linewidth',1.25,'maxheadsize',2e-2,'autoscale','off','color',color_holden15);
thisellipse = ellipse(sc*lonsc*gps_a.holden15.se,sc*gps_a.holden15.sn,0*gps_a.holden15.se,gps_a.holden15.lons + sc*lonsc*gps_a.holden15.ve,gps_a.holden15.lats + sc*gps_a.holden15.vn,color_holden15,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.holden15.lons(gps_a.holden15.outlier),gps_a.holden15.lats(gps_a.holden15.outlier),upper(gps_a.holden15.names(gps_a.holden15.outlier)))

quiver(gps_a.measures.lons,gps_a.measures.lats,sc*lonsc*gps_a.measures.ve,sc*gps_a.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.measures.se,sc*gps_a.measures.sn,0*gps_a.measures.se,gps_a.measures.lons + sc*lonsc*gps_a.measures.ve,gps_a.measures.lats + sc*gps_a.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.measures.lons(gps_a.measures.outlier),gps_a.measures.lats(gps_a.measures.outlier),upper(gps_a.measures.names(gps_a.measures.outlier)))

quiver(gps_a.pickle19.lons,gps_a.pickle19.lats,sc*lonsc*gps_a.pickle19.ve,sc*gps_a.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.pickle19.se,sc*gps_a.pickle19.sn,0*gps_a.pickle19.se,gps_a.pickle19.lons + sc*lonsc*gps_a.pickle19.ve,gps_a.pickle19.lats + sc*gps_a.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.pickle19.lons(gps_a.pickle19.outlier),gps_a.pickle19.lats(gps_a.pickle19.outlier),upper(gps_a.pickle19.names(gps_a.pickle19.outlier)))
%text(gps_a.pickle19.lons,gps_a.pickle19.lats,upper(gps_a.pickle19.names))

quiver(gps_a.beavan16.lons,gps_a.beavan16.lats,sc*lonsc*gps_a.beavan16.ve,sc*gps_a.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.beavan16.se,sc*gps_a.beavan16.sn,0*gps_a.beavan16.ve,gps_a.beavan16.lons + sc*lonsc*gps_a.beavan16.ve,gps_a.beavan16.lats + sc*gps_a.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.beavan16.lons(gps_a.beavan16.outlier),gps_a.beavan16.lats(gps_a.beavan16.outlier),upper(gps_a.beavan16.names(gps_a.beavan16.outlier)))

quiver(gps_a.vardic21.lons,gps_a.vardic21.lats,sc*lonsc*gps_a.vardic21.ve,sc*gps_a.vardic21.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_vardic21);
thisellipse = ellipse(sc*lonsc*gps_a.vardic21.se,sc*gps_a.vardic21.sn,0*gps_a.vardic21.se,gps_a.vardic21.lons + sc*lonsc*gps_a.vardic21.ve,gps_a.vardic21.lats + sc*gps_a.vardic21.vn,color_vardic21,30);
set(thisellipse,'linewidth',0.25)

text(177.5,-34.5 - 0.006,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275,-34.68975 - 0.006,'8 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075,-34.625,sc*lonsc*8,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.625 - 0.006,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.075,-34.75,sc*lonsc*8,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.75 - 0.006,'Pickle [2019+]','fontsize',13.5)
quiver(177.075,-34.875,sc*lonsc*8,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.875 - 0.006,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.075,-35,sc*lonsc*8,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35 - 0.006,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.075,-35.125,sc*lonsc*8,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35.125 - 0.006,'UNR [Kreemer et al., 2014]','fontsize',13.5)
quiver(177.075,-35.25,sc*lonsc*8,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35.25 - 0.006,'Holden et al. [2015]','fontsize',13.5)
quiver(177.075,-35.375,sc*lonsc*8,sc*0,0,'color',color_vardic21,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35.375 - 0.006,'Vardic et al. [2021]','fontsize',13.5)

text(177.5,-34.75 - 0.006 - 0.8125,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(177.44,-34.875 - 0.8125,7^2,-1,'filled')
text(177.5,-34.875 - 0.006 - 0.8125,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(177.426,-35 - 0.8125,10.5^2,-1,'linewidth',1.25)
text(177.5,-35 - 0.006 - 0.8125,['Pickle [2019+]'],'fontsize',13.5)
scatter(177.412,-35.125 - 0.8125,14^2,-1,'linewidth',1.25)
text(177.5,-35.125 - 0.006 - 0.8125,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(177.398,-35.25 - 0.8125,17.5^2,-1,'linewidth',1.25)
text(177.5,-35.25 - 0.006 - 0.8125,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(177.384,-35.375 - 0.8125,21^2,-1,'linewidth',1.25)
text(177.5,-35.375 - 0.006 - 0.8125,'Holden et al. [2015]','fontsize',13.5)
scatter(177.37,-35.5 - 0.8125,24.5^2,-1,'linewidth',1.25)
text(177.5,-35.5 - 0.006 - 0.8125,'Lamb and Smith [2013]','fontsize',13.5)
scatter(177.356,-35.625 - 0.8125,28^2,-1,'linewidth',1.25)
text(177.5,-35.625 - 0.006 - 0.8125,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(177.342,-35.75 - 0.8125,31.5^2,-1,'linewidth',1.25)
text(177.5,-35.75 - 0.006 - 0.8125,'Houlie & Stern [2017]','fontsize',13.5)
scatter(177.328,-35.875 - 0.8125,35^2,-1,'linewidth',1.25)
text(177.5,-35.875 - 0.006 - 0.8125,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(177.314,-36 - 0.8125,38.5^2,-1,'linewidth',1.25)
text(177.5,-36 - 0.006 - 0.8125,['Denys et al. [2020]'],'fontsize',13.5)
scatter(177.3,-36.125 - 0.8125,42^2,-1,'linewidth',1.25)
text(177.5,-36.125 - 0.006 - 0.8125,'ULR5 [2012]','fontsize',13.5)
scatter(177.286,-36.25 - 0.8125,45.5^2,-1,'linewidth',1.25)
text(177.5,-36.25 - 0.006 - 0.8125,['Vardic et al. [2021]'],'fontsize',13.5)
scatter(177.272,-36.375 - 0.8125,49^2,-1,'linewidth',1.25)
text(177.5,-36.375 - 0.006 - 0.8125,['ITRF2014 VLBI'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([172.5 178.75])
ylim([-37.75 -34.375])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0175;
ax_width = outerpos(3) - 0.0375;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned'],'fontsize',25)
%return
%print(2,['NZGNSS_aucknorth_' date '_aligned.png'],'-dpng');
%return
% 
%%
sc = 3.75e-2;
figure(3); clf; hold on; box on;
set(gcf,'position',[0 720 924 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_c.hamling22.lons,gps_c.hamling22.lats,7^2,gps_c.hamling22.vu,'filled')
scatter(gps_c.pickle19.lons,gps_c.pickle19.lats,10.5^2,gps_c.pickle19.vu,'linewidth',1.25)
scatter(gps_c.midas.lons,gps_c.midas.lats,14^2,gps_c.midas.vu,'linewidth',1.25)
scatter(gps_c.measures.lons,gps_c.measures.lats,17.5^2,gps_c.measures.vu,'linewidth',1.25)
scatter(gps_c.holden15.lons,gps_c.holden15.lats,21^2,gps_c.holden15.vu,'linewidth',1.25)
scatter(gps_c.lamb13.lons,gps_c.lamb13.lats,24.5^2,gps_c.lamb13.vu,'linewidth',1.25)
scatter(gps_c.beavan12.lons,gps_c.beavan12.lats,28^2,gps_c.beavan12.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,31.5^2,gps_c.houlie17.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,35^2,gps_c.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_c.denys20.lons,gps_c.denys20.lats,38.5^2,gps_c.denys20.vu,'linewidth',1.25)
scatter(gps_c.ulr5.lons,gps_c.ulr5.lats,42^2,gps_c.ulr5.vu,'linewidth',1.25)
scatter(gps_c.vardic21.lons,gps_c.vardic21.lats,45.5^2,gps_c.vardic21.vu,'linewidth',1.25)
scatter(gps_c.itrf14vlbi.lons,gps_c.itrf14vlbi.lats,49^2,gps_c.itrf14vlbi.vu,'linewidth',1.25)

quiver(gps_c.midas.lons,gps_c.midas.lats,sc*lonsc*gps_c.midas.ve,sc*gps_c.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.midas.se,sc*gps_c.midas.sn,0*gps_c.midas.se,gps_c.midas.lons + sc*lonsc*gps_c.midas.ve,gps_c.midas.lats + sc*gps_c.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.UNR.lons,gps_c.UNR.lats,sc*lonsc*gps_c.UNR.ve,sc*gps_c.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_c.UNR.se,sc*gps_c.UNR.sn,0*gps_c.UNR.se,gps_c.UNR.lons + sc*lonsc*gps_c.UNR.ve,gps_c.UNR.lats + sc*gps_c.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.holden15.lons,gps_c.holden15.lats,sc*lonsc*gps_c.holden15.ve,sc*gps_c.holden15.vn,0,'linewidth',1.25,'maxheadsize',2e-2,'autoscale','off','color',color_holden15);
thisellipse = ellipse(sc*lonsc*gps_c.holden15.se,sc*gps_c.holden15.sn,0*gps_c.holden15.se,gps_c.holden15.lons + sc*lonsc*gps_c.holden15.ve,gps_c.holden15.lats + sc*gps_c.holden15.vn,color_holden15,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.measures.lons,gps_c.measures.lats,sc*lonsc*gps_c.measures.ve,sc*gps_c.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.measures.se,sc*gps_c.measures.sn,0*gps_c.measures.se,gps_c.measures.lons + sc*lonsc*gps_c.measures.ve,gps_c.measures.lats + sc*gps_c.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.pickle19.lons,gps_c.pickle19.lats,sc*lonsc*gps_c.pickle19.ve,sc*gps_c.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.pickle19.se,sc*gps_c.pickle19.sn,0*gps_c.pickle19.se,gps_c.pickle19.lons + sc*lonsc*gps_c.pickle19.ve,gps_c.pickle19.lats + sc*gps_c.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.beavan16.lons,gps_c.beavan16.lats,sc*lonsc*gps_c.beavan16.ve,sc*gps_c.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.beavan16.se,sc*gps_c.beavan16.sn,0*gps_c.beavan16.ve,gps_c.beavan16.lons + sc*lonsc*gps_c.beavan16.ve,gps_c.beavan16.lats + sc*gps_c.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.vardic21.lons,gps_c.vardic21.lats,sc*lonsc*gps_c.vardic21.ve,sc*gps_c.vardic21.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_vardic21);
thisellipse = ellipse(sc*lonsc*gps_c.vardic21.se,sc*gps_c.vardic21.sn,0*gps_c.vardic21.se,gps_c.vardic21.lons + sc*lonsc*gps_c.vardic21.ve,gps_c.vardic21.lats + sc*gps_c.vardic21.vn,color_vardic21,30);
set(thisellipse,'linewidth',0.25)

text(177.5,-34.5 - 0.006,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275,-34.68975 - 0.006,'8 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075,-34.625,sc*lonsc*8,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.625 - 0.006,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.075,-34.75,sc*lonsc*8,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.75 - 0.006,'Pickle [2019+]','fontsize',13.5)
quiver(177.075,-34.875,sc*lonsc*8,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.875 - 0.006,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.075,-35,sc*lonsc*8,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35 - 0.006,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.075,-35.125,sc*lonsc*8,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35.125 - 0.006,'UNR [Kreemer et al., 2014]','fontsize',13.5)
quiver(177.075,-35.25,sc*lonsc*8,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35.25 - 0.006,'Holden et al. [2015]','fontsize',13.5)
quiver(177.075,-35.375,sc*lonsc*8,sc*0,0,'color',color_vardic21,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-35.375 - 0.006,'Vardic et al. [2021]','fontsize',13.5)

text(177.5,-34.75 - 0.006 - 0.8125,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(177.44,-34.875 - 0.8125,7^2,-1,'filled')
text(177.5,-34.875 - 0.006 - 0.8125,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(177.426,-35 - 0.8125,10.5^2,-1,'linewidth',1.25)
text(177.5,-35 - 0.006 - 0.8125,['Pickle [2019+]'],'fontsize',13.5)
scatter(177.412,-35.125 - 0.8125,14^2,-1,'linewidth',1.25)
text(177.5,-35.125 - 0.006 - 0.8125,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(177.398,-35.25 - 0.8125,17.5^2,-1,'linewidth',1.25)
text(177.5,-35.25 - 0.006 - 0.8125,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(177.384,-35.375 - 0.8125,21^2,-1,'linewidth',1.25)
text(177.5,-35.375 - 0.006 - 0.8125,'Holden et al. [2015]','fontsize',13.5)
scatter(177.37,-35.5 - 0.8125,24.5^2,-1,'linewidth',1.25)
text(177.5,-35.5 - 0.006 - 0.8125,'Lamb and Smith [2013]','fontsize',13.5)
scatter(177.356,-35.625 - 0.8125,28^2,-1,'linewidth',1.25)
text(177.5,-35.625 - 0.006 - 0.8125,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(177.342,-35.75 - 0.8125,31.5^2,-1,'linewidth',1.25)
text(177.5,-35.75 - 0.006 - 0.8125,'Houlie & Stern [2017]','fontsize',13.5)
scatter(177.328,-35.875 - 0.8125,35^2,-1,'linewidth',1.25)
text(177.5,-35.875 - 0.006 - 0.8125,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(177.314,-36 - 0.8125,38.5^2,-1,'linewidth',1.25)
text(177.5,-36 - 0.006 - 0.8125,['Denys et al. [2020]'],'fontsize',13.5)
scatter(177.3,-36.125 - 0.8125,42^2,-1,'linewidth',1.25)
text(177.5,-36.125 - 0.006 - 0.8125,'ULR5 [2012]','fontsize',13.5)
scatter(177.286,-36.25 - 0.8125,45.5^2,-1,'linewidth',1.25)
text(177.5,-36.25 - 0.006 - 0.8125,['Vardic et al. [2021]'],'fontsize',13.5)
scatter(177.272,-36.375 - 0.8125,49^2,-1,'linewidth',1.25)
text(177.5,-36.375 - 0.006 - 0.8125,['ITRF2014 VLBI'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([172.5 178.75])
ylim([-37.75 -34.375])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0175;
ax_width = outerpos(3) - 0.0375;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, outliers removed'],'fontsize',25)
%return
%print(3,['NZGNSS_aucknorth_' date '_aligned_culled.png'],'-dpng');
% %return
% 
%%
sc = 3.75e-2;
figure(4); clf; hold on; box on;
set(gcf,'position',[0 720 924 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),9^2,gps_m.vu(~isnan(gps_m.vu)),'filled','markeredgecolor','none')
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),14^2,gps_m.vu(~isnan(gps_m.vu)) - gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),19.5^2,gps_m.vu(~isnan(gps_m.vu)) + gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
quiver(gps_m.lons,gps_m.lats,sc/cosd(-41)*gps_m.ve,sc*gps_m.vn,0,'color','k','linewidth',1.25,'maxheadsize',1e-3);
gpsellipse = ellipse(sc/cosd(-41)*gps_m.se,sc*gps_m.sn,0*gps_m.sn,gps_m.lons + sc/cosd(-41)*gps_m.ve,gps_m.lats + sc*gps_m.vn,'k',30);
set(gpsellipse,'linewidth',0.25)

text(177.5,-34.5 - 0.006,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075,-34.625,sc*lonsc*8,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5,-34.625 - 0.006,'8 mm/yr','fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([172.5 178.75])
ylim([-37.75 -34.375])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0175;
ax_width = outerpos(3) - 0.0375;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, culled, averaged'],'fontsize',25)
%return
%print(4,['NZGNSS_aucknorth_' date '_averaged.png'],'-dpng');
%return
%%
sc = 7.5e-3;
figure(11); clf; hold on; box on;
set(gcf,'position',[0 720 1234 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_u.hamling22.lons,gps_u.hamling22.lats,7^2,gps_u.hamling22.vu,'filled')
scatter(gps_u.pickle19.lons,gps_u.pickle19.lats,10.5^2,gps_u.pickle19.vu,'linewidth',1.25)
scatter(gps_u.midas.lons,gps_u.midas.lats,14^2,gps_u.midas.vu,'linewidth',1.25)
scatter(gps_u.measures.lons,gps_u.measures.lats,17.5^2,gps_u.measures.vu,'linewidth',1.25)
scatter(gps_u.holden15.lons,gps_u.holden15.lats,21^2,gps_u.holden15.vu,'linewidth',1.25)
scatter(gps_u.lamb13.lons,gps_u.lamb13.lats,24.5^2,gps_u.lamb13.vu,'linewidth',1.25)
scatter(gps_u.beavan12.lons,gps_u.beavan12.lats,28^2,gps_u.beavan12.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,31.5^2,gps_u.houlie17.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,35^2,gps_u.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_u.denys20.lons,gps_u.denys20.lats,38.5^2,gps_u.denys20.vu,'linewidth',1.25)

% quiver(gps_u.lamb13.lons,gps_u.lamb13.lats,sc*lonsc*gps_u.lamb13.ve,sc*gps_u.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% thisellipse = ellipse(sc*lonsc*gps_u.lamb13.se,sc*gps_u.lamb13.sn,0*gps_u.lamb13.se,gps_u.lamb13.lons + sc*lonsc*gps_u.lamb13.ve,gps_u.lamb13.lats + sc*gps_u.lamb13.vn,color_lamb13,30);
% set(thisellipse,'linewidth',0.25)

quiver(gps_u.midas.lons,gps_u.midas.lats,sc*lonsc*gps_u.midas.ve,sc*gps_u.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.midas.se,sc*gps_u.midas.sn,0*gps_u.midas.se,gps_u.midas.lons + sc*lonsc*gps_u.midas.ve,gps_u.midas.lats + sc*gps_u.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.UNR.lons,gps_u.UNR.lats,sc*lonsc*gps_u.UNR.ve,sc*gps_u.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_u.UNR.se,sc*gps_u.UNR.sn,0*gps_u.UNR.se,gps_u.UNR.lons + sc*lonsc*gps_u.UNR.ve,gps_u.UNR.lats + sc*gps_u.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.holden15.lons,gps_u.holden15.lats,sc*lonsc*gps_u.holden15.ve,sc*gps_u.holden15.vn,0,'linewidth',1.25,'maxheadsize',2e-2,'autoscale','off','color',color_holden15);
thisellipse = ellipse(sc*lonsc*gps_u.holden15.se,sc*gps_u.holden15.sn,0*gps_u.holden15.se,gps_u.holden15.lons + sc*lonsc*gps_u.holden15.ve,gps_u.holden15.lats + sc*gps_u.holden15.vn,color_holden15,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.measures.lons,gps_u.measures.lats,sc*lonsc*gps_u.measures.ve,sc*gps_u.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.measures.se,sc*gps_u.measures.sn,0*gps_u.measures.se,gps_u.measures.lons + sc*lonsc*gps_u.measures.ve,gps_u.measures.lats + sc*gps_u.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.pickle19.lons,gps_u.pickle19.lats,sc*lonsc*gps_u.pickle19.ve,sc*gps_u.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.pickle19.se,sc*gps_u.pickle19.sn,0*gps_u.pickle19.se,gps_u.pickle19.lons + sc*lonsc*gps_u.pickle19.ve,gps_u.pickle19.lats + sc*gps_u.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.pickle19.lons,gps_u.pickle19.lats,upper(gps_u.pickle19.names))

quiver(gps_u.beavan16.lons,gps_u.beavan16.lats,sc*lonsc*gps_u.beavan16.ve,sc*gps_u.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.beavan16.se,sc*gps_u.beavan16.sn,0*gps_u.beavan16.ve,gps_u.beavan16.lons + sc*lonsc*gps_u.beavan16.ve,gps_u.beavan16.lats + sc*gps_u.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.beavan16.lons(gps_u.beavan16.outlier),gps_u.beavan16.lats(gps_u.beavan16.outlier),upper(gps_u.beavan16.names(gps_u.beavan16.outlier)))

text(177.5 - 3.25,-34.5 - 0.005 - 3.0875,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275 - 3.25,-34.68975 - 0.005 - 3.0875+0.075,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.085 - 3.25,-34.625 - 3.0875+0.05,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.625 - 0.005 - 3.0875+0.05,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.085 - 3.25,-34.75 - 3.0875+0.1,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.75 - 0.005 - 3.0875+0.1,'Pickle [2019+]','fontsize',13.5)
quiver(177.085 - 3.25,-34.875 - 3.0875+0.15,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.875 - 0.005 - 3.0875+0.15,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.085 - 3.25,-35 - 3.0875+0.2,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-35 - 0.005 - 3.0875+0.2,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.085 - 3.25,-35.125 - 3.0875+0.25,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-35.125 - 0.005 - 3.0875+0.25,'UNR [Kreemer et al., 2014]','fontsize',13.5)
quiver(177.085 - 3.25,-35.25 - 3.0875+0.3,sc*lonsc*40,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-35.25 - 0.005 - 3.0875+0.3,'Holden et al. [2015]','fontsize',13.5)

text(177.5 - 3.465,-38.15,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(174.215-0.215,-38.225,7^2,-1,'filled')
text(177.5 - 3.465,-38.225 - 0.005,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(174.207-0.215,-38.3,10.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.3 - 0.005,['Pickle [2019+]'],'fontsize',13.5)
scatter(174.199-0.215,-38.375,14^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.375 - 0.005,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(174.191-0.215,-38.45,17.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.45 - 0.005,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(174.183-0.215,-38.525,21^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.525 - 0.005,'Holden et al. [2015]','fontsize',13.5)
scatter(174.175-0.215,-38.6,24.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.6 - 0.005,'Lamb and Smith [2013]','fontsize',13.5)
scatter(174.167-0.215,-38.675,28^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.675 - 0.005,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(174.159-0.215,-38.75,31.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.75 - 0.005,'Houlie & Stern [2017]','fontsize',13.5)
scatter(174.151-0.215,-38.825,35^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.825 - 0.005,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(174.143-0.215,-38.9,38.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.9 - 0.005,['Denys et al. [2020]'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([173.75 178.75])
ylim([-39.5 -37.5])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0195;
ax_width = outerpos(3) - 0.0275;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities'],'fontsize',25)
%return
%print(11,['NZGNSS_north_' date '.png'],'-dpng');
%return
% 
%%
sc = 7.5e-3;
figure(12); clf; hold on; box on;
set(gcf,'position',[0 720 1234 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_a.hamling22.lons,gps_a.hamling22.lats,7^2,gps_a.hamling22.vu,'filled')
scatter(gps_a.pickle19.lons,gps_a.pickle19.lats,10.5^2,gps_a.pickle19.vu,'linewidth',1.25)
scatter(gps_a.midas.lons,gps_a.midas.lats,14^2,gps_a.midas.vu,'linewidth',1.25)
scatter(gps_a.measures.lons,gps_a.measures.lats,17.5^2,gps_a.measures.vu,'linewidth',1.25)
scatter(gps_a.holden15.lons,gps_a.holden15.lats,21^2,gps_a.holden15.vu,'linewidth',1.25)
scatter(gps_a.lamb13.lons,gps_a.lamb13.lats,24.5^2,gps_a.lamb13.vu,'linewidth',1.25)
scatter(gps_a.beavan12.lons,gps_a.beavan12.lats,28^2,gps_a.beavan12.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,31.5^2,gps_a.houlie17.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,35^2,gps_a.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_a.denys20.lons,gps_a.denys20.lats,38.5^2,gps_a.denys20.vu,'linewidth',1.25)
scatter(gps_a.ulr5.lons,gps_a.ulr5.lats,42^2,gps_a.ulr5.vu,'linewidth',1.25)

quiver(gps_a.midas.lons,gps_a.midas.lats,sc*lonsc*gps_a.midas.ve,sc*gps_a.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.midas.se,sc*gps_a.midas.sn,0*gps_a.midas.se,gps_a.midas.lons + sc*lonsc*gps_a.midas.ve,gps_a.midas.lats + sc*gps_a.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.midas.lons(gps_a.midas.outlier),gps_a.midas.lats(gps_a.midas.outlier),upper(gps_a.midas.names(gps_a.midas.outlier)))
%text(gps_a.midas.lons,gps_a.midas.lats,upper(gps_a.midas.names))

quiver(gps_a.UNR.lons,gps_a.UNR.lats,sc*lonsc*gps_a.UNR.ve,sc*gps_a.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_a.UNR.se,sc*gps_a.UNR.sn,0*gps_a.UNR.se,gps_a.UNR.lons + sc*lonsc*gps_a.UNR.ve,gps_a.UNR.lats + sc*gps_a.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.UNR.lons(gps_a.UNR.outlier),gps_a.UNR.lats(gps_a.UNR.outlier),upper(gps_a.UNR.names(gps_a.UNR.outlier)))

quiver(gps_a.holden15.lons,gps_a.holden15.lats,sc*lonsc*gps_a.holden15.ve,sc*gps_a.holden15.vn,0,'linewidth',1.25,'maxheadsize',2e-2,'autoscale','off','color',color_holden15);
thisellipse = ellipse(sc*lonsc*gps_a.holden15.se,sc*gps_a.holden15.sn,0*gps_a.holden15.se,gps_a.holden15.lons + sc*lonsc*gps_a.holden15.ve,gps_a.holden15.lats + sc*gps_a.holden15.vn,color_holden15,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.holden15.lons(gps_a.holden15.outlier),gps_a.holden15.lats(gps_a.holden15.outlier),upper(gps_a.holden15.names(gps_a.holden15.outlier)))

quiver(gps_a.measures.lons,gps_a.measures.lats,sc*lonsc*gps_a.measures.ve,sc*gps_a.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.measures.se,sc*gps_a.measures.sn,0*gps_a.measures.se,gps_a.measures.lons + sc*lonsc*gps_a.measures.ve,gps_a.measures.lats + sc*gps_a.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.measures.lons(gps_a.measures.outlier),gps_a.measures.lats(gps_a.measures.outlier),upper(gps_a.measures.names(gps_a.measures.outlier)))

quiver(gps_a.pickle19.lons,gps_a.pickle19.lats,sc*lonsc*gps_a.pickle19.ve,sc*gps_a.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.pickle19.se,sc*gps_a.pickle19.sn,0*gps_a.pickle19.se,gps_a.pickle19.lons + sc*lonsc*gps_a.pickle19.ve,gps_a.pickle19.lats + sc*gps_a.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.pickle19.lons(gps_a.pickle19.outlier),gps_a.pickle19.lats(gps_a.pickle19.outlier),upper(gps_a.pickle19.names(gps_a.pickle19.outlier)))

quiver(gps_a.beavan16.lons,gps_a.beavan16.lats,sc*lonsc*gps_a.beavan16.ve,sc*gps_a.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.beavan16.se,sc*gps_a.beavan16.sn,0*gps_a.beavan16.ve,gps_a.beavan16.lons + sc*lonsc*gps_a.beavan16.ve,gps_a.beavan16.lats + sc*gps_a.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.beavan16.lons(gps_a.beavan16.outlier),gps_a.beavan16.lats(gps_a.beavan16.outlier),upper(gps_a.beavan16.names(gps_a.beavan16.outlier)))
%text(gps_a.beavan16.lons,gps_a.beavan16.lats,upper(gps_a.beavan16.names))

text(177.5 - 3.25,-34.5 - 0.005 - 3.0875,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275 - 3.25,-34.68975 - 0.005 - 3.0875+0.075,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.085 - 3.25,-34.625 - 3.0875+0.05,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.625 - 0.005 - 3.0875+0.05,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.085 - 3.25,-34.75 - 3.0875+0.1,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.75 - 0.005 - 3.0875+0.1,'Pickle [2019+]','fontsize',13.5)
quiver(177.085 - 3.25,-34.875 - 3.0875+0.15,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.875 - 0.005 - 3.0875+0.15,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.085 - 3.25,-35 - 3.0875+0.2,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-35 - 0.005 - 3.0875+0.2,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.085 - 3.25,-35.125 - 3.0875+0.25,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-35.125 - 0.005 - 3.0875+0.25,'UNR [Kreemer et al., 2014]','fontsize',13.5)
quiver(177.085 - 3.25,-35.25 - 3.0875+0.3,sc*lonsc*40,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-35.25 - 0.005 - 3.0875+0.3,'Holden et al. [2015]','fontsize',13.5)

text(177.5 - 3.465,-38.15,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(174.215-0.215,-38.225,7^2,-1,'filled')
text(177.5 - 3.465,-38.225 - 0.005,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(174.207-0.215,-38.3,10.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.3 - 0.005,['Pickle [2019+]'],'fontsize',13.5)
scatter(174.199-0.215,-38.375,14^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.375 - 0.005,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(174.191-0.215,-38.45,17.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.45 - 0.005,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(174.183-0.215,-38.525,21^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.525 - 0.005,'Holden et al. [2015]','fontsize',13.5)
scatter(174.175-0.215,-38.6,24.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.6 - 0.005,'Lamb and Smith [2013]','fontsize',13.5)
scatter(174.167-0.215,-38.675,28^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.675 - 0.005,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(174.159-0.215,-38.75,31.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.75 - 0.005,'Houlie & Stern [2017]','fontsize',13.5)
scatter(174.151-0.215,-38.825,35^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.825 - 0.005,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(174.143-0.215,-38.9,38.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.9 - 0.005,['Denys et al. [2020]'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([173.75 178.75])
ylim([-39.5 -37.5])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0195;
ax_width = outerpos(3) - 0.0275;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned'],'fontsize',25)
%return
%print(12,['NZGNSS_north_' date '_aligned.png'],'-dpng');
%return
% 
%%
sc = 7.5e-3;
figure(13); clf; hold on; box on;
set(gcf,'position',[0 720 1234 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_c.hamling22.lons,gps_c.hamling22.lats,7^2,gps_c.hamling22.vu,'filled')
scatter(gps_c.pickle19.lons,gps_c.pickle19.lats,10.5^2,gps_c.pickle19.vu,'linewidth',1.25)
scatter(gps_c.midas.lons,gps_c.midas.lats,14^2,gps_c.midas.vu,'linewidth',1.25)
scatter(gps_c.measures.lons,gps_c.measures.lats,17.5^2,gps_c.measures.vu,'linewidth',1.25)
scatter(gps_c.holden15.lons,gps_c.holden15.lats,21^2,gps_c.holden15.vu,'linewidth',1.25)
scatter(gps_c.lamb13.lons,gps_c.lamb13.lats,24.5^2,gps_c.lamb13.vu,'linewidth',1.25)
scatter(gps_c.beavan12.lons,gps_c.beavan12.lats,28^2,gps_c.beavan12.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,31.5^2,gps_c.houlie17.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,35^2,gps_c.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_c.denys20.lons,gps_c.denys20.lats,38.5^2,gps_c.denys20.vu,'linewidth',1.25)
scatter(gps_c.ulr5.lons,gps_c.ulr5.lats,42^2,gps_c.ulr5.vu,'linewidth',1.25)

% quiver(gps_c.lamb13.lons,gps_c.lamb13.lats,sc*lonsc*gps_c.lamb13.ve,sc*gps_c.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% thisellipse = ellipse(sc*lonsc*gps_c.lamb13.se,sc*gps_c.lamb13.sn,0*gps_c.lamb13.se,gps_c.lamb13.lons + sc*lonsc*gps_c.lamb13.ve,gps_c.lamb13.lats + sc*gps_c.lamb13.vn,color_lamb13,30);
% set(thisellipse,'linewidth',0.25)

quiver(gps_c.midas.lons,gps_c.midas.lats,sc*lonsc*gps_c.midas.ve,sc*gps_c.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.midas.se,sc*gps_c.midas.sn,0*gps_c.midas.se,gps_c.midas.lons + sc*lonsc*gps_c.midas.ve,gps_c.midas.lats + sc*gps_c.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.UNR.lons,gps_c.UNR.lats,sc*lonsc*gps_c.UNR.ve,sc*gps_c.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_c.UNR.se,sc*gps_c.UNR.sn,0*gps_c.UNR.se,gps_c.UNR.lons + sc*lonsc*gps_c.UNR.ve,gps_c.UNR.lats + sc*gps_c.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.holden15.lons,gps_c.holden15.lats,sc*lonsc*gps_c.holden15.ve,sc*gps_c.holden15.vn,0,'linewidth',1.25,'maxheadsize',2e-2,'autoscale','off','color',color_holden15);
thisellipse = ellipse(sc*lonsc*gps_c.holden15.se,sc*gps_c.holden15.sn,0*gps_c.holden15.se,gps_c.holden15.lons + sc*lonsc*gps_c.holden15.ve,gps_c.holden15.lats + sc*gps_c.holden15.vn,color_holden15,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.measures.lons,gps_c.measures.lats,sc*lonsc*gps_c.measures.ve,sc*gps_c.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.measures.se,sc*gps_c.measures.sn,0*gps_c.measures.se,gps_c.measures.lons + sc*lonsc*gps_c.measures.ve,gps_c.measures.lats + sc*gps_c.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.pickle19.lons,gps_c.pickle19.lats,sc*lonsc*gps_c.pickle19.ve,sc*gps_c.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.pickle19.se,sc*gps_c.pickle19.sn,0*gps_c.pickle19.se,gps_c.pickle19.lons + sc*lonsc*gps_c.pickle19.ve,gps_c.pickle19.lats + sc*gps_c.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.beavan16.lons,gps_c.beavan16.lats,sc*lonsc*gps_c.beavan16.ve,sc*gps_c.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.beavan16.se,sc*gps_c.beavan16.sn,0*gps_c.beavan16.ve,gps_c.beavan16.lons + sc*lonsc*gps_c.beavan16.ve,gps_c.beavan16.lats + sc*gps_c.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)

text(177.5 - 3.25,-34.5 - 0.005 - 3.0875,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275 - 3.25,-34.68975 - 0.005 - 3.0875+0.075,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.085 - 3.25,-34.625 - 3.0875+0.05,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.625 - 0.005 - 3.0875+0.05,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.085 - 3.25,-34.75 - 3.0875+0.1,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.75 - 0.005 - 3.0875+0.1,'Pickle [2019+]','fontsize',13.5)
quiver(177.085 - 3.25,-34.875 - 3.0875+0.15,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.875 - 0.005 - 3.0875+0.15,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.085 - 3.25,-35 - 3.0875+0.2,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-35 - 0.005 - 3.0875+0.2,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.085 - 3.25,-35.125 - 3.0875+0.25,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-35.125 - 0.005 - 3.0875+0.25,'UNR [Kreemer et al., 2014]','fontsize',13.5)
quiver(177.085 - 3.25,-35.25 - 3.0875+0.3,sc*lonsc*40,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-35.25 - 0.005 - 3.0875+0.3,'Holden et al. [2015]','fontsize',13.5)

text(177.5 - 3.465,-38.15,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(174.215-0.215,-38.225,7^2,-1,'filled')
text(177.5 - 3.465,-38.225 - 0.005,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(174.207-0.215,-38.3,10.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.3 - 0.005,['Pickle [2019+]'],'fontsize',13.5)
scatter(174.199-0.215,-38.375,14^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.375 - 0.005,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(174.191-0.215,-38.45,17.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.45 - 0.005,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(174.183-0.215,-38.525,21^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.525 - 0.005,'Holden et al. [2015]','fontsize',13.5)
scatter(174.175-0.215,-38.6,24.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.6 - 0.005,'Lamb and Smith [2013]','fontsize',13.5)
scatter(174.167-0.215,-38.675,28^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.675 - 0.005,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(174.159-0.215,-38.75,31.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.75 - 0.005,'Houlie & Stern [2017]','fontsize',13.5)
scatter(174.151-0.215,-38.825,35^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.825 - 0.005,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(174.143-0.215,-38.9,38.5^2,-1,'linewidth',1.25)
text(177.5 - 3.465,-38.9 - 0.005,['Denys et al. [2020]'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([173.75 178.75])
ylim([-39.5 -37.5])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0195;
ax_width = outerpos(3) - 0.0275;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, outliers removed'],'fontsize',25)
%return
%print(13,['NZGNSS_north_' date '_aligned_culled.png'],'-dpng');
%return

%%
sc = 7.5e-3;
figure(14); clf; hold on; box on;
set(gcf,'position',[0 720 1234 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),9^2,gps_m.vu(~isnan(gps_m.vu)),'filled','markeredgecolor','none')
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),14^2,gps_m.vu(~isnan(gps_m.vu)) - gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),19.5^2,gps_m.vu(~isnan(gps_m.vu)) + gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
quiver(gps_m.lons,gps_m.lats,sc/cosd(-41)*gps_m.ve,sc*gps_m.vn,0,'color','k','linewidth',1.25,'maxheadsize',1e-3);
gpsellipse = ellipse(sc/cosd(-41)*gps_m.se,sc*gps_m.sn,0*gps_m.sn,gps_m.lons + sc/cosd(-41)*gps_m.ve,gps_m.lats + sc*gps_m.vn,'k',30);
set(gpsellipse,'linewidth',0.25)

text(177.5 - 3.25,-34.5 - 0.005 - 3.0875,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
%text(177.275 - 3.25,-34.68975 - 0.005 - 3.0875+0.075,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.085 - 3.25,-34.625 - 3.0875+0.05,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25,-34.625 - 0.005 - 3.0875+0.05,'40 mm/yr','fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([173.75 178.75])
ylim([-39.5 -37.5])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0195;
ax_width = outerpos(3) - 0.0275;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, culled, averaged'],'fontsize',25)
%return
%print(14,['NZGNSS_north_' date '_averaged.png'],'-dpng');


%%
sc = 7.5e-3;
figure(21); clf; hold on; box on;
set(gcf,'position',[0 720 1082 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_u.hamling22.lons,gps_u.hamling22.lats,7^2,gps_u.hamling22.vu,'filled')
scatter(gps_u.pickle19.lons,gps_u.pickle19.lats,10.5^2,gps_u.pickle19.vu,'linewidth',1.25)
scatter(gps_u.midas.lons,gps_u.midas.lats,14^2,gps_u.midas.vu,'linewidth',1.25)
scatter(gps_u.measures.lons,gps_u.measures.lats,17.5^2,gps_u.measures.vu,'linewidth',1.25)
%scatter(gps_u.holden15.lons,gps_u.holden15.lats,21^2,gps_u.holden15.vu,'linewidth',1.25)
scatter(gps_u.lamb13.lons,gps_u.lamb13.lats,24.5^2,gps_u.lamb13.vu,'linewidth',1.25)
scatter(gps_u.beavan12.lons,gps_u.beavan12.lats,28^2,gps_u.beavan12.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,31.5^2,gps_u.houlie17.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,35^2,gps_u.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_u.denys20.lons,gps_u.denys20.lats,38.5^2,gps_u.denys20.vu,'linewidth',1.25)
scatter(gps_u.ulr5.lons,gps_u.ulr5.lats,42^2,gps_u.ulr5.vu,'linewidth',1.25)

% quiver(gps_u.lamb13.lons,gps_u.lamb13.lats,sc*lonsc*gps_u.lamb13.ve,sc*gps_u.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% thisellipse = ellipse(sc*lonsc*gps_u.lamb13.se,sc*gps_u.lamb13.sn,0*gps_u.lamb13.se,gps_u.lamb13.lons + sc*lonsc*gps_u.lamb13.ve,gps_u.lamb13.lats + sc*gps_u.lamb13.vn,color_lamb13,30);
% set(thisellipse,'linewidth',0.25)

quiver(gps_u.midas.lons,gps_u.midas.lats,sc*lonsc*gps_u.midas.ve,sc*gps_u.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.midas.se,sc*gps_u.midas.sn,0*gps_u.midas.se,gps_u.midas.lons + sc*lonsc*gps_u.midas.ve,gps_u.midas.lats + sc*gps_u.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.midas.lons,gps_u.midas.lats,upper(gps_u.midas.names))

quiver(gps_u.UNR.lons,gps_u.UNR.lats,sc*lonsc*gps_u.UNR.ve,sc*gps_u.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_u.UNR.se,sc*gps_u.UNR.sn,0*gps_u.UNR.se,gps_u.UNR.lons + sc*lonsc*gps_u.UNR.ve,gps_u.UNR.lats + sc*gps_u.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.measures.lons,gps_u.measures.lats,sc*lonsc*gps_u.measures.ve,sc*gps_u.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.measures.se,sc*gps_u.measures.sn,0*gps_u.measures.se,gps_u.measures.lons + sc*lonsc*gps_u.measures.ve,gps_u.measures.lats + sc*gps_u.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.pickle19.lons,gps_u.pickle19.lats,sc*lonsc*gps_u.pickle19.ve,sc*gps_u.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.pickle19.se,sc*gps_u.pickle19.sn,0*gps_u.pickle19.se,gps_u.pickle19.lons + sc*lonsc*gps_u.pickle19.ve,gps_u.pickle19.lats + sc*gps_u.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.pickle19.lons,gps_u.pickle19.lats,upper(gps_u.pickle19.names))

quiver(gps_u.beavan16.lons,gps_u.beavan16.lats,sc*lonsc*gps_u.beavan16.ve,sc*gps_u.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.beavan16.se,sc*gps_u.beavan16.sn,0*gps_u.beavan16.ve,gps_u.beavan16.lons + sc*lonsc*gps_u.beavan16.ve,gps_u.beavan16.lats + sc*gps_u.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.beavan16.lons(gps_u.beavan16.outlier),gps_u.beavan16.lats(gps_u.beavan16.outlier),upper(gps_u.beavan16.names(gps_u.beavan16.outlier)))

text(177.5 - 5.75,-34.5 - 0.005 - 4.85,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275 - 5.75,-34.68975 - 0.005 - 4.85+0.0375,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075 - 5.75,-34.625 - 4.85+0.025,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.625 - 0.005 - 4.85+0.025,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.075 - 5.75,-34.75 - 4.85+0.05,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.75 - 0.005 - 4.85+0.05,'Pickle [2019+]','fontsize',13.5)
quiver(177.075 - 5.75,-34.875 - 4.85+0.075,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.875 - 0.005 - 4.85+0.075,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.075 - 5.75,-35 - 4.85+0.1,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-35 - 0.005 - 4.85+0.1,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.075 - 5.75,-35.125 - 4.85+0.125,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-35.125 - 0.005 - 4.85+0.125,'UNR [Kreemer et al., 2014]','fontsize',13.5)
%quiver(177.075 - 5.75,-35.25 - 4.85+0.15,sc*lonsc*40,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
%text(177.5 - 5.75,-35.25 - 0.005 - 4.85+0.15,'Holden et al. [2015]','fontsize',13.5)

text(173,-34.5 - 0.005 - 4.85,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(173 - 0.05,-34.5 - 4.95,7^2,-1,'filled')
text(173,-34.5 - 0.005 - 4.95,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(173 - 0.061,-34.5 - 5.05,10.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.05,['Pickle [2019+]'],'fontsize',13.5)
scatter(173 - 0.072,-34.5 - 5.15,14^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.15,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(173 - 0.083,-34.5 - 5.25,17.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.25,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(173 - 0.094,-34.5 - 5.35,21^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.35,'Holden et al. [2015]','fontsize',13.5)
scatter(173 - 0.105,-34.5 - 5.45,24.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.45,'Lamb and Smith [2013]','fontsize',13.5)
scatter(173 - 0.116,-34.5 - 5.55,28^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.55,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(173 - 0.127,-34.5 - 5.65,31.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.65,'Houlie & Stern [2017]','fontsize',13.5)
scatter(173 - 0.138,-34.5 - 5.75,35^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.75,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(173 - 0.149,-34.5 - 5.85,38.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.85,['Denys et al. [2020]'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([171.25 177.25])
ylim([-42 -39.25])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.019;
ax_width = outerpos(3) - 0.0325;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
title(['Compiled horizontal and vertical GNSS velocities'],'fontsize',25)
%return
%print(21,['NZGNSS_northsouth_' date '.png'],'-dpng');
%return

%%
sc = 7.5e-3;
figure(22); clf; hold on; box on;
set(gcf,'position',[0 720 1082 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_a.hamling22.lons,gps_a.hamling22.lats,7^2,gps_a.hamling22.vu,'filled')
scatter(gps_a.pickle19.lons,gps_a.pickle19.lats,10.5^2,gps_a.pickle19.vu,'linewidth',1.25)
scatter(gps_a.midas.lons,gps_a.midas.lats,14^2,gps_a.midas.vu,'linewidth',1.25)
scatter(gps_a.measures.lons,gps_a.measures.lats,17.5^2,gps_a.measures.vu,'linewidth',1.25)
scatter(gps_a.holden15.lons,gps_a.holden15.lats,21^2,gps_a.holden15.vu,'linewidth',1.25)
scatter(gps_a.lamb13.lons,gps_a.lamb13.lats,24.5^2,gps_a.lamb13.vu,'linewidth',1.25)
scatter(gps_a.beavan12.lons,gps_a.beavan12.lats,28^2,gps_a.beavan12.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,31.5^2,gps_a.houlie17.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,35^2,gps_a.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_a.denys20.lons,gps_a.denys20.lats,38.5^2,gps_a.denys20.vu,'linewidth',1.25)
scatter(gps_a.ulr5.lons,gps_a.ulr5.lats,42^2,gps_a.ulr5.vu,'linewidth',1.25)

quiver(gps_a.midas.lons,gps_a.midas.lats,sc*lonsc*gps_a.midas.ve,sc*gps_a.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.midas.se,sc*gps_a.midas.sn,0*gps_a.midas.se,gps_a.midas.lons + sc*lonsc*gps_a.midas.ve,gps_a.midas.lats + sc*gps_a.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_a.UNR.lons,gps_a.UNR.lats,sc*lonsc*gps_a.UNR.ve,sc*gps_a.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_a.UNR.se,sc*gps_a.UNR.sn,0*gps_a.UNR.se,gps_a.UNR.lons + sc*lonsc*gps_a.UNR.ve,gps_a.UNR.lats + sc*gps_a.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.UNR.lons(gps_a.UNR.outlier),gps_a.UNR.lats(gps_a.UNR.outlier),upper(gps_a.UNR.names(gps_a.UNR.outlier)))

quiver(gps_a.measures.lons,gps_a.measures.lats,sc*lonsc*gps_a.measures.ve,sc*gps_a.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.measures.se,sc*gps_a.measures.sn,0*gps_a.measures.se,gps_a.measures.lons + sc*lonsc*gps_a.measures.ve,gps_a.measures.lats + sc*gps_a.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.measures.lons(gps_a.measures.outlier),gps_a.measures.lats(gps_a.measures.outlier),upper(gps_a.measures.names(gps_a.measures.outlier)))

quiver(gps_a.pickle19.lons,gps_a.pickle19.lats,sc*lonsc*gps_a.pickle19.ve,sc*gps_a.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.pickle19.se,sc*gps_a.pickle19.sn,0*gps_a.pickle19.se,gps_a.pickle19.lons + sc*lonsc*gps_a.pickle19.ve,gps_a.pickle19.lats + sc*gps_a.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.pickle19.lons(gps_a.pickle19.outlier),gps_a.pickle19.lats(gps_a.pickle19.outlier),upper(gps_a.pickle19.names(gps_a.pickle19.outlier)))

quiver(gps_a.beavan16.lons,gps_a.beavan16.lats,sc*lonsc*gps_a.beavan16.ve,sc*gps_a.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.beavan16.se,sc*gps_a.beavan16.sn,0*gps_a.beavan16.ve,gps_a.beavan16.lons + sc*lonsc*gps_a.beavan16.ve,gps_a.beavan16.lats + sc*gps_a.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.beavan16.lons(gps_a.beavan16.outlier),gps_a.beavan16.lats(gps_a.beavan16.outlier),upper(gps_a.beavan16.names(gps_a.beavan16.outlier)))

text(177.5 - 5.75,-34.5 - 0.005 - 4.85,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275 - 5.75,-34.68975 - 0.005 - 4.85+0.0375,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075 - 5.75,-34.625 - 4.85+0.025,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.625 - 0.005 - 4.85+0.025,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.075 - 5.75,-34.75 - 4.85+0.05,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.75 - 0.005 - 4.85+0.05,'Pickle [2019+]','fontsize',13.5)
quiver(177.075 - 5.75,-34.875 - 4.85+0.075,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.875 - 0.005 - 4.85+0.075,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.075 - 5.75,-35 - 4.85+0.1,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-35 - 0.005 - 4.85+0.1,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.075 - 5.75,-35.125 - 4.85+0.125,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-35.125 - 0.005 - 4.85+0.125,'UNR [Kreemer et al., 2014]','fontsize',13.5)
%quiver(177.075 - 5.75,-35.25 - 4.85+0.15,sc*lonsc*40,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
%text(177.5 - 5.75,-35.25 - 0.005 - 4.85+0.15,'Holden et al. [2015]','fontsize',13.5)

text(173,-34.5 - 0.005 - 4.85,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(173 - 0.05,-34.5 - 4.95,7^2,-1,'filled')
text(173,-34.5 - 0.005 - 4.95,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(173 - 0.061,-34.5 - 5.05,10.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.05,['Pickle [2019+]'],'fontsize',13.5)
scatter(173 - 0.072,-34.5 - 5.15,14^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.15,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(173 - 0.083,-34.5 - 5.25,17.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.25,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(173 - 0.094,-34.5 - 5.35,21^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.35,'Holden et al. [2015]','fontsize',13.5)
scatter(173 - 0.105,-34.5 - 5.45,24.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.45,'Lamb and Smith [2013]','fontsize',13.5)
scatter(173 - 0.116,-34.5 - 5.55,28^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.55,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(173 - 0.127,-34.5 - 5.65,31.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.65,'Houlie & Stern [2017]','fontsize',13.5)
scatter(173 - 0.138,-34.5 - 5.75,35^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.75,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(173 - 0.149,-34.5 - 5.85,38.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.85,['Denys et al. [2020]'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([171.25 177.25])
ylim([-42 -39.25])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.019;
ax_width = outerpos(3) - 0.0325;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned'],'fontsize',25)
%return
%print(22,['NZGNSS_northsouth_' date '_aligned.png'],'-dpng');
%return

%%
sc = 7.5e-3;
figure(23); clf; hold on; box on;
set(gcf,'position',[0 720 1082 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_c.hamling22.lons,gps_c.hamling22.lats,7^2,gps_c.hamling22.vu,'filled')
scatter(gps_c.pickle19.lons,gps_c.pickle19.lats,10.5^2,gps_c.pickle19.vu,'linewidth',1.25)
scatter(gps_c.midas.lons,gps_c.midas.lats,14^2,gps_c.midas.vu,'linewidth',1.25)
scatter(gps_c.measures.lons,gps_c.measures.lats,17.5^2,gps_c.measures.vu,'linewidth',1.25)
scatter(gps_c.holden15.lons,gps_c.holden15.lats,21^2,gps_c.holden15.vu,'linewidth',1.25)
scatter(gps_c.lamb13.lons,gps_c.lamb13.lats,24.5^2,gps_c.lamb13.vu,'linewidth',1.25)
scatter(gps_c.beavan12.lons,gps_c.beavan12.lats,28^2,gps_c.beavan12.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,31.5^2,gps_c.houlie17.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,35^2,gps_c.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_c.denys20.lons,gps_c.denys20.lats,38.5^2,gps_c.denys20.vu,'linewidth',1.25)
scatter(gps_c.ulr5.lons,gps_c.ulr5.lats,42^2,gps_c.ulr5.vu,'linewidth',1.25)

% quiver(gps_c.lamb13.lons,gps_c.lamb13.lats,sc*lonsc*gps_c.lamb13.ve,sc*gps_c.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% thisellipse = ellipse(sc*lonsc*gps_c.lamb13.se,sc*gps_c.lamb13.sn,0*gps_c.lamb13.se,gps_c.lamb13.lons + sc*lonsc*gps_c.lamb13.ve,gps_c.lamb13.lats + sc*gps_c.lamb13.vn,color_lamb13,30);
% set(thisellipse,'linewidth',0.25)

quiver(gps_c.midas.lons,gps_c.midas.lats,sc*lonsc*gps_c.midas.ve,sc*gps_c.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.midas.se,sc*gps_c.midas.sn,0*gps_c.midas.se,gps_c.midas.lons + sc*lonsc*gps_c.midas.ve,gps_c.midas.lats + sc*gps_c.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.UNR.lons,gps_c.UNR.lats,sc*lonsc*gps_c.UNR.ve,sc*gps_c.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_c.UNR.se,sc*gps_c.UNR.sn,0*gps_c.UNR.se,gps_c.UNR.lons + sc*lonsc*gps_c.UNR.ve,gps_c.UNR.lats + sc*gps_c.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.measures.lons,gps_c.measures.lats,sc*lonsc*gps_c.measures.ve,sc*gps_c.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.measures.se,sc*gps_c.measures.sn,0*gps_c.measures.se,gps_c.measures.lons + sc*lonsc*gps_c.measures.ve,gps_c.measures.lats + sc*gps_c.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.pickle19.lons,gps_c.pickle19.lats,sc*lonsc*gps_c.pickle19.ve,sc*gps_c.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.pickle19.se,sc*gps_c.pickle19.sn,0*gps_c.pickle19.se,gps_c.pickle19.lons + sc*lonsc*gps_c.pickle19.ve,gps_c.pickle19.lats + sc*gps_c.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.measures.lons,gps_u.measures.lats,upper(gps_u.measures.names))

quiver(gps_c.beavan16.lons,gps_c.beavan16.lats,sc*lonsc*gps_c.beavan16.ve,sc*gps_c.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.beavan16.se,sc*gps_c.beavan16.sn,0*gps_c.beavan16.ve,gps_c.beavan16.lons + sc*lonsc*gps_c.beavan16.ve,gps_c.beavan16.lats + sc*gps_c.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_c.beavan16.lons(gps_c.beavan16.outlier),gps_c.beavan16.lats(gps_c.beavan16.outlier),upper(gps_c.beavan16.names(gps_c.beavan16.outlier)))

text(177.5 - 5.75,-34.5 - 0.005 - 4.85,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275 - 5.75,-34.68975 - 0.005 - 4.85+0.0375,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075 - 5.75,-34.625 - 4.85+0.025,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.625 - 0.005 - 4.85+0.025,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.075 - 5.75,-34.75 - 4.85+0.05,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.75 - 0.005 - 4.85+0.05,'Pickle [2019+]','fontsize',13.5)
quiver(177.075 - 5.75,-34.875 - 4.85+0.075,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.875 - 0.005 - 4.85+0.075,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.075 - 5.75,-35 - 4.85+0.1,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-35 - 0.005 - 4.85+0.1,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.075 - 5.75,-35.125 - 4.85+0.125,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-35.125 - 0.005 - 4.85+0.125,'UNR [Kreemer et al., 2014]','fontsize',13.5)
%quiver(177.075 - 5.75,-35.25 - 4.85+0.15,sc*lonsc*40,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
%text(177.5 - 5.75,-35.25 - 0.005 - 4.85+0.15,'Holden et al. [2015]','fontsize',13.5)

text(173,-34.5 - 0.005 - 4.85,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(173 - 0.05,-34.5 - 4.95,7^2,-1,'filled')
text(173,-34.5 - 0.005 - 4.95,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(173 - 0.061,-34.5 - 5.05,10.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.05,['Pickle [2019+]'],'fontsize',13.5)
scatter(173 - 0.072,-34.5 - 5.15,14^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.15,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(173 - 0.083,-34.5 - 5.25,17.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.25,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(173 - 0.094,-34.5 - 5.35,21^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.35,'Holden et al. [2015]','fontsize',13.5)
scatter(173 - 0.105,-34.5 - 5.45,24.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.45,'Lamb and Smith [2013]','fontsize',13.5)
scatter(173 - 0.116,-34.5 - 5.55,28^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.55,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(173 - 0.127,-34.5 - 5.65,31.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.65,'Houlie & Stern [2017]','fontsize',13.5)
scatter(173 - 0.138,-34.5 - 5.75,35^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.75,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(173 - 0.149,-34.5 - 5.85,38.5^2,-1,'linewidth',1.25)
text(173,-34.5 - 0.005 - 5.85,['Denys et al. [2020]'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([171.25 177.25])
ylim([-42 -39.25])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.019;
ax_width = outerpos(3) - 0.0325;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, outliers removed'],'fontsize',25)
%return
%print(23,['NZGNSS_northsouth_' date '_aligned_culled.png'],'-dpng');
%return

%%
sc = 7.5e-3;
figure(24); clf; hold on; box on;
set(gcf,'position',[0 720 1082 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),9^2,gps_m.vu(~isnan(gps_m.vu)),'filled','markeredgecolor','none')
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),14^2,gps_m.vu(~isnan(gps_m.vu)) - gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),19.5^2,gps_m.vu(~isnan(gps_m.vu)) + gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
gpsplot = quiver(gps_m.lons,gps_m.lats,sc/cosd(-41)*gps_m.ve,sc*gps_m.vn,0,'color','k','linewidth',1.25,'maxheadsize',1e-3);
gpsellipse = ellipse(sc/cosd(-41)*gps_m.se,sc*gps_m.sn,0*gps_m.sn,gps_m.lons + sc/cosd(-41)*gps_m.ve,gps_m.lats + sc*gps_m.vn,'k',30);
set(gpsellipse,'linewidth',0.25)

text(177.5 - 5.75,-34.5 - 0.005 - 4.85,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
%text(177.275 - 5.75,-34.68975 - 0.005 - 4.85+0.0375,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075 - 5.75,-34.625 - 4.85+0.025,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.75,-34.625 - 0.005 - 4.85+0.025,'40 mm/yr','fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([171.25 177.25])
ylim([-42 -39.25])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.019;
ax_width = outerpos(3) - 0.0325;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, outliers removed, averaged'],'fontsize',25)
%return
%print(24,['NZGNSS_northsouth_' date '_averaged.png'],'-dpng');
%return
%%
sc = 7.5e-3;
figure(31); clf; hold on; box on;
set(gcf,'position',[0 720 1234 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_u.hamling22.lons,gps_u.hamling22.lats,7^2,gps_u.hamling22.vu,'filled')
scatter(gps_u.pickle19.lons,gps_u.pickle19.lats,10.5^2,gps_u.pickle19.vu,'linewidth',1.25)
scatter(gps_u.midas.lons,gps_u.midas.lats,14^2,gps_u.midas.vu,'linewidth',1.25)
scatter(gps_u.measures.lons,gps_u.measures.lats,17.5^2,gps_u.measures.vu,'linewidth',1.25)
scatter(gps_u.holden15.lons,gps_u.holden15.lats,21^2,gps_u.holden15.vu,'linewidth',1.25)
scatter(gps_u.lamb13.lons,gps_u.lamb13.lats,24.5^2,gps_u.lamb13.vu,'linewidth',1.25)
scatter(gps_u.beavan12.lons,gps_u.beavan12.lats,28^2,gps_u.beavan12.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,31.5^2,gps_u.houlie17.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,35^2,gps_u.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_u.denys20.lons,gps_u.denys20.lats,38.5^2,gps_u.denys20.vu,'linewidth',1.25)
scatter(gps_u.ulr5.lons,gps_u.ulr5.lats,42^2,gps_u.ulr5.vu,'linewidth',1.25)
scatter(gps_u.beavan10.lons,gps_u.beavan10.lats,45.5^2,gps_u.beavan10.vu,'linewidth',1.25)

% quiver(gps_u.lamb13.lons,gps_u.lamb13.lats,sc*lonsc*gps_u.lamb13.ve,sc*gps_u.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% thisellipse = ellipse(sc*lonsc*gps_u.lamb13.se,sc*gps_u.lamb13.sn,0*gps_u.lamb13.se,gps_u.lamb13.lons + sc*lonsc*gps_u.lamb13.ve,gps_u.lamb13.lats + sc*gps_u.lamb13.vn,color_lamb13,30);
% set(thisellipse,'linewidth',0.25)

quiver(gps_u.midas.lons,gps_u.midas.lats,sc*lonsc*gps_u.midas.ve,sc*gps_u.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.midas.se,sc*gps_u.midas.sn,0*gps_u.midas.se,gps_u.midas.lons + sc*lonsc*gps_u.midas.ve,gps_u.midas.lats + sc*gps_u.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.UNR.lons,gps_u.UNR.lats,sc*lonsc*gps_u.UNR.ve,sc*gps_u.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_u.UNR.se,sc*gps_u.UNR.sn,0*gps_u.UNR.se,gps_u.UNR.lons + sc*lonsc*gps_u.UNR.ve,gps_u.UNR.lats + sc*gps_u.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.measures.lons,gps_u.measures.lats,sc*lonsc*gps_u.measures.ve,sc*gps_u.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.measures.se,sc*gps_u.measures.sn,0*gps_u.measures.se,gps_u.measures.lons + sc*lonsc*gps_u.measures.ve,gps_u.measures.lats + sc*gps_u.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.pickle19.lons,gps_u.pickle19.lats,sc*lonsc*gps_u.pickle19.ve,sc*gps_u.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.pickle19.se,sc*gps_u.pickle19.sn,0*gps_u.pickle19.se,gps_u.pickle19.lons + sc*lonsc*gps_u.pickle19.ve,gps_u.pickle19.lats + sc*gps_u.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.pickle19.lons,gps_u.pickle19.lats,upper(gps_u.pickle19.names))

quiver(gps_u.beavan16.lons,gps_u.beavan16.lats,sc*lonsc*gps_u.beavan16.ve,sc*gps_u.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.beavan16.se,sc*gps_u.beavan16.sn,0*gps_u.beavan16.ve,gps_u.beavan16.lons + sc*lonsc*gps_u.beavan16.ve,gps_u.beavan16.lats + sc*gps_u.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_u.beavan16.lons(gps_u.beavan16.outlier),gps_u.beavan16.lats(gps_u.beavan16.outlier),upper(gps_u.beavan16.names(gps_u.beavan16.outlier)))

text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(169.525,-41.955,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.075,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.075,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.15,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.15,'Pickle [2019+]','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.225,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.225,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.3,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.3,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.375,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.375,'UNR [Kreemer et al., 2014]','fontsize',13.5)
%quiver(177.075 - 3.25 - 4.5,-34.5 - 7.3375 - 0.45,sc*lonsc*40,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
%text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.45,'Holden et al. [2015]','fontsize',13.5)

text(177.5 - 3.25  - 4.5,-34.5 - 0.005 - 7.3375 - 0.475,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(177.465 - 3.25  - 4.5,-42.3875,7^2,-1,'filled')
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(177.457 - 3.25  - 4.5,-42.3875 - 0.075,10.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.075,['Pickle [2019+]'],'fontsize',13.5)
scatter(177.449 - 3.25  - 4.5,-42.3875 - 0.15,14^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.15,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(177.441 - 3.25  - 4.5,-42.3875 - 0.225,17.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.225,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(177.433 - 3.25  - 4.5,-42.3875 - 0.3,21^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.3,'Holden et al. [2015]','fontsize',13.5)
scatter(177.425 - 3.25  - 4.5,-42.3875 - 0.375,24.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.375,'Lamb and Smith [2013]','fontsize',13.5)
scatter(177.417 - 3.25  - 4.5,-42.3875 - 0.45,28^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.45,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(177.409 - 3.25  - 4.5,-42.3875 - 0.525,31.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.525,'Houlie & Stern [2017]','fontsize',13.5)
scatter(177.401 - 3.25  - 4.5,-42.3875 - 0.6,35^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.6,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(177.393 - 3.25  - 4.5,-42.3875 - 0.675,38.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.675,['Denys et al. [2020]'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([169.25 174.25])
ylim([-43.75 -41.75])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0195;
ax_width = outerpos(3) - 0.0275;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities'],'fontsize',25)
%return
%print(31,['NZGNSS_midsouth_' date '.png'],'-dpng');
%return

%%
sc = 7.5e-3;
figure(32); clf; hold on; box on;
set(gcf,'position',[0 720 1234 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_a.hamling22.lons,gps_a.hamling22.lats,7^2,gps_a.hamling22.vu,'filled')
scatter(gps_a.pickle19.lons,gps_a.pickle19.lats,10.5^2,gps_a.pickle19.vu,'linewidth',1.25)
scatter(gps_a.midas.lons,gps_a.midas.lats,14^2,gps_a.midas.vu,'linewidth',1.25)
scatter(gps_a.measures.lons,gps_a.measures.lats,17.5^2,gps_a.measures.vu,'linewidth',1.25)
scatter(gps_a.holden15.lons,gps_a.holden15.lats,21^2,gps_a.holden15.vu,'linewidth',1.25)
scatter(gps_a.lamb13.lons,gps_a.lamb13.lats,24.5^2,gps_a.lamb13.vu,'linewidth',1.25)
scatter(gps_a.beavan12.lons,gps_a.beavan12.lats,28^2,gps_a.beavan12.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,31.5^2,gps_a.houlie17.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,35^2,gps_a.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_a.denys20.lons,gps_a.denys20.lats,38.5^2,gps_a.denys20.vu,'linewidth',1.25)
scatter(gps_a.ulr5.lons,gps_a.ulr5.lats,42^2,gps_a.ulr5.vu,'linewidth',1.25)
scatter(gps_a.beavan10.lons,gps_a.beavan10.lats,45.5^2,gps_a.beavan10.vu,'linewidth',1.25)

quiver(gps_a.midas.lons,gps_a.midas.lats,sc*lonsc*gps_a.midas.ve,sc*gps_a.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.midas.se,sc*gps_a.midas.sn,0*gps_a.midas.se,gps_a.midas.lons + sc*lonsc*gps_a.midas.ve,gps_a.midas.lats + sc*gps_a.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.midas.lons(gps_a.midas.outlier),gps_a.midas.lats(gps_a.midas.outlier),upper(gps_a.midas.names(gps_a.midas.outlier)))

quiver(gps_a.UNR.lons,gps_a.UNR.lats,sc*lonsc*gps_a.UNR.ve,sc*gps_a.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_a.UNR.se,sc*gps_a.UNR.sn,0*gps_a.UNR.se,gps_a.UNR.lons + sc*lonsc*gps_a.UNR.ve,gps_a.UNR.lats + sc*gps_a.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.UNR.lons(gps_a.UNR.outlier),gps_a.UNR.lats(gps_a.UNR.outlier),upper(gps_a.UNR.names(gps_a.UNR.outlier)))

quiver(gps_a.measures.lons,gps_a.measures.lats,sc*lonsc*gps_a.measures.ve,sc*gps_a.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.measures.se,sc*gps_a.measures.sn,0*gps_a.measures.se,gps_a.measures.lons + sc*lonsc*gps_a.measures.ve,gps_a.measures.lats + sc*gps_a.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.measures.lons(gps_a.measures.outlier),gps_a.measures.lats(gps_a.measures.outlier),upper(gps_a.measures.names(gps_a.measures.outlier)))

quiver(gps_a.pickle19.lons,gps_a.pickle19.lats,sc*lonsc*gps_a.pickle19.ve,sc*gps_a.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.pickle19.se,sc*gps_a.pickle19.sn,0*gps_a.pickle19.se,gps_a.pickle19.lons + sc*lonsc*gps_a.pickle19.ve,gps_a.pickle19.lats + sc*gps_a.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.pickle19.lons(gps_a.pickle19.outlier),gps_a.pickle19.lats(gps_a.pickle19.outlier),upper(gps_a.pickle19.names(gps_a.pickle19.outlier)))

quiver(gps_a.beavan16.lons,gps_a.beavan16.lats,sc*lonsc*gps_a.beavan16.ve,sc*gps_a.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.beavan16.se,sc*gps_a.beavan16.sn,0*gps_a.beavan16.ve,gps_a.beavan16.lons + sc*lonsc*gps_a.beavan16.ve,gps_a.beavan16.lats + sc*gps_a.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.beavan16.lons(gps_a.beavan16.outlier),gps_a.beavan16.lats(gps_a.beavan16.outlier),upper(gps_a.beavan16.names(gps_a.beavan16.outlier)))

text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(169.525,-41.955,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.075,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.075,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.15,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.15,'Pickle [2019+]','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.225,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.225,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.3,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.3,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.375,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.375,'UNR [Kreemer et al., 2014]','fontsize',13.5)
%quiver(177.075 - 3.25 - 4.5,-34.5 - 7.3375 - 0.45,sc*lonsc*40,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
%text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.45,'Holden et al. [2015]','fontsize',13.5)

text(177.5 - 3.25  - 4.5,-34.5 - 0.005 - 7.3375 - 0.475,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(177.465 - 3.25  - 4.5,-42.3875,7^2,-1,'filled')
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(177.457 - 3.25  - 4.5,-42.3875 - 0.075,10.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.075,['Pickle [2019+]'],'fontsize',13.5)
scatter(177.449 - 3.25  - 4.5,-42.3875 - 0.15,14^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.15,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(177.441 - 3.25  - 4.5,-42.3875 - 0.225,17.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.225,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(177.433 - 3.25  - 4.5,-42.3875 - 0.3,21^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.3,'Holden et al. [2015]','fontsize',13.5)
scatter(177.425 - 3.25  - 4.5,-42.3875 - 0.375,24.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.375,'Lamb and Smith [2013]','fontsize',13.5)
scatter(177.417 - 3.25  - 4.5,-42.3875 - 0.45,28^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.45,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(177.409 - 3.25  - 4.5,-42.3875 - 0.525,31.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.525,'Houlie & Stern [2017]','fontsize',13.5)
scatter(177.401 - 3.25  - 4.5,-42.3875 - 0.6,35^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.6,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(177.393 - 3.25  - 4.5,-42.3875 - 0.675,38.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.675,['Denys et al. [2020]'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([169.25 174.25])
ylim([-43.75 -41.75])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0195;
ax_width = outerpos(3) - 0.0275;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned'],'fontsize',25)
%return
%print(32,['NZGNSS_midsouth_' date '_aligned.png'],'-dpng');
%return
%%
sc = 7.5e-3;
figure(33); clf; hold on; box on;
set(gcf,'position',[0 720 1234 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_c.hamling22.lons,gps_c.hamling22.lats,7^2,gps_c.hamling22.vu,'filled')
scatter(gps_c.pickle19.lons,gps_c.pickle19.lats,10.5^2,gps_c.pickle19.vu,'linewidth',1.25)
scatter(gps_c.midas.lons,gps_c.midas.lats,14^2,gps_c.midas.vu,'linewidth',1.25)
scatter(gps_c.measures.lons,gps_c.measures.lats,17.5^2,gps_c.measures.vu,'linewidth',1.25)
scatter(gps_c.holden15.lons,gps_c.holden15.lats,21^2,gps_c.holden15.vu,'linewidth',1.25)
scatter(gps_c.lamb13.lons,gps_c.lamb13.lats,24.5^2,gps_c.lamb13.vu,'linewidth',1.25)
scatter(gps_c.beavan12.lons,gps_c.beavan12.lats,28^2,gps_c.beavan12.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,31.5^2,gps_c.houlie17.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,35^2,gps_c.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_c.denys20.lons,gps_c.denys20.lats,38.5^2,gps_c.denys20.vu,'linewidth',1.25)
scatter(gps_c.ulr5.lons,gps_c.ulr5.lats,42^2,gps_c.ulr5.vu,'linewidth',1.25)
scatter(gps_c.beavan10.lons,gps_c.beavan10.lats,45.5^2,gps_c.beavan10.vu,'linewidth',1.25)

% quiver(gps_c.lamb13.lons,gps_c.lamb13.lats,sc*lonsc*gps_c.lamb13.ve,sc*gps_c.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% thisellipse = ellipse(sc*lonsc*gps_c.lamb13.se,sc*gps_c.lamb13.sn,0*gps_c.lamb13.se,gps_c.lamb13.lons + sc*lonsc*gps_c.lamb13.ve,gps_c.lamb13.lats + sc*gps_c.lamb13.vn,color_lamb13,30);
% set(thisellipse,'linewidth',0.25)

quiver(gps_c.midas.lons,gps_c.midas.lats,sc*lonsc*gps_c.midas.ve,sc*gps_c.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.midas.se,sc*gps_c.midas.sn,0*gps_c.midas.se,gps_c.midas.lons + sc*lonsc*gps_c.midas.ve,gps_c.midas.lats + sc*gps_c.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.UNR.lons,gps_c.UNR.lats,sc*lonsc*gps_c.UNR.ve,sc*gps_c.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_c.UNR.se,sc*gps_c.UNR.sn,0*gps_c.UNR.se,gps_c.UNR.lons + sc*lonsc*gps_c.UNR.ve,gps_c.UNR.lats + sc*gps_c.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.measures.lons,gps_c.measures.lats,sc*lonsc*gps_c.measures.ve,sc*gps_c.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.measures.se,sc*gps_c.measures.sn,0*gps_c.measures.se,gps_c.measures.lons + sc*lonsc*gps_c.measures.ve,gps_c.measures.lats + sc*gps_c.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.pickle19.lons,gps_c.pickle19.lats,sc*lonsc*gps_c.pickle19.ve,sc*gps_c.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.pickle19.se,sc*gps_c.pickle19.sn,0*gps_c.pickle19.se,gps_c.pickle19.lons + sc*lonsc*gps_c.pickle19.ve,gps_c.pickle19.lats + sc*gps_c.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_c.pickle19.lons,gps_c.pickle19.lats,upper(gps_c.pickle19.names))

quiver(gps_c.beavan16.lons,gps_c.beavan16.lats,sc*lonsc*gps_c.beavan16.ve,sc*gps_c.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.beavan16.se,sc*gps_c.beavan16.sn,0*gps_c.beavan16.ve,gps_c.beavan16.lons + sc*lonsc*gps_c.beavan16.ve,gps_c.beavan16.lats + sc*gps_c.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_c.beavan16.lons(gps_c.beavan16.outlier),gps_c.beavan16.lats(gps_c.beavan16.outlier),upper(gps_c.beavan16.names(gps_c.beavan16.outlier)))

text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(169.525,-41.955,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.075,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.075,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.15,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.15,'Pickle [2019+]','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.225,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.225,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.3,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.3,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.375,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.375,'UNR [Kreemer et al., 2014]','fontsize',13.5)
%quiver(177.075 - 3.25 - 4.5,-34.5 - 7.3375 - 0.45,sc*lonsc*40,sc*0,0,'color',color_holden15,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
%text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.45,'Holden et al. [2015]','fontsize',13.5)

text(177.5 - 3.25  - 4.5,-34.5 - 0.005 - 7.3375 - 0.475,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(177.465 - 3.25  - 4.5,-42.3875,7^2,-1,'filled')
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(177.457 - 3.25  - 4.5,-42.3875 - 0.075,10.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.075,['Pickle [2019+]'],'fontsize',13.5)
scatter(177.449 - 3.25  - 4.5,-42.3875 - 0.15,14^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.15,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(177.441 - 3.25  - 4.5,-42.3875 - 0.225,17.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.225,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(177.433 - 3.25  - 4.5,-42.3875 - 0.3,21^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.3,'Holden et al. [2015]','fontsize',13.5)
scatter(177.425 - 3.25  - 4.5,-42.3875 - 0.375,24.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.375,'Lamb and Smith [2013]','fontsize',13.5)
scatter(177.417 - 3.25  - 4.5,-42.3875 - 0.45,28^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.45,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(177.409 - 3.25  - 4.5,-42.3875 - 0.525,31.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.525,'Houlie & Stern [2017]','fontsize',13.5)
scatter(177.401 - 3.25  - 4.5,-42.3875 - 0.6,35^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.6,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(177.393 - 3.25  - 4.5,-42.3875 - 0.675,38.5^2,-1,'linewidth',1.25)
text(177.5 - 3.25  - 4.5,-42.3875 - 0.005 - 0.675,['Denys et al. [2020]'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([169.25 174.25])
ylim([-43.75 -41.75])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0195;
ax_width = outerpos(3) - 0.0275;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, outliers removed'],'fontsize',25)
%return
%print(33,['NZGNSS_midsouth_' date '_aligned_culled.png'],'-dpng');

%%
sc = 7.5e-3;
figure(34); clf; hold on; box on;
set(gcf,'position',[0 720 1234 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),9^2,gps_m.vu(~isnan(gps_m.vu)),'filled','markeredgecolor','none')
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),14^2,gps_m.vu(~isnan(gps_m.vu)) - gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),19.5^2,gps_m.vu(~isnan(gps_m.vu)) + gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
gpsplot = quiver(gps_m.lons,gps_m.lats,sc/cosd(-41)*gps_m.ve,sc*gps_m.vn,0,'color','k','linewidth',1.25,'maxheadsize',1e-3);
gpsellipse = ellipse(sc/cosd(-41)*gps_m.se,sc*gps_m.sn,0*gps_m.sn,gps_m.lons + sc/cosd(-41)*gps_m.ve,gps_m.lats + sc*gps_m.vn,'k',30);
set(gpsellipse,'linewidth',0.25)

text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
%text(169.525,-41.955,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.085 - 3.25 - 4.5,-34.5 - 7.3375 - 0.075,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 3.25 - 4.5,-34.5 - 0.005 - 7.3375 - 0.075,'40 mm/yr','fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([169.25 174.25])
ylim([-43.75 -41.75])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.0195;
ax_width = outerpos(3) - 0.0275;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, outliers removed, averaged'],'fontsize',25)
%return
%print(34,['NZGNSS_midsouth_' date '_averaged.png'],'-dpng');
%return
%%
sc = 7.5e-3;
figure(41); clf; hold on; box on;
set(gcf,'position',[0 720 1066 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_u.hamling22.lons,gps_u.hamling22.lats,7^2,gps_u.hamling22.vu,'filled')
scatter(gps_u.pickle19.lons,gps_u.pickle19.lats,10.5^2,gps_u.pickle19.vu,'linewidth',1.25)
scatter(gps_u.midas.lons,gps_u.midas.lats,14^2,gps_u.midas.vu,'linewidth',1.25)
scatter(gps_u.measures.lons,gps_u.measures.lats,17.5^2,gps_u.measures.vu,'linewidth',1.25)
scatter(gps_u.holden15.lons,gps_u.holden15.lats,21^2,gps_u.holden15.vu,'linewidth',1.25)
scatter(gps_u.lamb13.lons,gps_u.lamb13.lats,24.5^2,gps_u.lamb13.vu,'linewidth',1.25)
scatter(gps_u.beavan12.lons,gps_u.beavan12.lats,28^2,gps_u.beavan12.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,31.5^2,gps_u.houlie17.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,35^2,gps_u.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_u.denys20.lons,gps_u.denys20.lats,38.5^2,gps_u.denys20.vu,'linewidth',1.25)
scatter(gps_u.ulr5.lons,gps_u.ulr5.lats,42^2,gps_u.ulr5.vu,'linewidth',1.25)
scatter(gps_u.beavan10.lons,gps_u.beavan10.lats,45.5^2,gps_u.beavan10.vu,'linewidth',1.25)

% quiver(gps_u.lamb13.lons,gps_u.lamb13.lats,sc*lonsc*gps_u.lamb13.ve,sc*gps_u.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% thisellipse = ellipse(sc*lonsc*gps_u.lamb13.se,sc*gps_u.lamb13.sn,0*gps_u.lamb13.se,gps_u.lamb13.lons + sc*lonsc*gps_u.lamb13.ve,gps_u.lamb13.lats + sc*gps_u.lamb13.vn,color_lamb13,30);
% set(thisellipse,'linewidth',0.25)

quiver(gps_u.midas.lons,gps_u.midas.lats,sc*lonsc*gps_u.midas.ve,sc*gps_u.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.midas.se,sc*gps_u.midas.sn,0*gps_u.midas.se,gps_u.midas.lons + sc*lonsc*gps_u.midas.ve,gps_u.midas.lats + sc*gps_u.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.UNR.lons,gps_u.UNR.lats,sc*lonsc*gps_u.UNR.ve,sc*gps_u.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_u.UNR.se,sc*gps_u.UNR.sn,0*gps_u.UNR.se,gps_u.UNR.lons + sc*lonsc*gps_u.UNR.ve,gps_u.UNR.lats + sc*gps_u.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.measures.lons,gps_u.measures.lats,sc*lonsc*gps_u.measures.ve,sc*gps_u.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.measures.se,sc*gps_u.measures.sn,0*gps_u.measures.se,gps_u.measures.lons + sc*lonsc*gps_u.measures.ve,gps_u.measures.lats + sc*gps_u.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.pickle19.lons,gps_u.pickle19.lats,sc*lonsc*gps_u.pickle19.ve,sc*gps_u.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.pickle19.se,sc*gps_u.pickle19.sn,0*gps_u.pickle19.se,gps_u.pickle19.lons + sc*lonsc*gps_u.pickle19.ve,gps_u.pickle19.lats + sc*gps_u.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.beavan16.lons,gps_u.beavan16.lats,sc*lonsc*gps_u.beavan16.ve,sc*gps_u.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.beavan16.se,sc*gps_u.beavan16.sn,0*gps_u.beavan16.ve,gps_u.beavan16.lons + sc*lonsc*gps_u.beavan16.ve,gps_u.beavan16.lats + sc*gps_u.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.denys16.lons,gps_u.denys16.lats,sc*lonsc*gps_u.denys16.ve,sc*gps_u.denys16.vn,0,'color',color_denys16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.denys16.se,sc*gps_u.denys16.sn,0*gps_u.denys16.ve,gps_u.denys16.lons + sc*lonsc*gps_u.denys16.ve,gps_u.denys16.lats + sc*gps_u.denys16.vn,color_denys16,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_u.page18.lons,gps_u.page18.lats,sc*lonsc*gps_u.page18.ve,sc*gps_u.page18.vn,0,'color',color_page18,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_u.page18.se,sc*gps_u.page18.sn,0*gps_u.page18.ve,gps_u.page18.lons + sc*lonsc*gps_u.page18.ve,gps_u.page18.lats + sc*gps_u.page18.vn,color_page18,30);
set(thisellipse,'linewidth',0.25)


text(177.5 - 5.5,-34.5 - 0.006 - 10 + 0.5,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275 - 5.5,-34.68975 - 0.006 - 10 + 0.5,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075 - 5.5,-34.625 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-34.625 - 0.006 - 10 + 0.5,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.075 - 5.5,-34.75 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-34.75 - 0.006 - 10 + 0.5,'Pickle [2019+]','fontsize',13.5)
quiver(177.075 - 5.5,-34.875 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-34.875 - 0.006 - 10 + 0.5,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.075 - 5.5,-35 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35 - 0.006 - 10 + 0.5,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.075 - 5.5,-35.125 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35.125 - 0.006 - 10 + 0.5,'UNR [Kreemer et al., 2014]','fontsize',13.5)
quiver(177.075 - 5.5,-35.25 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_denys16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35.25 - 0.006 - 10 + 0.5,'Denys et al. [2016]','fontsize',13.5)
quiver(177.075 - 5.5,-35.375 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_page18,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35.375 - 0.006 - 10 + 0.5,'Page et al. [2018]','fontsize',13.5)

text(177.5 - 5.5,-34.75 - 0.006 - 10.8125 + 0.5,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(177.44 - 5.5,-34.875 - 10.8125 + 0.5,7^2,-1,'filled')
text(177.5 - 5.5,-34.875 - 0.006 - 10.8125 + 0.5,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(177.426 - 5.5,-35 - 10.8125 + 0.5,10.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35 - 0.006 - 10.8125 + 0.5,['Pickle [2019+]'],'fontsize',13.5)
scatter(177.412 - 5.5,-35.125 - 10.8125 + 0.5,14^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.125 - 0.006 - 10.8125 + 0.5,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(177.398 - 5.5,-35.25 - 10.8125 + 0.5,17.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.25 - 0.006 - 10.8125 + 0.5,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(177.384 - 5.5,-35.375 - 10.8125 + 0.5,21^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.375 - 0.006 - 10.8125 + 0.5,'Holden et al. [2015]','fontsize',13.5)
scatter(177.37 - 5.5,-35.5 - 10.8125 + 0.5,24.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.5 - 0.006 - 10.8125 + 0.5,'Lamb and Smith [2013]','fontsize',13.5)
scatter(177.356 - 5.5,-35.625 - 10.8125 + 0.5,28^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.625 - 0.006 - 10.8125 + 0.5,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(177.342 - 5.5,-35.75 - 10.8125 + 0.5,31.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.75 - 0.006 - 10.8125 + 0.5,'Houlie & Stern [2017]','fontsize',13.5)
scatter(177.328 - 5.5,-35.875 - 10.8125 + 0.5,35^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.875 - 0.006 - 10.8125 + 0.5,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(177.314 - 5.5,-36 - 10.8125 + 0.5,38.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-36 - 0.006 - 10.8125 + 0.5,['Denys et al. [2020]'],'fontsize',13.5)
scatter(177.3 - 5.5,-36.125 - 10.8125 + 0.5,42^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-36.125 - 0.006 - 10.8125 + 0.5,'ULR5 [2012]','fontsize',13.5)
scatter(177.296 - 5.5,-36.25 - 10.8125 + 0.5,42^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-36.25 - 0.006 - 10.8125 + 0.5,'Beavan et al. [2010]','fontsize',13.5)
%scatter(177.286 - 9.25,-36.25 - 0.8125,45.5^2,-1,'linewidth',1.25)
%text(177.5 - 9.25,-36.25 - 0.006 - 0.8125,['Vardic et al. [2021]'],'fontsize',13.5)
%scatter(177.272 - 9.25,-36.375 - 0.8125,49^2,-1,'linewidth',1.25)
%text(177.5 - 9.25,-36.375 - 0.006 - 0.8125,['ITRF2014 VLBI'],'fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([166 173.25])
ylim([-46.75 -43.375])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.01875;
ax_width = outerpos(3) - 0.0325;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities'],'fontsize',25)
%return
%print(41,['NZGNSS_lowersouth_' date '.png'],'-dpng');
%return
%%
sc = 7.5e-3;
figure(42); clf; hold on; box on;
set(gcf,'position',[0 720 1066 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_a.hamling22.lons,gps_a.hamling22.lats,7^2,gps_a.hamling22.vu,'filled')
scatter(gps_a.pickle19.lons,gps_a.pickle19.lats,10.5^2,gps_a.pickle19.vu,'linewidth',1.25)
scatter(gps_a.midas.lons,gps_a.midas.lats,14^2,gps_a.midas.vu,'linewidth',1.25)
scatter(gps_a.measures.lons,gps_a.measures.lats,17.5^2,gps_a.measures.vu,'linewidth',1.25)
scatter(gps_a.holden15.lons,gps_a.holden15.lats,21^2,gps_a.holden15.vu,'linewidth',1.25)
scatter(gps_a.lamb13.lons,gps_a.lamb13.lats,24.5^2,gps_a.lamb13.vu,'linewidth',1.25)
scatter(gps_a.beavan12.lons,gps_a.beavan12.lats,28^2,gps_a.beavan12.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,31.5^2,gps_a.houlie17.vu,'linewidth',1.25)
scatter(gps_a.houlie17.lons,gps_a.houlie17.lats,35^2,gps_a.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_a.denys20.lons,gps_a.denys20.lats,38.5^2,gps_a.denys20.vu,'linewidth',1.25)
scatter(gps_a.ulr5.lons,gps_a.ulr5.lats,42^2,gps_a.ulr5.vu,'linewidth',1.25)
scatter(gps_a.beavan10.lons,gps_a.beavan10.lats,45.5^2,gps_a.beavan10.vu,'linewidth',1.25)

quiver(gps_a.midas.lons,gps_a.midas.lats,sc*lonsc*gps_a.midas.ve,sc*gps_a.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.midas.se,sc*gps_a.midas.sn,0*gps_a.midas.se,gps_a.midas.lons + sc*lonsc*gps_a.midas.ve,gps_a.midas.lats + sc*gps_a.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.midas.lons(gps_a.midas.outlier),gps_a.midas.lats(gps_a.midas.outlier),upper(gps_a.midas.names(gps_a.midas.outlier)))

quiver(gps_a.UNR.lons,gps_a.UNR.lats,sc*lonsc*gps_a.UNR.ve,sc*gps_a.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_a.UNR.se,sc*gps_a.UNR.sn,0*gps_a.UNR.se,gps_a.UNR.lons + sc*lonsc*gps_a.UNR.ve,gps_a.UNR.lats + sc*gps_a.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.UNR.lons(gps_a.UNR.outlier),gps_a.UNR.lats(gps_a.UNR.outlier),upper(gps_a.UNR.names(gps_a.UNR.outlier)))

quiver(gps_a.measures.lons,gps_a.measures.lats,sc*lonsc*gps_a.measures.ve,sc*gps_a.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.measures.se,sc*gps_a.measures.sn,0*gps_a.measures.se,gps_a.measures.lons + sc*lonsc*gps_a.measures.ve,gps_a.measures.lats + sc*gps_a.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.measures.lons(gps_a.measures.outlier),gps_a.measures.lats(gps_a.measures.outlier),upper(gps_a.measures.names(gps_a.measures.outlier)))

quiver(gps_a.pickle19.lons,gps_a.pickle19.lats,sc*lonsc*gps_a.pickle19.ve,sc*gps_a.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.pickle19.se,sc*gps_a.pickle19.sn,0*gps_a.pickle19.se,gps_a.pickle19.lons + sc*lonsc*gps_a.pickle19.ve,gps_a.pickle19.lats + sc*gps_a.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.pickle19.lons(gps_a.pickle19.outlier),gps_a.pickle19.lats(gps_a.pickle19.outlier),upper(gps_a.pickle19.names(gps_a.pickle19.outlier)))

quiver(gps_a.beavan16.lons,gps_a.beavan16.lats,sc*lonsc*gps_a.beavan16.ve,sc*gps_a.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.beavan16.se,sc*gps_a.beavan16.sn,0*gps_a.beavan16.ve,gps_a.beavan16.lons + sc*lonsc*gps_a.beavan16.ve,gps_a.beavan16.lats + sc*gps_a.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.beavan16.lons(gps_a.beavan16.outlier),gps_a.beavan16.lats(gps_a.beavan16.outlier),upper(gps_a.beavan16.names(gps_a.beavan16.outlier)))

quiver(gps_a.denys16.lons,gps_a.denys16.lats,sc*lonsc*gps_a.denys16.ve,sc*gps_a.denys16.vn,0,'color',color_denys16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.denys16.se,sc*gps_a.denys16.sn,0*gps_a.denys16.ve,gps_a.denys16.lons + sc*lonsc*gps_a.denys16.ve,gps_a.denys16.lats + sc*gps_a.denys16.vn,color_denys16,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.denys16.lons(gps_a.denys16.outlier),gps_a.denys16.lats(gps_a.denys16.outlier),upper(gps_a.denys16.names(gps_a.denys16.outlier)))

quiver(gps_a.page18.lons,gps_a.page18.lats,sc*lonsc*gps_a.page18.ve,sc*gps_a.page18.vn,0,'color',color_page18,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_a.page18.se,sc*gps_a.page18.sn,0*gps_a.page18.ve,gps_a.page18.lons + sc*lonsc*gps_a.page18.ve,gps_a.page18.lats + sc*gps_a.page18.vn,color_page18,30);
set(thisellipse,'linewidth',0.25)
%text(gps_a.page18.lons(gps_a.page18.outlier),gps_a.page18.lats(gps_a.page18.outlier),upper(gps_a.page18.names(gps_a.page18.outlier)))

text(177.5 - 5.5,-34.5 - 0.006 - 10 + 0.5,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275 - 5.5,-34.68975 - 0.006 - 10 + 0.5,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075 - 5.5,-34.625 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-34.625 - 0.006 - 10 + 0.5,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.075 - 5.5,-34.75 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-34.75 - 0.006 - 10 + 0.5,'Pickle [2019+]','fontsize',13.5)
quiver(177.075 - 5.5,-34.875 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-34.875 - 0.006 - 10 + 0.5,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.075 - 5.5,-35 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35 - 0.006 - 10 + 0.5,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.075 - 5.5,-35.125 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35.125 - 0.006 - 10 + 0.5,'UNR [Kreemer et al., 2014]','fontsize',13.5)
quiver(177.075 - 5.5,-35.25 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_denys16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35.25 - 0.006 - 10 + 0.5,'Denys et al. [2016]','fontsize',13.5)
quiver(177.075 - 5.5,-35.375 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_page18,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35.375 - 0.006 - 10 + 0.5,'Page et al. [2018]','fontsize',13.5)

text(177.5 - 5.5,-34.75 - 0.006 - 10.8125 + 0.5,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(177.44 - 5.5,-34.875 - 10.8125 + 0.5,7^2,-1,'filled')
text(177.5 - 5.5,-34.875 - 0.006 - 10.8125 + 0.5,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(177.426 - 5.5,-35 - 10.8125 + 0.5,10.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35 - 0.006 - 10.8125 + 0.5,['Pickle [2019+]'],'fontsize',13.5)
scatter(177.412 - 5.5,-35.125 - 10.8125 + 0.5,14^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.125 - 0.006 - 10.8125 + 0.5,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(177.398 - 5.5,-35.25 - 10.8125 + 0.5,17.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.25 - 0.006 - 10.8125 + 0.5,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(177.384 - 5.5,-35.375 - 10.8125 + 0.5,21^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.375 - 0.006 - 10.8125 + 0.5,'Holden et al. [2015]','fontsize',13.5)
scatter(177.37 - 5.5,-35.5 - 10.8125 + 0.5,24.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.5 - 0.006 - 10.8125 + 0.5,'Lamb and Smith [2013]','fontsize',13.5)
scatter(177.356 - 5.5,-35.625 - 10.8125 + 0.5,28^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.625 - 0.006 - 10.8125 + 0.5,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(177.342 - 5.5,-35.75 - 10.8125 + 0.5,31.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.75 - 0.006 - 10.8125 + 0.5,'Houlie & Stern [2017]','fontsize',13.5)
scatter(177.328 - 5.5,-35.875 - 10.8125 + 0.5,35^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.875 - 0.006 - 10.8125 + 0.5,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(177.314 - 5.5,-36 - 10.8125 + 0.5,38.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-36 - 0.006 - 10.8125 + 0.5,['Denys et al. [2020]'],'fontsize',13.5)
scatter(177.3 - 5.5,-36.125 - 10.8125 + 0.5,42^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-36.125 - 0.006 - 10.8125 + 0.5,'ULR5 [2012]','fontsize',13.5)
scatter(177.296 - 5.5,-36.25 - 10.8125 + 0.5,42^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-36.25 - 0.006 - 10.8125 + 0.5,'Beavan et al. [2010]','fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([166 173.25])
ylim([-46.75 -43.375])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.01875;
ax_width = outerpos(3) - 0.0325;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned'],'fontsize',25)
%return
%print(42,['NZGNSS_lowersouth_' date '_aligned.png'],'-dpng');

%%
sc = 7.5e-3;
figure(43); clf; hold on; box on;
set(gcf,'position',[0 720 1066 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_c.hamling22.lons,gps_c.hamling22.lats,7^2,gps_c.hamling22.vu,'filled')
scatter(gps_c.pickle19.lons,gps_c.pickle19.lats,10.5^2,gps_c.pickle19.vu,'linewidth',1.25)
scatter(gps_c.midas.lons,gps_c.midas.lats,14^2,gps_c.midas.vu,'linewidth',1.25)
scatter(gps_c.measures.lons,gps_c.measures.lats,17.5^2,gps_c.measures.vu,'linewidth',1.25)
scatter(gps_c.holden15.lons,gps_c.holden15.lats,21^2,gps_c.holden15.vu,'linewidth',1.25)
scatter(gps_c.lamb13.lons,gps_c.lamb13.lats,24.5^2,gps_c.lamb13.vu,'linewidth',1.25)
scatter(gps_c.beavan12.lons,gps_c.beavan12.lats,28^2,gps_c.beavan12.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,31.5^2,gps_c.houlie17.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,35^2,gps_c.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_c.denys20.lons,gps_c.denys20.lats,38.5^2,gps_c.denys20.vu,'linewidth',1.25)
scatter(gps_c.ulr5.lons,gps_c.ulr5.lats,42^2,gps_c.ulr5.vu,'linewidth',1.25)
scatter(gps_c.beavan10.lons,gps_c.beavan10.lats,45.5^2,gps_c.beavan10.vu,'linewidth',1.25)

% quiver(gps_c.lamb13.lons,gps_c.lamb13.lats,sc*lonsc*gps_c.lamb13.ve,sc*gps_c.lamb13.vn,0,'linewidth',1.25,'maxheadsize',5e-3,'autoscale','off','color',color_lamb13);
% thisellipse = ellipse(sc*lonsc*gps_c.lamb13.se,sc*gps_c.lamb13.sn,0*gps_c.lamb13.se,gps_c.lamb13.lons + sc*lonsc*gps_c.lamb13.ve,gps_c.lamb13.lats + sc*gps_c.lamb13.vn,color_lamb13,30);
% set(thisellipse,'linewidth',0.25)

quiver(gps_c.midas.lons,gps_c.midas.lats,sc*lonsc*gps_c.midas.ve,sc*gps_c.midas.vn,0,'color',color_midas,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.midas.se,sc*gps_c.midas.sn,0*gps_c.midas.se,gps_c.midas.lons + sc*lonsc*gps_c.midas.ve,gps_c.midas.lats + sc*gps_c.midas.vn,color_midas,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.UNR.lons,gps_c.UNR.lats,sc*lonsc*gps_c.UNR.ve,sc*gps_c.UNR.vn,0,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off','color',color_unr);
thisellipse = ellipse(sc*lonsc*gps_c.UNR.se,sc*gps_c.UNR.sn,0*gps_c.UNR.se,gps_c.UNR.lons + sc*lonsc*gps_c.UNR.ve,gps_c.UNR.lats + sc*gps_c.UNR.vn,color_unr,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.measures.lons,gps_c.measures.lats,sc*lonsc*gps_c.measures.ve,sc*gps_c.measures.vn,0,'color',color_measures,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.measures.se,sc*gps_c.measures.sn,0*gps_c.measures.se,gps_c.measures.lons + sc*lonsc*gps_c.measures.ve,gps_c.measures.lats + sc*gps_c.measures.vn,color_measures,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.pickle19.lons,gps_c.pickle19.lats,sc*lonsc*gps_c.pickle19.ve,sc*gps_c.pickle19.vn,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.pickle19.se,sc*gps_c.pickle19.sn,0*gps_c.pickle19.se,gps_c.pickle19.lons + sc*lonsc*gps_c.pickle19.ve,gps_c.pickle19.lats + sc*gps_c.pickle19.vn,color_pickle19,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.beavan16.lons,gps_c.beavan16.lats,sc*lonsc*gps_c.beavan16.ve,sc*gps_c.beavan16.vn,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.beavan16.se,sc*gps_c.beavan16.sn,0*gps_c.beavan16.ve,gps_c.beavan16.lons + sc*lonsc*gps_c.beavan16.ve,gps_c.beavan16.lats + sc*gps_c.beavan16.vn,color_beavan16,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.denys16.lons,gps_c.denys16.lats,sc*lonsc*gps_c.denys16.ve,sc*gps_c.denys16.vn,0,'color',color_denys16,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.denys16.se,sc*gps_c.denys16.sn,0*gps_c.denys16.ve,gps_c.denys16.lons + sc*lonsc*gps_c.denys16.ve,gps_c.denys16.lats + sc*gps_c.denys16.vn,color_denys16,30);
set(thisellipse,'linewidth',0.25)

quiver(gps_c.page18.lons,gps_c.page18.lats,sc*lonsc*gps_c.page18.ve,sc*gps_c.page18.vn,0,'color',color_page18,'linewidth',1.25,'maxheadsize',1e-3,'autoscale','off');
thisellipse = ellipse(sc*lonsc*gps_c.page18.se,sc*gps_c.page18.sn,0*gps_c.page18.ve,gps_c.page18.lons + sc*lonsc*gps_c.page18.ve,gps_c.page18.lats + sc*gps_c.page18.vn,color_page18,30);
set(thisellipse,'linewidth',0.25)

text(177.5 - 5.5,-34.5 - 0.006 - 10 + 0.5,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
text(177.275 - 5.5,-34.68975 - 0.006 - 10 + 0.5,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.075 - 5.5,-34.625 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-34.625 - 0.006 - 10 + 0.5,'Beavan et al. [2016]','fontsize',13.5)
quiver(177.075 - 5.5,-34.75 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_pickle19,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-34.75 - 0.006 - 10 + 0.5,'Pickle [2019+]','fontsize',13.5)
quiver(177.075 - 5.5,-34.875 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_measures,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-34.875 - 0.006 - 10 + 0.5,'NASA MEASURES (2025)','fontsize',13.5)
quiver(177.075 - 5.5,-35 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_midas,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35 - 0.006 - 10 + 0.5,'UNR MIDAS (2024)','fontsize',13.5)
quiver(177.075 - 5.5,-35.125 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_unr,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35.125 - 0.006 - 10 + 0.5,'UNR [Kreemer et al., 2014]','fontsize',13.5)
quiver(177.075 - 5.5,-35.25 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_denys16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35.25 - 0.006 - 10 + 0.5,'Denys et al. [2016]','fontsize',13.5)
quiver(177.075 - 5.5,-35.375 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_page18,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(177.5 - 5.5,-35.375 - 0.006 - 10 + 0.5,'Page et al. [2018]','fontsize',13.5)

text(177.5 - 5.5,-34.75 - 0.006 - 10.8125 + 0.5,'Verticals:','fontsize',13.5,'HorizontalAlignment','center')

scatter(177.44 - 5.5,-34.875 - 10.8125 + 0.5,7^2,-1,'filled')
text(177.5 - 5.5,-34.875 - 0.006 - 10.8125 + 0.5,['Hamling et al. [2022]'],'fontsize',13.5)
scatter(177.426 - 5.5,-35 - 10.8125 + 0.5,10.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35 - 0.006 - 10.8125 + 0.5,['Pickle [2019+]'],'fontsize',13.5)
scatter(177.412 - 5.5,-35.125 - 10.8125 + 0.5,14^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.125 - 0.006 - 10.8125 + 0.5,['UNR MIDAS (2024)'],'fontsize',13.5)
scatter(177.398 - 5.5,-35.25 - 10.8125 + 0.5,17.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.25 - 0.006 - 10.8125 + 0.5,['NASA MEASURES (2025)'],'fontsize',13.5)
scatter(177.384 - 5.5,-35.375 - 10.8125 + 0.5,21^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.375 - 0.006 - 10.8125 + 0.5,'Holden et al. [2015]','fontsize',13.5)
scatter(177.37 - 5.5,-35.5 - 10.8125 + 0.5,24.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.5 - 0.006 - 10.8125 + 0.5,'Lamb and Smith [2013]','fontsize',13.5)
scatter(177.356 - 5.5,-35.625 - 10.8125 + 0.5,28^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.625 - 0.006 - 10.8125 + 0.5,['Beavan & Litchfield [2012]'],'fontsize',13.5)
scatter(177.342 - 5.5,-35.75 - 10.8125 + 0.5,31.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.75 - 0.006 - 10.8125 + 0.5,'Houlie & Stern [2017]','fontsize',13.5)
scatter(177.328 - 5.5,-35.875 - 10.8125 + 0.5,35^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-35.875 - 0.006 - 10.8125 + 0.5,'Houlie & Stern, unadjusted','fontsize',13.5)
scatter(177.314 - 5.5,-36 - 10.8125 + 0.5,38.5^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-36 - 0.006 - 10.8125 + 0.5,['Denys et al. [2020]'],'fontsize',13.5)
scatter(177.3 - 5.5,-36.125 - 10.8125 + 0.5,42^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-36.125 - 0.006 - 10.8125 + 0.5,'ULR5 [2012]','fontsize',13.5)
scatter(177.296 - 5.5,-36.25 - 10.8125 + 0.5,42^2,-1,'linewidth',1.25)
text(177.5 - 5.5,-36.25 - 0.006 - 10.8125 + 0.5,'Beavan et al. [2010]','fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([166 173.25])
ylim([-46.75 -43.375])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.01875;
ax_width = outerpos(3) - 0.0325;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, outliers removed'],'fontsize',25)
%return
%print(43,['NZGNSS_lowersouth_' date '_aligned_culled.png'],'-dpng');

%%
sc = 7.5e-3;
figure(44); clf; hold on; box on;
set(gcf,'position',[0 720 1066 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
plot(coastlines.X,coastlines.Y,'color',[0 0.5 0.75],'linewidth',0.5); plot(lakes.X,lakes.Y,'color',[0 0.5 0.75],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),9^2,gps_m.vu(~isnan(gps_m.vu)),'filled','markeredgecolor','none')
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),14^2,gps_m.vu(~isnan(gps_m.vu)) - gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),19.5^2,gps_m.vu(~isnan(gps_m.vu)) + gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
gpsplot = quiver(gps_m.lons,gps_m.lats,sc/cosd(-41)*gps_m.ve,sc*gps_m.vn,0,'color','k','linewidth',1.25,'maxheadsize',1e-3);
gpsellipse = ellipse(sc/cosd(-41)*gps_m.se,sc*gps_m.sn,0*gps_m.sn,gps_m.lons + sc/cosd(-41)*gps_m.ve,gps_m.lats + sc*gps_m.vn,'k',30);
set(gpsellipse,'linewidth',0.25)

text(178 - 5.5,-35.5 - 0.006 - 10 + 0.5,'Horizontals:','fontsize',13.5,'HorizontalAlignment','center')
%text(177.275 - 5.5,-34.68975 - 0.006 - 10 + 0.5,'40 mm/yr','fontsize',13.5,'HorizontalAlignment','center')
quiver(177.575 - 5.5,-35.625 - 10 + 0.5,sc*lonsc*40,sc*0,0,'color',color_beavan16,'linewidth',1.25,'maxheadsize',0.25,'autoscale','off');
text(178 - 5.5,-35.625 - 0.006 - 10 + 0.5,'40 mm/yr','fontsize',13.5)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([166 173.25])
ylim([-46.75 -43.375])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.01875;
ax_width = outerpos(3) - 0.0325;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Compiled horizontal and vertical GNSS velocities, aligned, outliers removed, averaged'],'fontsize',25)
%return
%print(44,['NZGNSS_lowersouth_' date '_sillremoved.png'],'-dpng');
return

%%
sc = 7.5e-3;
figure(51); clf; hold on; box on;
set(gcf,'position',[0 720 1308 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
c2 = load(['southpacificcoastlines.mat']);
plot(c2.coastlines(:,1),c2.coastlines(:,2),'color',[0.625 0.775 0.85],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_u.hamling22.lons,gps_u.hamling22.lats,7^2,gps_u.hamling22.vu,'filled')
scatter(gps_u.pickle19.lons,gps_u.pickle19.lats,10.5^2,gps_u.pickle19.vu,'linewidth',1.25)
scatter(gps_u.midas.lons,gps_u.midas.lats,14^2,gps_u.midas.vu,'linewidth',1.25)
scatter(gps_u.measures.lons,gps_u.measures.lats,17.5^2,gps_u.measures.vu,'linewidth',1.25)
scatter(gps_u.holden15.lons,gps_u.holden15.lats,21^2,gps_u.holden15.vu,'linewidth',1.25)
scatter(gps_u.lamb13.lons,gps_u.lamb13.lats,24.5^2,gps_u.lamb13.vu,'linewidth',1.25)
scatter(gps_u.beavan12.lons,gps_u.beavan12.lats,28^2,gps_u.beavan12.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,31.5^2,gps_u.houlie17.vu,'linewidth',1.25)
scatter(gps_u.houlie17.lons,gps_u.houlie17.lats,35^2,gps_u.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_u.denys20.lons,gps_u.denys20.lats,38.5^2,gps_u.denys20.vu,'linewidth',1.25)
scatter(gps_u.ulr5.lons,gps_u.ulr5.lats,42^2,gps_u.ulr5.vu,'linewidth',1.25)
scatter(gps_u.vardic21.lons,gps_u.vardic21.lats,45.5^2,gps_u.vardic21.vu,'linewidth',1.25)
scatter(gps_u.itrf14vlbi.lons,gps_u.itrf14vlbi.lats,49^2,gps_u.itrf14vlbi.vu,'linewidth',1.25)
scatter(gps_u.beavan10.lons,gps_u.beavan10.lats,52.5^2,gps_u.beavan10.vu,'linewidth',1.25)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([110 230])
ylim([-55 -10])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.02;
ax_width = outerpos(3) - 0.02;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Regional vertical GNSS velocities'],'fontsize',25)
%return
print(51,['NZGNSS_regional_' date '.png'],'-dpng');
%return
% 
%%
sc = 7.5e-3;
figure(53); clf; hold on; box on;
set(gcf,'position',[0 720 1308 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
c2 = load(['southpacificcoastlines.mat']);
plot(c2.coastlines(:,1),c2.coastlines(:,2),'color',[0.625 0.775 0.85],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_c.hamling22.lons,gps_c.hamling22.lats,7^2,gps_c.hamling22.vu,'filled')
scatter(gps_c.pickle19.lons,gps_c.pickle19.lats,10.5^2,gps_c.pickle19.vu,'linewidth',1.25)
scatter(gps_c.midas.lons,gps_c.midas.lats,14^2,gps_c.midas.vu,'linewidth',1.25)
scatter(gps_c.measures.lons,gps_c.measures.lats,17.5^2,gps_c.measures.vu,'linewidth',1.25)
scatter(gps_c.holden15.lons,gps_c.holden15.lats,21^2,gps_c.holden15.vu,'linewidth',1.25)
scatter(gps_c.lamb13.lons,gps_c.lamb13.lats,24.5^2,gps_c.lamb13.vu,'linewidth',1.25)
scatter(gps_c.beavan12.lons,gps_c.beavan12.lats,28^2,gps_c.beavan12.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,31.5^2,gps_c.houlie17.vu,'linewidth',1.25)
scatter(gps_c.houlie17.lons,gps_c.houlie17.lats,35^2,gps_c.houlie17_unadj.vu,'linewidth',1.25)
scatter(gps_c.denys20.lons,gps_c.denys20.lats,38.5^2,gps_c.denys20.vu,'linewidth',1.25)
scatter(gps_c.ulr5.lons,gps_c.ulr5.lats,42^2,gps_c.ulr5.vu,'linewidth',1.25)
scatter(gps_c.vardic21.lons,gps_c.vardic21.lats,45.5^2,gps_c.vardic21.vu,'linewidth',1.25)
scatter(gps_c.itrf14vlbi.lons,gps_c.itrf14vlbi.lats,49^2,gps_c.itrf14vlbi.vu,'linewidth',1.25)
scatter(gps_c.beavan10.lons,gps_c.beavan10.lats,52.5^2,gps_c.beavan10.vu,'linewidth',1.25)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([110 230])
ylim([-55 -10])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.02;
ax_width = outerpos(3) - 0.02;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Regional vertical GNSS velocities, aligned, outliers removed using horizontals'],'fontsize',25)
%return
print(53,['NZGNSS_regional_' date '_aligned_culled.png'],'-dpng');
%return

%%
sc = 7.5e-3;
figure(54); clf; hold on; box on;
set(gcf,'position',[0 720 1308 720])
%set(gcf,'position',[0 1025 1450 1025])
set(gca,'DataAspectRatio',[1*lonsc 1 111.1/5])
c2 = load(['southpacificcoastlines.mat']);
plot(c2.coastlines(:,1),c2.coastlines(:,2),'color',[0.625 0.775 0.85],'linewidth',0.5)
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.9 0.8 0.7],'linewidth',0.25)
%plot([174.75 176 176],[-37.85 -37.85 -37],'--','color','k')

scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),9^2,gps_m.vu(~isnan(gps_m.vu)),'filled','markeredgecolor','none')
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),14^2,gps_m.vu(~isnan(gps_m.vu)) - gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)
scatter(gps_m.lons(~isnan(gps_m.vu)),gps_m.lats(~isnan(gps_m.vu)),19.5^2,gps_m.vu(~isnan(gps_m.vu)) + gps_m.su(~isnan(gps_m.vu)),'linewidth',1.25)

colormap(jet6_1000);
clim([-2 2])
hcb = colorbar;
hcb.Location = "southoutside";
%hcb.Label.String = "Uplift (mm/yr)";
hcb.TickLabels{5} = "Uplift (mm/yr)";

xlim([110 230])
ylim([-55 -10])
fontname(gcf,'Minion Pro')
ax = gca;
ax.FontSize = 12.5;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) - 0.02;
ax_width = outerpos(3) - 0.02;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%title(['Compiled Compiled horizontal and vertical GNSS velocities'],'fontsize',20)
title(['Regional vertical GNSS velocities, aligned, outliers removed, averaged'],'fontsize',25)
%return
print(54,['NZGNSS_regional_' date '_aligned_culled_averaged.png'],'-dpng');
%return



