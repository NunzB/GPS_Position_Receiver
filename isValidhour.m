function res = isValidhour( hour )

if (length(hour)~=8)
    res=0;
    return;
end

HHstr = hour(1:2);
HH=str2double(HHstr);

if (HH<0 || HH>23 || isnan(HH)==1)
    res=0;
    return;
end

MMstr = hour(4:5);
MM=str2double(MMstr);
if (isempty(MM)==1 || (60<MM)&& (MM<0) || isnan(MM)==1)
   res=0;
   return;
end

SSstr = hour(7:8);
SS=str2double(SSstr);
if (isempty(SS)==1 || (60<SS)&& (SS<0) || isnan(SS)==1)
   res=0;
   return;
end

if (hour(3)~=':'||hour(6)~=':')
    res=0;
    return;
end
res=1;
end

