function mypostcallback(~,~,app) 
zdir = get(app.zm); % get the handle
if(strcmp(zdir.Direction,'out')) %is it going out?
    set(app.zm,'ActionPostCallback',[]); %disable the call back to zoom out
    zoom('out'); %zoom to original spot
    set(app.zm,'ActionPostCallback',{@mypostcallback,S}); %re-enable the callback
end
end