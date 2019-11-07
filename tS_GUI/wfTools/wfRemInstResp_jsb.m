function [corrTraces] = wfRemInstResp_jsb(Traces, new_response)
% This removes the instrument response from a bunch of traces.
% First it gets the responses if needed, then corrects (one at a time), and
% returns a traces array of structures of the same size as the input, but
% with the instrument response removed. It also adds a comment that the
% response was removed on such and scuh date.
% USAGE: [corrTraces] = wfRemInstResp(Traces)
% where both corrTraces and Traces are arrays of trace structures as would
% be returned by irisFetch.traces.

corrTraces=Traces;
corrComment=['instrument response removed using wfInstResp by ' getenv('username') ' on ' datestr(now)];
for k=1:length(Traces)
    
    %note: todo, here I could add something to check if sacpz is populated.
    %If it is, I can just read that information into resp and not use
    %wfGetResp
    
    resp=wfGetResp(Traces(k));
    corrTraces(k) = wfInstCorr1trace(Traces(k), resp);
    

end


%once all is said and done, add the comments.
%note, this is not all done in the same for loop to avoid problems if the
%original traces structures don't already have a 'comments' field.

for k=1:length(Traces)
    if isfield(corrTraces,'comments')
        C=corrTraces(k).comments;
        C=[C;corrComment];
        corrTraces(k).comments=C;
    else
        corrTraces(k).comments=corrComment;
    end
end