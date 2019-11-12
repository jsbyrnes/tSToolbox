%make the kernel thingie
hSize=7;
S=fspecial('gaussian',hSize,1.2); %S is the kernel
S=S/max(S(:));

mid=round(hSize/2);
% S(mid,mid)=0;
% S=S/(sum(S(:)));
% 
% S(mid,mid)=-1;

%define the matrix
SM=zeros(size(xMat));
nCols=size(SM,2);
nRows=size(SM,1);

S0=S; %this is the untrimmed version
SM0=SM;

Gsmth=nan(nModel);

gri=0; %initialize counter for row of Gsmth

for cix=1:size(xMat,2)
    for rix=1:size(xMat,1)
        
        S=S0; % start with a fresh, untrimmed S each time.
        SM=SM0; %start with a blank slate as well.
        
        %left, right, above, below
        width=mid-1; %how many elements to each side of the central point
        lft=cix-width;
        rgt=cix+width;
        abv=rix-width;
        blw=rix+width;
        
        %to handle edges and corners:
        lft=max(lft,1); rgt=min(rgt,size(SM,2));
        abv=max(abv,1); blw=min(blw,size(SM,1));
        
        %trim S accordingly:
        %trim left side if necessary
        margin=mid-cix;
        if margin>0; S=S(:,1+margin:end); end
        %trim right side if necessary
        margin=cix+width-nCols;
        if margin>0; S=S(:,1:end-margin); end
        %trim top if necessary
        margin=mid-rix;
        if margin>0; S=S(1+margin:end,:); end
        %trim bottom if necessary
        margin=rix+width-nRows;
        if margin>0; S=S(1:end-margin,:); end
        
        %place the thing
        try
        SM(abv:blw,lft:rgt)=S;
        catch
            keyboard
        end
        
        %normalize so the whole thing sums up to zero.
        Mtemp=SM;
        ix=Mtemp==1;
        Mtemp(ix)=0;
        Msum=sum(Mtemp(:));
        Mtemp=Mtemp/Msum;
        Mtemp(ix)=-1;
        SM=Mtemp;
        
        
%         figure(1); clf; imagesc(SM);
%         colormap bone
%         caxis([-.15 .15]);
%         colorbar

        SMvec=SM(:);
        gri=gri+1; %smoothness matrix row index. note inner loop is on rows to make this consistent.
        Gsmth(gri,:)=SMvec';
        
    end
end
