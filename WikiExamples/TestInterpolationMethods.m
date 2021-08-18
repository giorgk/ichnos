%% Generate sample data for 1D
p = [8,13;27,28;51,52;69,14;85,13;90,75];
%% Radial Basis functions
s = 40;
RBF_gauss =@(x) exp(-(x/s).^2);
RBF_mquad =@(x) sqrt(1+(x/s).^2);
RBF_invquad =@(x) 1./(1+(x/s).^2);
for ii = 1:size(p,1)
    for jj = 1:size(p,1)
        A(ii,jj) = RBF_gauss(pdist([p(ii,1);p(jj,1)]));
    end
end
Wrbf = A\p(:,2);
%% plot RBF
xx = -10:0.1:10;
clf
hold on
plot(xx,RBF_gauss(xx),'.-')
%plot(xx,RBF_mquad(xx),'.-')
%plot(xx,RBF_invquad(xx),'.-')
%% Inverse Distance Weighting
X = 0:1:100;
V = [];
VV = [];
pow = 3;
idst = (1:size(p,1))';
rexp = 1/(2*sqrt(6/100));
ww = (100/6)./(100*p(:,2)/sum(p(:,2)));
for x = X
    
        
    
    %idx = find(p(:,1) < x);
    %if isempty(idx)
    %    idst = 1;
    %else
    %    idst = [idx(end); idx(end) + 1];
    %    idst(idst>=7,:) = [];
    %end
        
    dst = pdist2(x, p(idst,1),'euclidean')';
    robs = sum(dst)/length(dst);% mean
    Rs = robs/rexp;
    %dst = sqrt((x - p(idst,1)).^2);
    idx = find(dst < 0.01);
    if isempty(idx)
        
        w = 1./(dst.^pow);
        %w = w.*(dst/sum(dst));
        V =[V; sum(w.*p(idst,2))/sum(w)];
    else
        V = [V; p(idst(idx),2)];
    end
    tmp = RBF_gauss(dst);
    %ww = tmp./sum(tmp);
    VV = [VV; sum(Wrbf.*tmp)];
end

clf
plot(p(:,1),p(:,2),'.-')
hold on
plot(X,V,'.')
%plot(X,VV,'.')