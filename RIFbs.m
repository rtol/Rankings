%% citations by paper
citpap = readtable("citpap.csv");
P = table2array(citpap(1:end,2:end));

NPap = size(P,1);
NJrn = 319;
%% citations by paper, no self-citations
citpapself = readtable("citpapself.csv");
Pself = table2array(citpapself(1:end,2:end));

%% citations by paper, paper sorted in blocks
citpaps = readtable("citpaps.csv");
PS = table2array(citpaps(1:end,2:end));

%% this creates and saves the matrix that maps papers to journals; about 33 minutes
Q = zeros(NJrn,NPap);

tic
for i=1:NJrn
    disp(i)
    for j=1:NPap
        Q(i,j) = strcmp(strcat(citpap.CITED_ISSN{j}(1:4),citpap.CITED_ISSN{j}(6:9)),strcat(citpap.Properties.VariableNames{i+1}(2:5),citpap.Properties.VariableNames{i+1}(7:10)));
    end
end
toc

save('concordance',"Q")

%% ditto without self-citations
Qself = zeros(NJrn,NPap);

tic
for i=1:NJrn
    disp(i)
    for j=1:NPap
        Qself(i,j) = strcmp(strcat(citpapself.Var1{j}(1:4),citpapself.Var1{j}(6:9)),strcat(citpapself.Properties.VariableNames{i+1}(2:5),citpapself.Properties.VariableNames{i+1}(7:10)));
    end
end
toc

save('concordself',"Qself")

%% ditto for clustered journals
QS = zeros(NJrn,NPap);

tic
for i=1:NJrn
    disp(i)
    for j=1:NPap
        QS(i,j) = strcmp(strcat(citpaps.CITED_ISSN{j}(1:4),citpaps.CITED_ISSN{j}(6:9)),strcat(citpaps.Properties.VariableNames{i+1}(2:5),citpaps.Properties.VariableNames{i+1}(7:10)));
    end
end
toc

save('concordsort',"QS")

%% this is much quicker
load("concordance.mat")
load("concordsort.mat")
load("concordself.mat")

%% 
papers = readtable("papers.csv"); %number of papers per journal
citations = readtable("citations.csv"); %citation matrix created in Excel to double check

C = Q*P;
CC = table2array(citations(1:end,2:end));
test = [min(min(CC-C)) max(max(CC-C))] %should be zero
A = diag(papers.papers);

RIFbest = RIF(Q,P,A);

%%
Cself = Qself*Pself;
RIFself = RIF(Qself,Pself,A);

%%
papersorted = readtable("papersorted.csv");
AS = diag(papersorted.Papers);

RIFbests = RIF(QS,PS,AS);

%% this creates the blocks for sampling within journals
From = zeros(NPap,1);
To = zeros(NPap,1);

cnt = 1;
for i=1:NJrn
    for j=1:A(i,i)
        From(cnt+j-1,1)=cnt;
        To(cnt+j-1,1)=cnt+A(i,i)-1;
    end
    cnt = cnt + A(i,i);
end
%% this creates the blocks for sampling within reordered journals
FromS = zeros(NPap,1);
ToS = zeros(NPap,1);

cnt = 1;
for i=1:NJrn
    for j=1:AS(i,i)
        FromS(cnt+j-1,1)=cnt;
        ToS(cnt+j-1,1)=cnt+AS(i,i)-1;
    end
    cnt = cnt + AS(i,i);
end

%% this creates the blocks for sampling within journal clusters
NCluster = max(papersorted.Cluster);
AC = zeros(NCluster,1);
for i=1:NCluster
    AC(i) = sum(papersorted.Papers(papersorted.Cluster==i));
end

FromC = zeros(NPap,1);
ToC = zeros(NPap,1);

cnt = 1;
for i=1:NCluster
    for j=1:AC(i)
        FromC(cnt+j-1,1)=cnt;
        ToC(cnt+j-1,1)=cnt+AC(i)-1;
    end
    cnt = cnt + AC(i);
end

%% Monte Carlo analysis, reshuffling papers within journals; 5 minutes
NBS = 1000;
U = rand([NPap NBS]);
BS = int32(From + (To-From).*U);

tic
aux = zeros(NJrn,NBS);
auc = zeros(NJrn,NBS);
for i=1:NBS
    aux(:,i) = RIF(Q,P(BS(:,i),:),A);
    auc(:,i) = RIF(Qself,Pself(BS(:,i),:),A);
end
toc

RIFmean = sum(abs(aux),2)/NBS;
RIFvar = sum(aux.*aux,2)/NBS - RIFmean.*RIFmean;
RIFstd = RIFvar.^0.5;
RIFse = RIFstd/sqrt(NBS);

RIFmeanself = sum(abs(auc),2)/NBS;
RIFvarself = sum(auc.*auc,2)/NBS - RIFmeanself.*RIFmeanself;
RIFstdself = RIFvarself.^0.5;
RIFseself = RIFstdself/sqrt(NBS);

rank = zeros(NJrn,NBS);
rankself = zeros(NJrn,NBS);
for i=1:NJrn
    rank(i,:) = 0.5 + sum(abs(aux(:,:))>abs(aux(i,:)),1) + 0.5*sum(abs(aux(:,:))==abs(aux(i,:)),1);
    rankself(i,:) = 0.5 + sum(abs(auc(:,:))>abs(auc(i,:)),1) + 0.5*sum(abs(auc(:,:))==abs(auc(i,:)),1);
end
rankup = prctile(rank',97.5);
ranklo = prctile(rank',2.5);
rankupself = prctile(rankself',97.5);
rankloself = prctile(rankself',2.5);

%% Monte Carlo analysis, reshuffling papers within and between journals
NBSS = 500; %vary this between 100 and 500
BSS = int32(FromS + (ToS-FromS).*U(:,1:NBSS));
BSC = int32(FromC + (ToC-FromC).*U(:,NBSS+1:NBS));
BSSC = [BSS BSC];

tic
auxS = zeros(NJrn,NBS);
for i=1:NBS
    auxS(:,i) = RIF(QS,PS(BSSC(:,i),:),AS);
end
toc

RIFmeanS = sum(abs(auxS),2)/NBS;
RIFvarS = sum(auxS.*auxS,2)/NBS - RIFmeanS.*RIFmeanS;
RIFstdS = RIFvarS.^0.5;
RIFseS = RIFstdS/sqrt(NBS);

rankS = zeros(NJrn,NBS);
for i=1:NJrn
    rankS(i,:) = 0.5 + sum(abs(auxS(:,:))>abs(auxS(i,:)),1) + 0.5*sum(abs(auxS(:,:))==abs(auxS(i,:)),1);
end
rankupS = prctile(rankS',97.5);
rankloS = prctile(rankS',2.5);

%% Monte Carlo analysis, reshuffling papers, no clusters
NBSS = 500; %vary this between 100 and 500
BSS = int32(FromS + (ToS-FromS).*U(:,1:NBSS));
BSC = int32(1+(NPap-1).*U(:,NBSS+1:NBS));
BSSC = [BSS BSC];

tic
auxR = zeros(NJrn,NBS);
for i=1:NBS
    auxR(:,i) = RIF(QS,PS(BSSC(:,i),:),AS);
end
toc

RIFmeanR = sum(abs(auxR),2)/NBS;
RIFvarR = sum(auxR.*auxR,2)/NBS - RIFmeanR.*RIFmeanR;
RIFstdR = RIFvarR.^0.5;
RIFseR = RIFstdR/sqrt(NBS);

rankR = zeros(NJrn,NBS);
for i=1:NJrn
    rankR(i,:) = 0.5 + sum(abs(auxR(:,:))>abs(auxR(i,:)),1) + 0.5*sum(abs(auxR(:,:))==abs(auxR(i,:)),1);
end
rankupR = prctile(rankR',97.5);
rankloR = prctile(rankR',2.5);

%% creating tables
[RIFsort, I] = sort(RIFbest,'descend');
jrn = papers.journal(I);
C00 = [RIFmean(I) ranklo(I)' rankup(I)'];

[RIFsorts, IS] = sort(RIFbests,'descend');
jrns = papersorted.name(IS);
C50 = [RIFmeanS(IS) rankloS(IS)' rankupS(IS)'];

R50 = [RIFmeanR(IS) rankloR(IS)' rankupR(IS)'];
%%
%[RIFsorts RIFmeanS(IS) RIFstdS(IS)];
scatter(RIFsorts, RIFmeanS(IS))
xlabel("Mean RIF, bootstrap within journals")
ylabel("Mean RIF, bootstrap within and between journals")
ylim([0 0.6])

%%
scatter(RIFsorts, RIFmeanR(IS))
xlabel("Mean RIF, bootstrap within journals")
ylabel("Mean RIF, bootstrap within and between journals")
ylim([0 0.6])
%%
scatter(RIFsorts, [C50(:,1)'; R50(:,1)'])
xlabel("Mean RIF, bootstrap within journals")
ylabel("Mean RIF, bootstrap within and between journals")
ylim([0 0.6])
legend({'clusters', 'no clusters'},'Location','northwest')