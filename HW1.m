% Homework 1. Due before class on 9/5/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
x = 3; y = 5; % integers
x = '3'; y= '5'; %strings
x = 3; y = '5'; %mixed

%your code goes here

x_num=isnumeric(x);
y_num=isnumeric(y);

if x_num==0
   disp('X is NOT a number')
   fprintf('converting X to number \n')
   x=str2num(x);
else
    disp('X is a number')
end
if y_num==0
   disp('Y is NOT a number')
   fprintf('converting Y to number \n')
   y=str2num(y);
else
    disp('Y is a number')
end
fprintf('\n x and y are now numbers \n')
result=x+y;
%output your answer
result

%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 

N=500;
A='A'; T='T'; G='G'; C='C';
DNA=[A,T,G,C];
rand_seq=datasample(DNA,N);

%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.

N=500;
A='A'; T='T'; G='G'; C='C';
DNA=[A,T,G,C];
rand_seq=datasample(DNA,N);

start=strfind(rand_seq,'ATG');
end1=strfind(rand_seq,'TAG');
end2=strfind(rand_seq,'TGA');
end3=strfind(rand_seq,'TAA');

nstart=length(start); %columnas
nend1=length(end1); %rows
nend2=length(end2); %rows
nend3=length(end3); %rows

ORF1=bsxfun(@minus, end1, reshape(start(:),1,1,[])); %need to reshape this
rORF1=reshape(ORF1,[nend1,nstart]);

r1entre3=rORF1./3;
r1=rORF1./3;
reads1=(~mod(r1,1)); %so it can read all the substractions that are divisible by 3: this will allow me to check distances between the start codon and the stop codons at the same frame
r1(r1 < 1 | ~reads1)=0;
r1(r1<1)=NaN;
minr1=min(r1);
minr1(isnan(minr1))=[];
theORF1=max(minr1);
ORFs1=sum(length(minr1));

ORF2=bsxfun(@minus, end2, reshape(start(:),1,1,[]));
rORF2=reshape(ORF2,[nend2,nstart]);
r2entre3=rORF2./3;
r2=rORF2./3;
reads2=(~mod(r2,1));
r2(r2 < 1 | ~reads2)=0;

r2(r2<1)=NaN;
minr2=min(r2);
minr2(isnan(minr2))=[];
theORF2=max(minr2);
ORFs2=sum(length(minr2));

ORF3=bsxfun(@minus, end3, reshape(start(:),1,1,[]));
rORF3=reshape(ORF3,[nend3,nstart]);
r3entre3=rORF3./3;
r3=rORF3./3;
reads3=(~mod(r3,1));
r3(r3 < 1 | ~reads3)=0;

r3(r3<1)=NaN;
minr3=min(r3);
minr3(isnan(minr3))=[];
theORF3=max(minr3);
ORFs3=sum(length(minr3));

MAXORF=3*max(max(theORF1,theORF2),theORF3); 
fprintf('\n The sequence for the longest ORF is: %d bp', MAXORF)

%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.

N=500;
A='A'; T='T'; G='G'; C='C';
DNA=[A,T,G,C];

a=1000; 
probs=0;
pstore=zeros(1,a);

for loop=1:a;

rand_seq=datasample(DNA,N);
    
start=strfind(rand_seq,'ATG');
end1=strfind(rand_seq,'TAG');
end2=strfind(rand_seq,'TGA');
end3=strfind(rand_seq,'TAA');

nstart=length(start); %columnas
nend1=length(end1); %rows
nend2=length(end2); %rows
nend3=length(end3); %rows

ORF1=bsxfun(@minus, end1, reshape(start(:),1,1,[])); 
rORF1=reshape(ORF1,[nend1,nstart]);

r1entre3=rORF1./3;
r1=rORF1./3;
reads1=(~mod(r1,1)); %so it can read all the substractions that are divisible by 3: this will allow me to check that the start and stop codons are in the same frame
r1(r1 < 1 | ~reads1)=0;
r1(r1<1)=NaN;
minr1=min(r1);
minr1(isnan(minr1))=[];
theORF1=max(minr1);
ORFs1=sum(length(minr1));

fiftyORFs1=(minr1*3)>50;
totalfiftyORFs1=sum(sum(fiftyORFs1));


ORF2=bsxfun(@minus, end2, reshape(start(:),1,1,[]));
rORF2=reshape(ORF2,[nend2,nstart]);
r2entre3=rORF2./3;
r2=rORF2./3;
reads2=(~mod(r2,1));
r2(r2 < 1 | ~reads2)=0;

r2(r2<1)=NaN;
minr2=min(r2);
minr2(isnan(minr2))=[];
theORF2=max(minr2);
ORFs2=sum(length(minr2));

fiftyORFs2=(minr2*3)>50;
totalfiftyORFs2=sum(sum(fiftyORFs2));

ORF3=bsxfun(@minus, end3, reshape(start(:),1,1,[]));
rORF3=reshape(ORF3,[nend3,nstart]);
r3entre3=rORF3./3;
r3=rORF3./3;
reads3=(~mod(r3,1));
r3(r3 < 1 | ~reads3)=0;

r3(r3<1)=NaN;
minr3=min(r3);
minr3(isnan(minr3))=[];
theORF3=max(minr3);
ORFs3=sum(length(minr3));

fiftyORFs3=(minr3*3)>50;
totalfiftyORFs3=sum(sum(fiftyORFs3));

totalfifty=totalfiftyORFs1+totalfiftyORFs2+totalfiftyORFs3;


prob=totalfifty/N;
pstore(loop)=prob;

end

result=sum(pstore)/a;
fprintf('\n The probability of finding an ORF greater than 50 pb in a 500 bp is: \n\n ')
result

%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a function of the sequence length. 

N=500;
A='A'; T='T'; G='G'; C='C';
DNA=[A,T,G,C];

a=1000; %las veces que se repita
probs=0;
pstore=zeros(1,a);
pnewstore=zeros(1,N);

for newloop=1:N
for loop=1:a

rand_seq=datasample(DNA,N);
    
start=strfind(rand_seq,'ATG');
end1=strfind(rand_seq,'TAG');
end2=strfind(rand_seq,'TGA');
end3=strfind(rand_seq,'TAA');

nstart=length(start); %columnas
nend1=length(end1); %rows
nend2=length(end2); %rows
nend3=length(end3); %rows

ORF1=bsxfun(@minus, end1, reshape(start(:),1,1,[])); 
rORF1=reshape(ORF1,[nend1,nstart]);


r1entre3=rORF1./3;
r1=rORF1./3;
reads1=(~mod(r1,1)); %so it can read all the substractions that are divisible by 3: this will allow me to check that the start and stop codons are in the same frame
r1(r1 < 1 | ~reads1)=0;
r1(r1<1)=NaN;
minr1=min(r1);
minr1(isnan(minr1))=[];
theORF1=max(minr1);
ORFs1=sum(length(minr1));

fiftyORFs1=(minr1*3)>50;
totalfiftyORFs1=sum(sum(fiftyORFs1));


ORF2=bsxfun(@minus, end2, reshape(start(:),1,1,[]));
rORF2=reshape(ORF2,[nend2,nstart]);
r2entre3=rORF2./3;
r2=rORF2./3;
reads2=(~mod(r2,1));
r2(r2 < 1 | ~reads2)=0;

r2(r2<1)=NaN;
minr2=min(r2);
minr2(isnan(minr2))=[];
theORF2=max(minr2);
ORFs2=sum(length(minr2));

fiftyORFs2=(minr2*3)>50;
totalfiftyORFs2=sum(sum(fiftyORFs2));

ORF3=bsxfun(@minus, end3, reshape(start(:),1,1,[]));
rORF3=reshape(ORF3,[nend3,nstart]);
r3entre3=rORF3./3;
r3=rORF3./3;
reads3=(~mod(r3,1));
r3(r3 < 1 | ~reads3)=0;

r3(r3<1)=NaN;
minr3=min(r3);
minr3(isnan(minr3))=[];
theORF3=max(minr3);
ORFs3=sum(length(minr3));

fiftyORFs3=(minr3*3)>50;
totalfiftyORFs3=sum(sum(fiftyORFs3));

totalfifty=totalfiftyORFs1+totalfiftyORFs2+totalfiftyORFs3;

prob=totalfifty/N;
pstore(loop)=prob;
end
result=sum(pstore)/a;
%fprintf('\n The probability of finding an ORF greater than 50 pb in a 500 bp is: \n\n ')
result;

pnewstore(newloop)=result;
end

newresult=sum(pnewstore)/N;
newresult
plot(1:N,pnewstore)

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

%Miguel Angel: the shorter you make your sequence the less probability
%you'll get for having Opening reading frames (ORFs) greater than 50. Same
%goes for the other way around. The bigger the length the higher the
%chances are of encountering ORFs>50bps. The distance between the minimum
%and maximum values should increase greatly. 
%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 


filename='qPCRdata.txt';
delimiterIn= '\t';
headerslinesIN=2;
A=importdata(filename,delimiterIn,headerslinesIN);
B=struct2cell(A);
C=B{1,1};

rowA=C((1:12),1);
rowB=C((13:24),1);
rowC=C((25:36),1);
rowD=C((37:48),1);
rowE=C((49:60),1);
rowF=C((61:72),1);

% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 

plate=zeros(8,12);

plate(1,1:12)=rowA';
plate(2,1:12)=rowB';
plate(3,1:12)=rowC';
plate(4,1:12)=rowD';
plate(5,1:12)=rowE';
plate(6,1:12)=rowF';
plate

% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

averageplate=zeros(8,4);

averageplate(1,1)=mean(plate(1,1:3)); averageplate(2,1)=mean(plate(2,1:3));
averageplate(1,2)=mean(plate(1,4:6)); averageplate(2,2)=mean(plate(2,4:6));
averageplate(1,3)=mean(plate(1,7:9)); averageplate(2,3)=mean(plate(2,7:9));
averageplate(1,4)=mean(plate(1,10:12)); averageplate(2,4)=mean(plate(2,10:12));

averageplate(3,1)=mean(plate(3,1:3)); averageplate(4,1)=mean(plate(4,1:3));
averageplate(3,2)=mean(plate(3,4:6)); averageplate(4,2)=mean(plate(4,4:6));
averageplate(3,3)=mean(plate(3,7:9)); averageplate(4,3)=mean(plate(4,7:9));
averageplate(3,4)=mean(plate(3,10:12)); averageplate(4,4)=mean(plate(4,10:12));

averageplate(5,1)=mean(plate(5,1:3)); averageplate(6,1)=mean(plate(6,1:3));
averageplate(5,2)=mean(plate(5,4:6)); averageplate(6,2)=mean(plate(6,4:6));
averageplate(5,3)=mean(plate(5,7:9)); averageplate(6,3)=mean(plate(6,7:9));
averageplate(5,4)=mean(plate(5,10:12)); averageplate(6,4)=mean(plate(6,10:12));

averageplate;

%I just named the genes alpha, beta and gamma
genealpha=zeros(1:10);
genebeta=zeros(1:10);
genegamma=zeros(1:10);

genealpha=averageplate(1:8,1);
genebeta=averageplate(1:8,2);
genegamma=averageplate(1:8,3);
genenormal=averageplate(1:8,4);

condition1=averageplate(1,1:4);
condition2=averageplate(2,1:4);
condition3=averageplate(3,1:4);
condition4=averageplate(4,1:4);
condition5=averageplate(5,1:4);
condition6=averageplate(6,1:4);

normalizationalpha=2.^(genealpha(2:8)-genealpha(1:7)-(genenormal(2:8)-genenormal(1:7)));
normalizationbeta=2.^(genebeta(2:8)-genebeta(1:7)-(genenormal(2:8)-genenormal(1:7)));
normalizationgamma=2.^(genegamma(2:8)-genegamma(1:7)-(genenormal(2:8)-genenormal(1:7)));

normplate=zeros(6,3);
normplate(1:6,1)=normalizationalpha(1:6);
normplate(1:6,2)=normalizationbeta(1:6);
normplate(1:6,3)=normalizationgamma(1:6);
normplate;

plot(1:6,normplate)


%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

%Miguel Angel: my code from part 2 is already vectorized. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


