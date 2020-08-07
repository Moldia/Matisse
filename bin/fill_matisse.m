%Function: fill_matisse
%2020, Nilsson Lab, Scilifelab
%Objective: Combine files containing the XY position of reads in different
%samples into a single object
%Input: List containing objects containing the basic spatial information of
%genes derived from ISS, here callsed SUPEROBJECT. Every object inside the 
%main object represents sample
%Output: Reference Mastisse object
function[SPOTS] =fill_matisse(SUPEROBJECT)

SPOTS=matisseMOD

samplenames=fieldnames(SUPEROBJECT);
scales={};
images={};
NAME=[];
POS=[];
SAMPLE=[];
for i=1:size(samplenames,1)
% start processing
images{i}=SUPEROBJECT.(samplenames{i}).image;
scales{i}=SUPEROBJECT.(samplenames{i}).scale;
NAME=[NAME;SUPEROBJECT.(samplenames{i}).name];
POS=[POS;SUPEROBJECT.(samplenames{i}).pos];
SAMPLE=[SAMPLE;repelem(samplenames(i),size(SUPEROBJECT.(samplenames{i}).pos,1))'];


end

SPOTS.location=POS;
SPOTS.spotname=NAME;
SPOTS.sample=SAMPLE;
SPOTS.images=images;
SPOTS.scale=scales;

end