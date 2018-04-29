% Defining surface accesible residue positions for HOG1.
%
% Also constructing a list of all D/E/N/Q residues which are sector
% connected and surface accessible - following the ideas from Ferrell's
% evolution of phospho-regulation paper. These will be interesting
% positions to prioritize for optimizing regulation.
%
% Feb. 2017, K. Reynolds

%% Define Surface accesible positions.

% Relative solvent accessible surface areas were calculated using
% Michel Sanner's MSMS with a probe size of 1.4A, exlcuding all water and
% heteroatoms, from the structure HOG1.pdb (a swissmodel homology model). 
%(see NIH strucTools server: http://helixweb.nih.gov/structbio/basic.html)
% Import these from HOG1_a_res.area, remove ending NaN.

% Total surface areas are taken from: 
% http://www.imb-jena.de/IMAGE_AA.html, 
% they cite C.Chothia, J. Mol. Biol., 105(1975)1-14
% The areas are computed for the residue in a gly-X-gly tripeptide.
resType = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',...
        'T','V','W','Y'];
sa = [115,135,150,190,210,75,195,175,200,170,185,160,145,180,225,115,...
    140,155,255,230];

% Read in the pdb file, calculate a distance matrix, and assign residue areas.
pdb=pdbread('HOG1.pdb');
[dist,seq_pdb,ats]=pdb2dist(pdb,'A');
for k=1:numel(seq_pdb)
    resArea(k) = sa((strfind(resType,seq_pdb(k))));
end
dlmwrite('HOG1distmat.txt',dist);
fid = fopen('HOG1distats.txt', 'w');
for k=1:length(ats)
    fprintf(fid,'%s,', ats{k});
end
fclose(fid);

% Imported the surface area values from the file 2B9H_a_res.area.
% Using a cutoff of 20% RSA to define solvent exposed positions:
% Momen-Roknabadi A. BMC Bioinfo 208 9:357 
fracAccess = access./resArea';
surfPos = Resi((fracAccess > 0.2)); 
surfPosIx = find((fracAccess > 0.2));
sprintf('%i surface positions, %0.2f of the structure', numel(surfPos), numel(surfPos)/numel(Resi))
'Surface positions:'
sprintf('%g+', surfPos)

fid = fopen('HOG1_SurfPos.txt', 'w');
for k=1:length(surfPos)
    fprintf(fid,'%i ', surfPos(k));
end
fclose(fid);

%% Surface positions that are also acidic residues (D/E) 
surfDE = zeros(numel(Resi),1);
for k=1:numel(surfPosIx)
    tmpPos = surfPosIx(k);
    if ((seq_pdb(tmpPos)=='D') || (seq_pdb(tmpPos) == 'E'))%|| ...
        %(seq_pdb(tmpPos)=='N') || (seq_pdb(tmpPos) == 'Q'))
        surfDE(tmpPos) = 1;
    end
end

sprintf('%i acidic (D/E) surface positions:', sum(surfDE))
sprintf('%g+', Resi(find(surfDE)))


fid = fopen('HOG1_SurfDE.txt', 'w');
surfDEResi = Resi(find(surfDE));
for k=1:length(surfDEResi)
    fprintf(fid,'%i ', surfDEResi(k));
end
fclose(fid);

%% Scan for PKA motifs around all surface accessible D/E sites

surfDEIx = find(surfDE);
muts = zeros(numel(surfDEIx),1);
newhd = cell(numel(surfDEIx)-1,1);

%starting this loop at the second surfDE because the first is the
%N-terminal Glu6 (and thus can't have a PKA motif)
for k=2:numel(surfDEIx)
    motif(k,:) = seq_pdb(surfDEIx(k)-3:surfDEIx(k)+1);
    newmotif(k,:) = ['RR',motif(k,3),'S',motif(k,5)];
    for j=1:5
        if (motif(k,j) ~= newmotif(k,j))
            muts(k) = muts(k) + 1;
        end 
    end
    newseq(k-1,:) = seq_pdb;
    newhd{k-1} = ['Hog1 pos ', num2str(surfDEResi(k))];
    newseq(k-1,surfDEIx(k)-3:surfDEIx(k)+1) = newmotif(k,:);
    sprintf('Old motif: %s, New motif: %s, Num mutations: %i', ...
        motif(k,:), newmotif(k,:), muts(k))
end

fastawrite('HOG1_mutsites.fasta', newhd, newseq);