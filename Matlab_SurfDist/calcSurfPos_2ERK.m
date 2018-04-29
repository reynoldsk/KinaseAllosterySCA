% Defining surface accesible residue positions for 2ERK, 

% Also constructing a list of all D/E residues - following the ideas
% from Ferrell's evolution of phospho-regulation paper.
%
% April 2017, K. Reynolds

%% Define Surface accesible positions.

% Relative solvent accessible surface areas were calculated using
% Michel Sanner's MSMS with a probe size of 1.4A, exlcuding all water and
% heteroatoms, from the structure 2ERK (rat ERK2). 
%(see NIH strucTools server: http://helixweb.nih.gov/structbio/basic.html)
% Import these from 2ERK_a_res.area, remove ending NaN.

% Total surface areas are taken from: 
% http://www.imb-jena.de/IMAGE_AA.html, 
% they cite C.Chothia, J. Mol. Biol., 105(1975)1-14
% The areas are computed for the residue in a gly-X-gly tripeptide.
resType = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',...
        'T','V','W','Y'];
sa = [115,135,150,190,210,75,195,175,200,170,185,160,145,180,225,115,...
    140,155,255,230];

% Read in the pdb file, calculate a distance matrix, and assign residue areas.
pdb=pdbread('2ERK.pdb');
[dist,seq_pdb,ats]=pdb2dist(pdb,'A');
for k=1:numel(seq_pdb)
    resArea(k) = sa((strfind(resType,seq_pdb(k))));
end

% Imported the surface area values from the file 2B9H_a_res.area.
% Using a cutoff of 20% RSA to define solvent exposed positions:
% Momen-Roknabadi A. BMC Bioinfo 208 9:357 
fracAccess = access./resArea';
surfPos = Resi((fracAccess > 0.2)); 
surfPosIx = find((fracAccess > 0.2));
sprintf('%i surface positions, %0.2f of the structure', numel(surfPos), numel(surfPos)/numel(Resi))
'Surface positions:'
sprintf('%g+', surfPos)
sprintf('%g ', surfPos)
%% Identify surface accessible D/E positions 
surfDE = [];
for k=1:length(surfPos)
    if ((seq_pdb(surfPosIx(k)) == 'D') | (seq_pdb(surfPosIx(k)) == 'E'))
        surfDE = [surfDE, surfPos(k)];
    end
end

'Surface D/E residues:'
sprintf('%g+', surfDE)
sprintf('%g ', surfDE)