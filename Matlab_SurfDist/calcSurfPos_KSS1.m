% Defining surface accesible residue positions for KSS1, and
% calculating which of these are sector or sector-connected positions.
%
% Also constructing a list of all D/E/N/Q residues which are sector
% connected and surface accessible - following the ideas from Ferrell's
% evolution of phospho-regulation paper. These will be interesting
% positions to prioritize for optimizing regulation.
%
% May 2015, K. Reynolds

%% Define Surface accesible positions.

% Relative solvent accessible surface areas were calculated using
% Michel Sanner's MSMS with a probe size of 1.4A, exlcuding all water and
% heteroatoms, from the structure KSS1 (a homology model). 
%(see NIH strucTools server: http://helixweb.nih.gov/structbio/basic.html)
% Import these from KSS1_a_res.area, remove ending NaN.

% Total surface areas are taken from: 
% http://www.imb-jena.de/IMAGE_AA.html, 
% they cite C.Chothia, J. Mol. Biol., 105(1975)1-14
% The areas are computed for the residue in a gly-X-gly tripeptide.
resType = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S',...
        'T','V','W','Y'];
sa = [115,135,150,190,210,75,195,175,200,170,185,160,145,180,225,115,...
    140,155,255,230];

% Read in the pdb file, calculate a distance matrix, and assign residue areas.
pdb=pdbread('KSS1.pdb');
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

%% Compare to sector and sector-connected positions

% The list of sector positions as defined in:
% ~/Documents/Collaborations/ResnekovPincus_PTM/Kinases/1501_KR_Kinases/FullKinAln.ipynb
cutoffs = [0.95];
sec{1} = [13,20,22,23,24,25,27,30,39,40,41,42,56,58,59,63,69,71,73,75,78,80,...
    91,92,93,94,95,97,99,100,115,128,130,131,133,134,135,139,141,142,143,...
    145,146,148,149,150,158,161,162,163,164,165,185,187,188,194,195,203,...
    208,210,211,213,214,218,219,284,287,288];

%surface and sector positions:
for k=1:numel(sec)
    sec_surf = intersect(sec{k}, surfPos);
    sprintf('Surface and sector positions:')
    sprintf('%g+', sec_surf)
    sprintf('Sector contacting surface positions:')
    sec_contact = zeros(numel(surfPos));
    for l=1:numel(surfPos)
        contacts = Resi(dist(surfPosIx(l),:) < 4.0);
        if ~isempty(intersect(contacts,sec{k}))
             sec_contact(l) = 1;
        end
    end
    secContSurf = surfPos(find(sec_contact)); 
    sprintf('%g+', secContSurf)
    sprintf('%i sec-contacting pos, %.2f of surface', numel(secContSurf), numel(secContSurf)/numel(surfPos))
end
    

%% All sector contacting positions (without regard to solvent accessibility)

for k=1:numel(sec)
    sprintf('Sector contacting positions - cutoff %.2f:', cutoffs(k))
    allsec_contact = zeros(numel(Resi),1);
    for l=1:numel(Resi)
        contacts = Resi(dist(l,:) < 4.0);
        if ~isempty(intersect(contacts,sec{k}))
             allsec_contact(l) = 1;
        end
    end
    secContAll = Resi(find(allsec_contact)); 
    secNotContAll = Resi(find(~(allsec_contact)));
    sprintf('%i sec-contacting pos, %.2f of protein', numel(secContAll), numel(secContAll)/numel(Resi))
    sprintf('%g,', secContAll)
    sprintf('%i non-sec-contacting pos, %.2f of protein', numel(secNotContAll), numel(secNotContAll)/numel(Resi))
    sprintf('%g,', secNotContAll)
end

%% Surface positions that are also acidic residues (D/E) 
acidicConn = zeros(numel(Resi),1);
acidicNonConn = zeros(numel(Resi),1);
for k=1:numel(surfPosIx)
    tmpPos = surfPosIx(k);
    if ((seq_pdb(tmpPos)=='D') || (seq_pdb(tmpPos) == 'E'))%|| ...
        %(seq_pdb(tmpPos)=='N') || (seq_pdb(tmpPos) == 'Q'))
        contacts = Resi(dist(tmpPos,:) < 4.0);
        if ~isempty(intersect(contacts,sec{1}))
            acidicConn(tmpPos) = 1;
        else
            acidicNonConn(tmpPos) = 1;
        end
    end
end

sprintf('%i acidic, sector connected, surface positions:', sum(acidicConn))
sprintf('%g+', Resi(find(acidicConn)))
sprintf('%i acidic, NOT connected, surface positions:', sum(acidicNonConn))
sprintf('%g+', Resi(find(acidicNonConn)))

%% Scan for PKA motifs around the acidic positions of interest

acConnIx = find(acidicConn);
muts = zeros(numel(acConnIx),1);
newhd = cell(numel(acConnIx),1);
for k=1:numel(acConnIx)
    motif(k,:) = seq_pdb(acConnIx(k)-3:acConnIx(k)+1);
    newmotif(k,:) = ['RR',motif(k,3),'S',motif(k,5)];
    for j=1:5
        if (motif(k,j) ~= newmotif(k,j))
            muts(k) = muts(k) + 1;
        end 
    end
    newseq(k,:) = seq_pdb;
    newhd{k} = ['Kss1 pos ', num2str(Resi(acConnIx(k)))];
    newseq(k,acConnIx(k)-3:acConnIx(k)+1) = newmotif(k,:);
    sprintf('Old motif: %s, New motif: %s, Num mutations: %i', ...
        motif(k,:), newmotif(k,:), muts(k))
end

fastawrite('KSS1_mutsites.fasta', newhd, newseq);