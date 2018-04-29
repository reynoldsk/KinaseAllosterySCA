function [dist,seq_pdb,ats,coord]=pdb2dist(pdb,chain)
% 
% Usage: [dist,seq_pdb,ats,coord]=pdb2dist(pdb,chain)
%
% This function extracts information from a PDB structure file, and is
% called within other functions to determine the contact graph of a
% reference protein.  The function takes in a PDB file (pdb) and a chain ID
% (chain), and returns: 
%
% (1) dist, the distance between each pair of residues (here calculated as
% a minimal distance between the constituent atoms), (2) seq_pdb, the
% protein sequence as given by the PDB file, (3) ats, the labels for the
% positions according to the PDB file, and (4) coord, the xyz coordinates
% for each atom (format is coord{pos}(atom,1:3)).
%
% OR 2012-01
% RR 2012-08

[coord,seq_pdb,ats]=get_coord(pdb,chain);
dist=distmat(coord);

% Coordinates of atoms, sequence, and position labels from the pdb:
function [coord,sequence,ats]=get_coord(pdb,chain)
cod3lett={'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'};
cod1lett='ARNDCEQGHILKMFPSTWYV';
indices=strmatch(chain,{pdb.Model.Atom.chainID});
pos=1; at=1; i=indices(1);
resSeq(pos)=pdb.Model.Atom(i).resSeq;
iCode{pos}=pdb.Model.Atom(i).iCode;
ats{pos}=[num2str(resSeq(pos)) iCode{pos}];
sequence=cod1lett(strcmp(cod3lett,pdb.Model.Atom(indices(i)).resName));
coord{pos}(at,1)=pdb.Model.Atom(i).X;
coord{pos}(at,2)=pdb.Model.Atom(i).Y;
coord{pos}(at,3)=pdb.Model.Atom(i).Z;
for n=2:numel(indices)
    i=indices(n);
    if ((pdb.Model.Atom(i).resSeq~=resSeq(pos)) || (~strcmp(pdb.Model.Atom(i).iCode,iCode{pos})))
        pos=pos+1; at=0;
        resSeq(pos)=pdb.Model.Atom(i).resSeq;
        iCode{pos}=pdb.Model.Atom(i).iCode;
        ats{pos}=[num2str(resSeq(pos)) iCode{pos}];
        sequence=[sequence cod1lett(strcmp(cod3lett,pdb.Model.Atom(indices(i)).resName))];
    end
    at=at+1;
    coord{pos}(at,1)=pdb.Model.Atom(i).X;
    coord{pos}(at,2)=pdb.Model.Atom(i).Y;
    coord{pos}(at,3)=pdb.Model.Atom(i).Z;   
end

% Distances between positions from the coordinates of atoms:
function dist=distmat(coord)
N_pos=numel(coord);
dist=zeros(N_pos,N_pos);
for i1=1:N_pos
    N_at1=size(coord{i1},1);
    for i2=(i1+1):N_pos
        N_at2=size(coord{i2},1);
        dmin=inf;
        for a1=1:N_at1
            for a2=1:N_at2
                d=0;
                for i=1:3
                    d=d+(coord{i1}(a1,i)-coord{i2}(a2,i)).^2;
                end
                if d<dmin, dmin=d; end
            end
        end
        dist(i1,i2)=sqrt(dmin);
    end
end
dist=dist+dist';