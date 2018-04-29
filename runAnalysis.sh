echo "CMGC alignment, kpos=4, p=0.95 Calculations:" > Outputs/CMGC_095_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaProcessMSA.py Inputs/CMGC_KinaseMSA.fasta -s 2PK9 -c A --output CMGC_095_2PK9 > Outputs/CMGC_095_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaCore.py Outputs/CMGC_095_2PK9.db >> Outputs/CMGC_095_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaSectorID.py -p 0.95 -k 4 Outputs/CMGC_095_2PK9.db >> Outputs/CMGC_095_2PK9.log
echo "CMGC alignment, kpos=4, p=0.96 Calculations:" > Outputs/CMGC_096_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaProcessMSA.py Inputs/CMGC_KinaseMSA.fasta -s 2PK9 -c A --output CMGC_096_2PK9 > Outputs/CMGC_096_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaCore.py Outputs/CMGC_096_2PK9.db >> Outputs/CMGC_096_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaSectorID.py -p 0.96 -k 4 Outputs/CMGC_096_2PK9.db >> Outputs/CMGC_096_2PK9.log
echo "CMGC alignment, kpos=4, p=0.95 Calculations:" > Outputs/CMGC_097_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaProcessMSA.py Inputs/CMGC_KinaseMSA.fasta -s 2PK9 -c A --output CMGC_097_2PK9 > Outputs/CMGC_097_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaCore.py Outputs/CMGC_097_2PK9.db >> Outputs/CMGC_097_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaSectorID.py -p 0.97 -k 4 Outputs/CMGC_097_2PK9.db >> Outputs/CMGC_097_2PK9.log
echo "CMGC alignment, kpos=4, p=0.95 Calculations:" > Outputs/CMGC_098_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaProcessMSA.py Inputs/CMGC_KinaseMSA.fasta -s 2PK9 -c A --output CMGC_098_2PK9 > Outputs/CMGC_098_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaCore.py Outputs/CMGC_098_2PK9.db >> Outputs/CMGC_098_2PK9.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaSectorID.py -p 0.98 -k 4 Outputs/CMGC_098_2PK9.db >> Outputs/CMGC_098_2PK9.log


echo "Alignment truncated to 2PK9, kpos=4, p=0.95 Calculations:" > Outputs/masterAln_095_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaProcessMSA.py Inputs/masterAln.an -s 2PK9 -c A -t --output masterAln_095_2PK9t > Outputs/masterAln_095_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaCore.py Outputs/masterAln_095_2PK9t.db >> Outputs/masterAln_095_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaSectorID.py Outputs/masterAln_095_2PK9t.db -k 4 -p 0.95 >> Outputs/masterAln_095_2PK9t.log
echo "Alignment truncated to 2PK9, kpos=4, p=0.96 Calculations:" > Outputs/masterAln_096_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaProcessMSA.py Inputs/masterAln.an -s 2PK9 -c A -t --output masterAln_096_2PK9t >> Outputs/masterAln_096_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaCore.py Outputs/masterAln_096_2PK9t.db >> Outputs/masterAln_096_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaSectorID.py Outputs/masterAln_096_2PK9t.db -k 4 -p 0.96 >> Outputs/masterAln_096_2PK9t.log
echo "Alignment truncated to 2PK9, kpos=4, p=0.97 Calculations:" > Outputs/masterAln_097_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaProcessMSA.py Inputs/masterAln.an -s 2PK9 -c A -t --output masterAln_097_2PK9t >> Outputs/masterAln_097_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaCore.py Outputs/masterAln_097_2PK9t.db >> Outputs/masterAln_097_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaSectorID.py Outputs/masterAln_097_2PK9t.db -k 4 -p 0.97 >> Outputs/masterAln_097_2PK9t.log
echo "Alignment truncated to 2PK9, kpos=4, p=0.98 Calculations:" > Outputs/masterAln_098_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaProcessMSA.py Inputs/masterAln.an -s 2PK9 -c A -t --output masterAln_098_2PK9t >> Outputs/masterAln_098_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaCore.py Outputs/masterAln_098_2PK9t.db >> Outputs/masterAln_098_2PK9t.log
~/Documents/PySCA/pysca-toolbox-dist6_3/scaSectorID.py Outputs/masterAln_098_2PK9t.db -k 4 -p 0.98 >> Outputs/masterAln_098_2PK9t.log
