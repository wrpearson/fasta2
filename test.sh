#!/bin/csh -f
echo "starting fasta - protein"
fasta -q mgstm1.aa q > test_m1.ok2
echo "done"
echo "starting fastx"
fastx -q mgstm1.e05 q > test_m1.xk2
echo "done"
echo "starting fastx (reverse)"
fastx -q -i mgstm1.rev q > test_m1.xk2r
echo "done"
echo "starting fasta - DNA"
fasta -q mgstm1.seq %PRMB > test_m1.ok4
echo "done"
echo "starting ssearch"
ssearch -q mgstm1.aa q > test_m1.ss
echo "done"
echo "starting tfasta"
tfasta -q mgstm1.aa %PRMB > test_m1.tk2
echo "done"
echo "starting tfastx"
tfastx -q mgstm1.aa %PRMB > test_m1.tx2
echo "done"
