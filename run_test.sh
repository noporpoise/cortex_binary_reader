
if [ $# -ne 3 ]
then
  echo "usage: ./run_test.sh <ctx_bin_path> <suffix> <kmer>"
  exit;
fi

CTX_PATH=$1
SUFFIX=$2
KMER=$3


mkdir -p test
FIRST=test/seq1.k$KMER.$SUFFIX.ctx
SECOND=test/seq2.k$KMER.$SUFFIX.ctx
JOINT=test/joint.k$KMER.c2.$SUFFIX.ctx

CMD1=$CTX_PATH/cortex_var_31_c1
CMD2=$CTX_PATH/cortex_var_31_c2

if [ $KMER -gt 31 ]
then
  CMD1=$CTX_PATH/cortex_var_63_c1
  CMD2=$CTX_PATH/cortex_var_63_c2
fi

if [ ! -e $CMD1 ]
then
  echo "Cannot find $CMD1"
  exit -1
fi

if [ ! -e $CMD2 ]
then
  echo "Cannot find $CMD2"
  exit -1
fi

rm -rf $FIRST $SECOND $JOINT

$CMD1 --kmer_size $KMER --se_list data/seq1.falist --dump_binary $FIRST &> $FIRST.log
$CMD1 --kmer_size $KMER --se_list data/seq2.falist --dump_binary $SECOND &> $SECOND.log
$CMD2 --kmer_size $KMER --colour_list data/joint.colours --dump_binary $JOINT &> $JOINT.log

./cortex_bin_reader $FIRST
echo
echo "=========================="
echo
./cortex_bin_reader $SECOND
echo
echo "=========================="
echo
./cortex_bin_reader $JOINT
