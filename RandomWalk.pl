use strict;
use lib 'C:\Strawberry\perl\lib';
use ModPolymatic();
use Polymatic();
use Getopt::Long();
use Math::Trig;
# Import constants pi2, pip2, pip4 (2*pi, pi/2, pi/4).
use Math::Trig ':pi';


my $str = 'first probe';
my $N=5;
my $M=2;
my $C=1.76;
my $b=0.97;
my ($i,$k);
my @randpos;
my $Q;
my $fi;
my $l;
my @nullpos=(0,0,0);

my %chains= ('header' => $str,
            'N' => $N,
            'M' => $M,
            'Cinf' => $C,
            'b' => $b );
my $O_max=2*acos(sqrt(($chains{'Cinf'}-1)/($chains{'Cinf'}+1)));
my $l_max=$chains{'b'}*sqrt((1-cos(pi-$O_max))**2+sin($O_max)**2);


#make first atom
for ($i=1; $i<=$chains{'M'}; $i++) {
@{$chains{'mols'}[$i][1]} = @nullpos;
}
#make second atom
for ($i=1; $i<=$chains{'M'}; $i++) {
$fi = 2*pi*rand();
$Q = pi*rand();
@randpos = ($chains{'b'}*sin($Q)*cos($fi),$chains{'b'}*sin($Q)*sin($fi),$chains{'b'}*cos($Q));
@{$chains{'mols'}[$i][2]} = @randpos;
}




#make atoms from 3 to N
for ($k=3; $k<=$N; $k++){
for ($i=1; $i<=$M; $i++) {
my $fuse=0;
do{
$fuse++;

$fi = 2*pi*rand();
$Q = pi*rand();
#print "rand $fi $Q\n";
#print "$chains{'mols'}[$i][2][0] $chains{'mols'}[$i][2][1] $chains{'mols'}[$i][2][2] \n";
$randpos[0] = $chains{'mols'}[$i][$k-1][0] + $chains{'b'}*sin($Q)*cos($fi);
$randpos[1] = $chains{'mols'}[$i][$k-1][1] + $chains{'b'}*sin($Q)*sin($fi);
$randpos[2] = $chains{'mols'}[$i][$k-1][2] + $chains{'b'}*cos($Q);
#print  "$randpos[0] $randpos[1] $randpos[2] \n";
$l=sqrt(($chains{'mols'}[$i][$k-2][0]-$randpos[0])**2+($chains{'mols'}[$i][$k-2][1]-$randpos[1])**2+($chains{'mols'}[$i][$k-2][2]-$randpos[2])**2);
print $k," ",$i," ",$fuse," ",$l," ", $l_max,"\n";
if ($fuse>1000) { print "fuse is over N= $k M= $i"; last;}
if ($l>1.94) {print "ERROR";}
} until ($l>$l_max) ;
#print "$i $fi $Q \n";
@{$chains{'mols'}[$i][$k]} = @randpos;
}
}
#print "Scitical params = ", rad2deg($O_max)," ", $l_max,"\n";
out_in_console(%chains);
#
#output in console
sub out_in_console{
for ($i=1; $i<=$_{'M'}; $i++){
 my $tempR2 = ($_{'mols'}[$i][1][0]-$_{'mols'}[$i][$N][0])**2+($_{'mols'}[$i][1][1]-$_{'mols'}[$i][$N][1])**2+($_{'mols'}[$i][1][2]-$_{'mols'}[$i][$N][2])**2;
 print "Molecular $i R**2 = $tempR2 \n";
    for ($k=1;$k<=$_{'N'};$k++){
        my @temp=@{$_{mols}[$i][$k]};
        print "@temp \n";
        }
}
}