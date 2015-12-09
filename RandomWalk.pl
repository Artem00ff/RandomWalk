use strict;
use lib 'C:\Strawberry\perl\lib';
use ModPolymatic();
use Polymatic();
use Getopt::Long();
use Math::Trig;
# Import constants pi2, pip2, pip4 (2*pi, pi/2, pi/4).
use Math::Trig ':pi';

######################## system params ###############
my $str = 'first probe';
my $N=5;
my $M=2;
my $C=1.76;
my $b=0.97;
my %chains;
#my ($i,$k,$Q,$fi,$l);
################## make walk - chains generation ###################
sub MakeWalk(%){
#making system hash
my %chains= ('header' => $str,
            'N' => $N,
            'M' => $M,
            'Cinf' => $C,
            'b' => $b );

my ($i,$k,$Q,$fi,$l,@randpos);
my @nullpos=(0,0,0);

#calculating of critical parameters
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
return %chains;
}
############# end make walk #################
#print "Scitical params = ", rad2deg($O_max)," ", $l_max,"\n";
###################################################################################
#output in console
# %chains transmited as link
sub out_in_console{
  my %inhash = %{shift()};
  my ($i,$k);
     for ($i=1; $i<=$inhash{'M'}; $i++){
         my $tempR2 = ($inhash{'mols'}[$i][1][0]-$inhash{'mols'}[$i][$N][0])**2+($inhash{'mols'}[$i][1][1]-$inhash{'mols'}[$i][$N][1])**2+($inhash{'mols'}[$i][1][2]-$inhash{'mols'}[$i][$N][2])**2;
         print "Molecular $i R**2 = $tempR2 \n";
               for ($k=1;$k<=$inhash{'N'};$k++){
               my @temp=@{$inhash{mols}[$i][$k]};
               print "@temp \n";
               }

}
}
###################################################################################

###############   main  #############################
%chains=MakeWalk();
out_in_console(\%chains);
###############  end main ###########################

