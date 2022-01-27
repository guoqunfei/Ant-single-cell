#!/usr/bin/perl -w
use strict;
open IN,"$ARGV[0]";
my$g=<IN>;
chomp$g;
my@title=split(/\t/,$g);
my%exp;
my%spe;
my%cell;
while(<IN>){
	chomp;
	my@A=split(/\t/);
	for(my$i=1;$i<=$#A;$i++){
		my@B=split(/_/,$title[$i-1]);		
		$exp{$A[0]}{$title[$i-1]}=$A[$i];
		$spe{$B[0]}{$title[$i-1]}++;
		$cell{$title[$i-1]}=$B[0];
	}
}
close IN;
print "$g\n";
foreach my$key(keys %exp){
	print "$key";
	my%sd;
	foreach my$spe(keys %spe){
		my@temp;
		foreach my$cell(keys %{$spe{$spe}}){
			push @temp, $exp{$key}{$cell};
		}
		my@sta=sta(@temp);
		$sd{$spe}{mean}=$sta[0];
		$sd{$spe}{sd}=$sta[1];
	}
	foreach my$ele(@title){
		my@B=split(/_/,$ele);
		if($sd{$B[0]}{sd}==0){$sd{$B[0]}{sd}=0.01}
		my$exp=($exp{$key}{$ele}-$sd{$B[0]}{mean})/$sd{$B[0]}{sd};
		print "\t$exp";
	}
	print "\n";
}
sub sta{
	my$num=0;
	my$sum=0;
	foreach my$ele(@_){
		$sum+=$ele;
		$num++;
	}
	my$mean=$sum/$num;
	my$sd=0;
	foreach my$ele(@_){
		$sd+=($ele-$mean)**2;
	}
	$sd=($sd/($num-1))**0.5;
	return($mean,$sd);
}
